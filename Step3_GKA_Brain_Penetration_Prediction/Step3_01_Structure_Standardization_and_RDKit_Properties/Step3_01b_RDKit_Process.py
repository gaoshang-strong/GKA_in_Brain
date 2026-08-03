#!/usr/bin/env python
"""
Step3_01b：RDKit 结构标准化与描述符计算

读 Step3_01_Molecule_Input.csv（1,274 个分子 = 787 GKA 候选 + 487 入脑对照），
逐个分子做标准化并算描述符，写出 Step3_01_RDKit_Processed.csv。

标准化流程（顺序固定，每一步都记录改动）：
  1. Cleanup          去除不规范价态、消毒芳香性、断开金属配位
  2. LargestFragment  去盐/去反离子，保留最大有机片段
  3. Uncharger        中性化
  4. canonical SMILES

⚠ 第 3 步不是可选的。SwissADME FAQ 原文：
   "Is it preferable to input the neutral form of the molecule? Yes, definitively.
    The SMILES entry is taken as given and not neutralized"
   模型基本在中性化合物上训练，喂离子型结构会引入显著预测偏差。

⚠ 不做删除式去冗余：787 个 GKA 候选一个不删，去盐后出现相同结构的只打标
   （dup_group / is_dup_representative），行全部保留。骨架也只作分组标注。

⚠ CNS MPO 的 6 个参数这里只能算出 4 个（MW / cLogP / TPSA / HBD）。
   logD7.4 与 pKa RDKit 算不出来，留给 ADMETlab（见 Step3_03）。

用法：
  micromamba run -n GKA_in_Brain python Step3_01b_RDKit_Process.py
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / "Step3_01_Molecule_Input.csv"

# 含羧酸的分子在生理 pH 带负电，用 cLogP 代 logD 会高估膜通透性。
# 这批分子里羧酸占比不低，必须逐个标出来，供下游区别对待。
SMARTS_FLAGS = {
    "has_carboxylic_acid": "[CX3](=O)[OX2H1]",
    "has_basic_amine": "[NX3;H2,H1,H0;!$(N[!#6]);!$(N=*);!$(N#*);!$(NC=O);!$(Nc)]",
    "has_sulfonamide": "[SX4](=O)(=O)[NX3]",
    "has_tetrazole": "c1nnn[nH]1",
}


# ---------------------------------------------------------------- 标准化


def build_standardizer():
    """构造标准化器。复用同一组对象，避免逐分子重复初始化。"""
    return (rdMolStandardize.LargestFragmentChooser(),
            rdMolStandardize.Uncharger())


def standardize(smiles: str, lfc, unch) -> dict:
    """返回标准化结果与逐步改动记录。任何失败都如实返回，不抛出。"""
    res = {
        "parse_ok": False, "std_ok": False, "std_smiles": None,
        "n_fragments_in": np.nan, "was_desalted": False,
        "charge_before": np.nan, "charge_after": np.nan,
        "was_neutralized": False, "std_changed_structure": False,
        "std_note": "",
    }
    if not isinstance(smiles, str) or not smiles.strip():
        res["std_note"] = "输入 SMILES 为空"
        return res

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        res["std_note"] = "RDKit 无法解析 SMILES"
        return res
    res["parse_ok"] = True
    res["n_fragments_in"] = len(Chem.GetMolFrags(mol))
    res["charge_before"] = Chem.GetFormalCharge(mol)

    notes = []
    try:
        m = rdMolStandardize.Cleanup(mol)
        frags_after_cleanup = len(Chem.GetMolFrags(m))
        m2 = lfc.choose(m)
        if len(Chem.GetMolFrags(m)) > 1:
            res["was_desalted"] = True
            notes.append(f"去盐：{frags_after_cleanup} 个片段取最大")
        m3 = unch.uncharge(m2)
        res["charge_after"] = Chem.GetFormalCharge(m3)
        if res["charge_after"] != Chem.GetFormalCharge(m2):
            res["was_neutralized"] = True
            notes.append(f"中性化：{Chem.GetFormalCharge(m2)} → {res['charge_after']}")
        smi = Chem.MolToSmiles(m3)
        res["std_ok"] = True
        res["std_smiles"] = smi
        # 与原输入的规范化形式比，判断结构是否真的变了
        res["std_changed_structure"] = smi != Chem.MolToSmiles(mol)
        # Cleanup 单独造成的归一化（价态/芳香性/互变写法）不在上面两项里，
        # 不单独记一笔的话 std_note 会与 std_changed_structure 自相矛盾
        if res["std_changed_structure"] and not notes:
            notes.append("Cleanup 归一化（未去盐、未改电荷）")
        res["std_note"] = "；".join(notes) if notes else "无改动"
    except Exception as exc:                                  # noqa: BLE001
        res["std_note"] = f"标准化异常：{type(exc).__name__}: {exc}"
    return res


# ---------------------------------------------------------------- 描述符


def descriptors(smiles: str | None, patterns: dict) -> dict:
    """标准化后结构的 RDKit 描述符。CNS MPO 相关的单独标注来源。"""
    keys = ["mw", "mw_exact", "heavy_atoms", "tpsa", "tpsa_sandp", "clogp",
            "molmr", "hbd", "hba", "rtb", "rtb_nonstrict",
            "rings", "aromatic_rings", "fsp3",
            "heteroatoms", "formal_charge", "qed", "lipinski_violations",
            "murcko_scaffold", "murcko_generic", "is_acyclic",
            "inchi", "inchikey", "inchikey14", "n_stereocenters"]
    out = {k: None for k in keys}
    out.update({k: False for k in patterns})
    if not smiles:
        return out
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return out

    out["mw"] = Descriptors.MolWt(m)
    out["mw_exact"] = Descriptors.ExactMolWt(m)
    out["heavy_atoms"] = m.GetNumHeavyAtoms()
    # 两种口径都留：RDKit 默认不计 S/P 的极性贡献，SwissADME/BOILED-Egg 计入。
    # 实测 193/199 与 SwissADME 对齐的是 includeSandP=True 那个，选错会改 CNS MPO 分数。
    out["tpsa"] = Descriptors.TPSA(m)
    out["tpsa_sandp"] = Descriptors.TPSA(m, includeSandP=True)
    out["clogp"] = Crippen.MolLogP(m)
    out["molmr"] = Crippen.MolMR(m)
    out["hbd"] = Descriptors.NumHDonors(m)
    out["hba"] = Descriptors.NumHAcceptors(m)
    # 同理：SwissADME 用非严格定义，实测 199/199 与 NonStrict 完全一致
    out["rtb"] = Descriptors.NumRotatableBonds(m)
    out["rtb_nonstrict"] = rdMolDescriptors.CalcNumRotatableBonds(
        m, rdMolDescriptors.NumRotatableBondsOptions.NonStrict)
    out["rings"] = rdMolDescriptors.CalcNumRings(m)
    out["aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(m)
    out["fsp3"] = rdMolDescriptors.CalcFractionCSP3(m)
    out["heteroatoms"] = rdMolDescriptors.CalcNumHeteroatoms(m)
    out["formal_charge"] = Chem.GetFormalCharge(m)
    try:
        out["qed"] = QED.qed(m)
    except Exception:                                          # noqa: BLE001
        out["qed"] = None
    out["lipinski_violations"] = int(
        (out["mw"] > 500) + (out["clogp"] > 5)
        + (out["hbd"] > 5) + (out["hba"] > 10)
    )

    scaf = MurckoScaffold.GetScaffoldForMol(m)
    out["murcko_scaffold"] = Chem.MolToSmiles(scaf)
    try:
        out["murcko_generic"] = Chem.MolToSmiles(
            MurckoScaffold.MakeScaffoldGeneric(scaf))
    except Exception:                                          # noqa: BLE001
        out["murcko_generic"] = None
    out["is_acyclic"] = out["murcko_scaffold"] == ""

    try:
        out["inchi"] = Chem.MolToInchi(m)
        out["inchikey"] = Chem.MolToInchiKey(m)
        out["inchikey14"] = out["inchikey"][:14] if out["inchikey"] else None
    except Exception:                                          # noqa: BLE001
        pass
    out["n_stereocenters"] = len(Chem.FindMolChiralCenters(
        m, includeUnassigned=True, useLegacyImplementation=False))

    for name, sma in patterns.items():
        out[name] = m.HasSubstructMatch(sma)
    return out


# ---------------------------------------------------------------- 自检


def selfcheck(df: pd.DataFrame) -> list[str]:
    p = []
    n_fail = int((~df.parse_ok).sum())
    if n_fail:
        p.append(f"{n_fail} 个分子 SMILES 无法解析（已保留并标注，未丢行）")
    n_std_fail = int((df.parse_ok & ~df.std_ok).sum())
    if n_std_fail:
        p.append(f"{n_std_fail} 个分子解析成功但标准化失败")

    # 输入行数必须原样保留——本步不允许删行
    if df.mol_id.duplicated().any():
        p.append("mol_id 出现重复，编号逻辑有问题")

    # 中性化后仍带电的分子：SwissADME 要求中性输入，这些要单独知会下游
    charged = df[df.std_ok & (df.formal_charge != 0)]
    if len(charged):
        p.append(f"{len(charged)} 个分子中性化后仍带净电荷（季铵盐等无法中性化），"
                 f"例：{', '.join(str(x) for x in charged.compound_name.head(5))}")

    # 阳性对照必须全部走完标准化
    pc = df[df.gka_is_positive_control.fillna(False).astype(bool)]
    if len(pc) and not pc.std_ok.all():
        p.append(f"有 GKA 阳性对照未通过标准化：{list(pc[~pc.std_ok].source_id)}")
    return p


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_IN)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    started = datetime.now()
    if not args.input.exists():
        sys.exit(f"[FATAL] 找不到输入：{args.input}（先跑 Step3_01a）")

    print("[1/5] 读入输入表 …")
    df = pd.read_csv(args.input)
    print(f"      {len(df)} 个分子：{df['set'].value_counts().to_dict()}")

    print("[2/5] 标准化（Cleanup → 去盐 → 中性化）…")
    lfc, unch = build_standardizer()
    std = pd.DataFrame([standardize(s, lfc, unch) for s in df.input_smiles])
    print(f"      解析成功 {std.parse_ok.sum()}/{len(std)}；"
          f"标准化成功 {std.std_ok.sum()}")
    print(f"      去盐 {int(std.was_desalted.sum())} 个；"
          f"中性化 {int(std.was_neutralized.sum())} 个；"
          f"结构有改动 {int(std.std_changed_structure.sum())} 个")

    print("[3/5] 计算描述符 …")
    patterns = {k: Chem.MolFromSmarts(v) for k, v in SMARTS_FLAGS.items()}
    desc = pd.DataFrame([descriptors(s, patterns) for s in std.std_smiles])

    out = pd.concat([df.reset_index(drop=True), std, desc], axis=1)

    print("[4/5] 标注重复结构与骨架簇（只打标，不删行）…")
    ok = out.inchikey.notna()
    out["dup_group"] = out.inchikey.where(ok)
    out["dup_group_size"] = out.groupby("dup_group").inchikey.transform("size")
    out["is_dup_representative"] = ok & ~out.duplicated("dup_group", keep="first")
    out["scaffold_cluster_size"] = (
        out.groupby(out.murcko_scaffold.where(ok))
        .murcko_scaffold.transform("size")
    )
    n_dup = int((out.dup_group_size.fillna(1) > 1).sum())
    print(f"      去盐后结构相同的行 {n_dup} 个"
          f"（分属 {out.loc[out.dup_group_size.fillna(1) > 1, 'dup_group'].nunique()} 组），"
          f"全部保留")

    print("[5/5] 自检并写出 …")
    problems = selfcheck(out)
    if problems:
        for x in problems:
            print(f"  ⚠ {x}")
    else:
        print("  ✓ 全部通过")

    p_out = args.outdir / "Step3_01_RDKit_Processed.csv"
    out.to_csv(p_out, index=False)
    print(f"      {p_out.name}  ({out.shape[0]} × {out.shape[1]})")

    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "input": str(args.input),
        "n_in": int(len(df)),
        "n_out": int(len(out)),
        "parse_ok": int(std.parse_ok.sum()),
        "std_ok": int(std.std_ok.sum()),
        "desalted": int(std.was_desalted.sum()),
        "neutralized": int(std.was_neutralized.sum()),
        "structure_changed": int(std.std_changed_structure.sum()),
        "duplicate_rows_after_std": n_dup,
        "unique_scaffolds": int(out.murcko_scaffold.nunique()),
        "cns_mpo_note": "RDKit 只能给 MW/cLogP/TPSA/HBD 四项；logD7.4 与 pKa 需 ADMETlab",
        "selfcheck_problems": problems,
    }
    (args.outdir / "Step3_01b_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n{len(out)} 行进、{len(out)} 行出（不删行）；"
          f"唯一骨架 {out.murcko_scaffold.nunique()} 种")


if __name__ == "__main__":
    main()
