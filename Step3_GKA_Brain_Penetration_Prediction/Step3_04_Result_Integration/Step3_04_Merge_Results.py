#!/usr/bin/env python
"""
Step3_04：把两个网页工具的结果并进一张表

输入三份：
  Step3_01_RDKit_Processed.csv          1,274 行，本地标准化结构 + RDKit 描述符（骨架表）
  Step3_02_SwissADME_Input/             提交清单 + 回填清单 + 7 批原始结果 + 人工排除记录
  Step3_03_ADMETlab_Input/              提交清单 + 回填清单 + 16 批原始结果

输出一张 1,274 行的整合表，外加三张追溯用的小表。

三条硬规则（都来自 CLAUDE.md 已踩过的坑）：

1. **回填前必须逐行比对「提交的结构」与「返回的结构」的 InChIKey，不符的整行剔除。**
   SwissADME 出过一次返回樟脑的事故（batch6 的 B3D_0441），名称是对的、只有结构错了，
   光核对名称发现不了。ADMETlab 无 ID 字段、只能按行序回填，更需要结构校验兜底。
   本脚本对两个工具都做这件事，剔除的行写进 Step3_04_Verification_Failures.csv。

2. **加列不删行。** 1,274 行进、1,274 行出。没拿到结果的分子留在表里，
   `*_ok = False` 并写明原因，不填任何默认值。

3. **不做任何判定。** 阈值、排序、流程验收都属 Step3_05。
   本步只把数值取全、取对、摆到同一张表上。
   CNS MPO 算**分**（拐点已从原文 PDF 取得，见 Step3_04_CNS_MPO.py，
   六条曲线用原文算例自检），但 **不套 `≥4` 之类的阈值**——那是 Step3_05 的事。

锚点（每批重复提交的同一组 20 个对照分子）在两个工具上都有多次结果，
按「数值取中位数、文本取众数」合并，并把重复间的分歧量写进 Step3_04_Anchor_Drift.csv。

用法：
  micromamba run -n GKA_in_Brain python Step3_04_Merge_Results.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Step3_04_CNS_MPO as cns_mpo          # noqa: E402  拐点取自 cn100008c.pdf

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
STEP3 = HERE.parent
DEFAULT_RDKIT = (STEP3 / "Step3_01_Structure_Standardization_and_RDKit_Properties"
                 / "Step3_01_RDKit_Processed.csv")
DEFAULT_SA_DIR = STEP3 / "Step3_02_SwissADME_Input"
DEFAULT_AD_DIR = STEP3 / "Step3_03_ADMETlab_Input"

N_SA_BATCHES = 7
N_AD_BATCHES = 16

# ADMETlab 的 molstr 是一整段 SVG 结构图（7.6k–36k 字符/分子），
# 不是预测值；结构已由 std_smiles / inchi 完整表达，不进整合表。
ADMET_DROP_COLS = {"molstr"}

# CNS MPO 的六项输入分别取哪一列。
# 主口径的取法与理由见 Step3_04_CNS_MPO.py 的模块文档；三处需要显式选择：
#   · ClogP  —— 原文用 BioByte，本项目没有；取 ADMETlab `logP`，与同源的 logD / pKa 一致
#   · TPSA   —— 原文 ref 9 是 Ertl 2000，通行实现只计 N/O，故取 `tpsa` 而非 `tpsa_sandp`
#               （⚠ 与 BOILED-Egg 相反，那个用含 S/P 的）
#   · MW     —— 取平均分子量 `mw`，**不是** ADMETlab 那个单同位素质量
CNSMPO_INPUTS = {
    "MW": "mw",
    "cLogP": "admetlab_logp",
    "cLogD7.4": "admetlab_logd",
    "TPSA": "tpsa",
    "HBD": "hbd",
    "pKa(most basic)": "admetlab_pka_basic",
}

# 换一个口径就重算一版，让 Step3_05 看得见分数对口径有多敏感。
# key = 列名后缀，value = 覆盖掉主口径的哪一项
CNSMPO_VARIANTS = {
    "tpsa_sandp": {"tpsa": "tpsa_sandp"},                     # TPSA 含 S/P
    "logp_rdkit": {"clogp": "clogp"},                          # RDKit Crippen
    "logp_swiss": {"clogp": "swissadme_consensus_log_p"},      # SwissADME 四法均值
}

# 结构自检：787 个 GKA 候选里应当一个不少的 12 个身份对照（不用于入脑阈值）
N_GKA_POSITIVE_CONTROLS = 12


# --------------------------------------------------------------------------- #
# 工具函数
# --------------------------------------------------------------------------- #
def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def slug(name: str) -> str:
    """列名 → snake_case。'#Heavy atoms' → n_heavy_atoms，'t0.5' → t0_5。"""
    s = name.replace("#", "n_")
    s = re.sub(r"[^0-9a-zA-Z]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_").lower()
    return s


def inchikey(smiles) -> str | None:
    if not isinstance(smiles, str) or not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m is not None else None


def to_numeric_where_possible(df: pd.DataFrame, skip: set[str]) -> pd.DataFrame:
    """能整列转成数值的就转，转不动的保持原样（分类标签、原子序号列表等）。"""
    out = df.copy()
    for c in out.columns:
        if c in skip or out[c].dtype.kind in "fiu":
            continue
        conv = pd.to_numeric(out[c], errors="coerce")
        nonnull = out[c].notna()
        if nonnull.any() and conv[nonnull].notna().all():
            out[c] = conv
    return out


def consolidate(g: pd.DataFrame, value_cols: list[str]) -> tuple[dict, dict]:
    """把同一结构的多次结果（锚点）合成一行。

    数值取中位数、文本取众数；同时返回逐列的分歧量，供漂移表使用。
    """
    row, disagree = {}, {}
    for c in value_cols:
        s = g[c]
        vals = s.dropna()
        if vals.empty:
            row[c] = np.nan
            continue
        if s.dtype.kind in "fiu":
            row[c] = float(np.median(vals))
            spread = float(vals.max() - vals.min())
            if spread > 0:
                disagree[c] = {"n": int(len(vals)), "spread": spread,
                               "min": float(vals.min()), "max": float(vals.max())}
        else:
            vc = vals.value_counts()
            row[c] = vc.index[0]
            if len(vc) > 1:
                disagree[c] = {"n": int(len(vals)), "values": vc.to_dict()}
    return row, disagree


# --------------------------------------------------------------------------- #
# 解析 + 校验：SwissADME
# --------------------------------------------------------------------------- #
def parse_swissadme(d: Path, n_batches: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """读 7 批原始结果，按 (batch, 名称) 对回提交清单，逐行核对结构。

    返回 (通过校验的长表, 校验失败记录, 每批出处)。
    SwissADME 的名称字段就是 mol_id（Step3_02 特意这么提交的），
    但**名称对不代表结构对**——必须比 InChIKey。
    """
    man = pd.read_csv(d / "Step3_02_Submission_Manifest.csv")
    sub_map = man.set_index(["batch", "submitted_name"]).std_smiles.to_dict()

    frames, prov, fails = [], [], []
    for b in range(1, n_batches + 1):
        p = d / f"swissadme_batch{b}.csv"
        if not p.exists():
            sys.exit(f"[FATAL] 缺少 SwissADME 第 {b} 批结果：{p}")
        r = pd.read_csv(p)
        r.insert(0, "batch", b)
        r.insert(1, "mol_id", r["Molecule"])
        frames.append(r)
        prov.append({
            "tool": "SwissADME", "batch": b, "result_file": p.name,
            "n_submitted": int((man.batch == b).sum()), "n_returned": int(len(r)),
            "file_mtime": datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "sha256": sha256(p),
        })
    res = pd.concat(frames, ignore_index=True)

    # "n/d" 是工具自己给的「没算出来」，转成空值；哪些列整列 n/d 记进 summary
    res = res.replace("n/d", np.nan)

    res["submitted_smiles"] = [sub_map.get((b, n)) for b, n in zip(res.batch, res.mol_id)]
    res["returned_ik"] = res["Canonical SMILES"].map(inchikey)
    res["submitted_ik"] = res.submitted_smiles.map(inchikey)

    ok = res.submitted_ik.notna() & (res.returned_ik == res.submitted_ik)
    dup = res.duplicated(["batch", "mol_id"], keep=False)
    for _, r in res[~ok | dup].iterrows():
        fails.append({
            "tool": "SwissADME", "batch": int(r.batch), "line_no": np.nan,
            "mol_id": r.mol_id,
            "reason": ("名称在该批提交清单中不存在" if pd.isna(r.submitted_smiles)
                       else "同一批内名称重复" if dup.loc[r.name]
                       else "返回结构与提交结构不符"),
            "submitted_smiles": r.submitted_smiles,
            "returned_smiles": r["Canonical SMILES"],
            "submitted_inchikey": r.submitted_ik, "returned_inchikey": r.returned_ik,
        })
    res = res[ok & ~dup].copy()

    # 每批返回的行数与提交数对不上，要留痕（batch3 少 1 = 人工排除的 B3D_0012）
    for row in prov:
        got = int((res.batch == row["batch"]).sum())
        row["n_accepted"] = got
        row["n_rejected"] = row["n_returned"] - got
    return res, pd.DataFrame(fails), pd.DataFrame(prov)


# --------------------------------------------------------------------------- #
# 解析 + 校验：ADMETlab
# --------------------------------------------------------------------------- #
def parse_admetlab(d: Path, n_batches: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """读 16 批原始结果，按 (batch, 行序) 对回提交清单，再用 raw_smiles 独立校验。

    ADMETlab 不接受名称字段，回填只能靠行序——行序一旦错位后面全错且不报错。
    好在它把输入原样放在 raw_smiles 列返回，于是有了独立于行序的第二条校验路径。
    """
    man = pd.read_csv(d / "Step3_03_Submission_Manifest.csv")

    frames, prov, fails = [], [], []
    for b in range(1, n_batches + 1):
        p = d / f"ADMetlab_batch{b}.csv"
        if not p.exists():
            sys.exit(f"[FATAL] 缺少 ADMETlab 第 {b} 批结果：{p}")
        r = pd.read_csv(p)
        r = pd.concat([pd.DataFrame({"batch": b, "line_no": range(1, len(r) + 1)}), r], axis=1)
        n_sub = int((man.batch == b).sum())
        if len(r) != n_sub:
            fails.append({
                "tool": "ADMETlab", "batch": b, "line_no": np.nan, "mol_id": "",
                "reason": f"该批返回 {len(r)} 行，提交 {n_sub} 行——行数不符，行序回填不可信",
                "submitted_smiles": "", "returned_smiles": "",
                "submitted_inchikey": "", "returned_inchikey": "",
            })
        frames.append(r)
        prov.append({
            "tool": "ADMETlab", "batch": b, "result_file": p.name,
            "n_submitted": n_sub, "n_returned": int(len(r)),
            "file_mtime": datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "sha256": sha256(p),
        })
    res = pd.concat(frames, ignore_index=True)

    res = res.merge(man[["batch", "line_no", "mol_id", "std_smiles"]],
                    on=["batch", "line_no"], how="left")
    res = res.rename(columns={"std_smiles": "submitted_smiles"})
    res["returned_ik"] = res.raw_smiles.map(inchikey)
    res["submitted_ik"] = res.submitted_smiles.map(inchikey)

    ok = res.submitted_ik.notna() & (res.returned_ik == res.submitted_ik)
    for _, r in res[~ok].iterrows():
        fails.append({
            "tool": "ADMETlab", "batch": int(r.batch), "line_no": int(r.line_no),
            "mol_id": r.mol_id if isinstance(r.mol_id, str) else "",
            "reason": ("该行序在提交清单中不存在" if pd.isna(r.submitted_smiles)
                       else "返回结构（raw_smiles）与提交结构不符"),
            "submitted_smiles": r.submitted_smiles, "returned_smiles": r.raw_smiles,
            "submitted_inchikey": r.submitted_ik, "returned_inchikey": r.returned_ik,
        })
    res = res[ok].copy()

    for row in prov:
        got = int((res.batch == row["batch"]).sum())
        row["n_accepted"] = got
        row["n_rejected"] = row["n_returned"] - got
    return res, pd.DataFrame(fails), pd.DataFrame(prov)


# --------------------------------------------------------------------------- #
# 合并一个工具的结果
# --------------------------------------------------------------------------- #
def collapse(res: pd.DataFrame, value_cols: list[str], prefix: str,
             tool: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """把长表按 mol_id 收成一行一分子，并产出锚点漂移表。"""
    rows, drift = [], []
    for mid, g in res.groupby("mol_id", sort=True):
        vals, disagree = consolidate(g, value_cols)
        rec = {"mol_id": mid, "n_replicates": int(len(g)),
               "batches": json.dumps(sorted(int(b) for b in g.batch))}
        rec.update({f"{prefix}{slug(c)}": v for c, v in vals.items()})
        rows.append(rec)
        if len(g) > 1:
            worst = max((d.get("spread", np.inf) for d in disagree.values()), default=0.0)
            drift.append({
                "tool": tool, "mol_id": mid, "n_replicates": int(len(g)),
                "batches": json.dumps(sorted(int(b) for b in g.batch)),
                "n_cols_compared": len(value_cols),
                "n_cols_disagree": len(disagree),
                "max_abs_spread": worst if disagree else 0.0,
                "detail": json.dumps(disagree, ensure_ascii=False, default=float),
            })
    return pd.DataFrame(rows), pd.DataFrame(drift)


def add_cns_mpo(out: pd.DataFrame) -> pd.DataFrame:
    """算 CNS MPO：主口径 + 三个换口径的敏感性版本。

    任一项缺失总分即为 NaN，**不补零**——补零等于把「没测到」说成「最差」。
    """
    def cols_for(override: dict[str, str] | None = None) -> dict[str, str]:
        m = {"mw": CNSMPO_INPUTS["MW"], "clogp": CNSMPO_INPUTS["cLogP"],
             "clogd": CNSMPO_INPUTS["cLogD7.4"], "tpsa": CNSMPO_INPUTS["TPSA"],
             "hbd": CNSMPO_INPUTS["HBD"], "pka": CNSMPO_INPUTS["pKa(most basic)"]}
        if override:
            m.update(override)
        return m

    m = cols_for()
    out["cnsmpo_inputs_complete"] = out[list(m.values())].notna().all(axis=1)
    t0 = cns_mpo.score_frame(**{k: out[v] for k, v in m.items()})
    for p in cns_mpo.PROPERTIES:
        out[f"cnsmpo_t0_{p}"] = t0[f"t0_{p}"].values
    out["cnsmpo_score"] = t0["score"].values

    variant_cols = ["cnsmpo_score"]
    for name, override in CNSMPO_VARIANTS.items():
        mv = cols_for(override)
        if any(c not in out.columns for c in mv.values()):
            continue
        c = f"cnsmpo_score_{name}"
        out[c] = cns_mpo.score_frame(**{k: out[v] for k, v in mv.items()})["score"].values
        variant_cols.append(c)
    # 分数对口径有多敏感，逐分子给出，供 Step3_05 判断哪些结论经得起换口径
    out["cnsmpo_score_variant_spread"] = (out[variant_cols].max(axis=1)
                                          - out[variant_cols].min(axis=1))
    return out


def attach(base: pd.DataFrame, join_map: pd.DataFrame, collapsed: pd.DataFrame,
           prefix: str, excluded: dict[str, str]) -> pd.DataFrame:
    """按回填清单把结果贴到 1,274 行上。

    join_map 的 submitted_as 指出「本行的结构由哪个 mol_id 的那次提交覆盖」——
    去盐后结构相同的分子（MK-0941 / Globalagliatin 的盐与游离碱）共用一次提交。
    """
    jm = join_map.set_index("mol_id")
    src = base.mol_id.map(jm.submitted_as)
    out = base[["mol_id"]].copy()
    out[f"{prefix}result_from"] = src
    merged = out.merge(collapsed, left_on=f"{prefix}result_from", right_on="mol_id",
                       how="left", suffixes=("", "_r")).drop(columns=["mol_id_r"])
    merged = merged.rename(columns={"n_replicates": f"{prefix}n_replicates",
                                    "batches": f"{prefix}batches"})
    got = merged[f"{prefix}n_replicates"].notna()
    merged[f"{prefix}ok"] = got
    merged[f"{prefix}n_replicates"] = merged[f"{prefix}n_replicates"].astype("Int64")

    reason = pd.Series("", index=merged.index, dtype=object)
    for i, mid in enumerate(merged.mol_id):
        if got.iloc[i]:
            continue
        if mid in excluded:
            reason.iloc[i] = excluded[mid]
        elif pd.isna(src.iloc[i]):
            reason.iloc[i] = "未提交（回填清单中无 submitted_as）"
        else:
            reason.iloc[i] = f"提交了但结果未通过结构校验或未返回（submitted_as={src.iloc[i]}）"
    merged[f"{prefix}missing_reason"] = reason
    merged[f"{prefix}shared_result_with"] = np.where(
        got & (merged[f"{prefix}result_from"] != merged.mol_id),
        merged[f"{prefix}result_from"], "")
    return merged.drop(columns=["mol_id"])


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rdkit", type=Path, default=DEFAULT_RDKIT)
    ap.add_argument("--swissadme-dir", type=Path, default=DEFAULT_SA_DIR)
    ap.add_argument("--admetlab-dir", type=Path, default=DEFAULT_AD_DIR)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    started = datetime.now()
    if not args.rdkit.exists():
        sys.exit(f"[FATAL] 找不到骨架表：{args.rdkit}（先跑 Step3_01b）")

    print("[1/7] 读骨架表 …")
    base = pd.read_csv(args.rdkit)
    n_in = len(base)
    print(f"      {n_in} 行 × {base.shape[1]} 列"
          f"（{' / '.join(f'{k} {v}' for k, v in base.set.value_counts().items())}）")

    # 人工排除留痕：与「工具没返回」要分得开，否则对账会误报
    excl_path = args.swissadme_dir / "Step3_02_Manual_Exclusions.csv"
    excluded: dict[str, str] = {}
    if excl_path.exists():
        ex = pd.read_csv(excl_path)
        excluded = {r.mol_id: f"人工排除：{r.reason}（{r.excluded_at} / {r.excluded_by}）"
                    for r in ex.itertuples()}
        print(f"      人工排除记录 {len(excluded)} 条：{', '.join(excluded)}")

    print("[2/7] 解析 SwissADME 结果并逐行核对结构 …")
    sa, sa_fail, sa_prov = parse_swissadme(args.swissadme_dir, N_SA_BATCHES)
    print(f"      返回 {int(sa_prov.n_returned.sum())} 行，"
          f"通过校验 {len(sa)}，剔除 {len(sa_fail)}")
    for _, r in sa_fail.iterrows():
        print(f"      ⚠ 剔除 batch{r.batch} {r.mol_id}：{r.reason}")

    print("[3/7] 解析 ADMETlab 结果并逐行核对结构 …")
    ad, ad_fail, ad_prov = parse_admetlab(args.admetlab_dir, N_AD_BATCHES)
    print(f"      返回 {int(ad_prov.n_returned.sum())} 行，"
          f"通过校验 {len(ad)}，剔除 {len(ad_fail)}")
    for _, r in ad_fail.iterrows():
        print(f"      ⚠ 剔除 batch{r.batch} 第 {r.line_no} 行：{r.reason}")

    print("[4/7] 合并锚点重复 …")
    sa_meta = {"batch", "mol_id", "Molecule", "submitted_smiles",
               "returned_ik", "submitted_ik"}
    sa_value_cols = [c for c in sa.columns if c not in sa_meta]
    sa = to_numeric_where_possible(sa, skip=sa_meta | {"Canonical SMILES", "Formula"})
    sa_flat, sa_drift = collapse(sa, sa_value_cols, "swissadme_", "SwissADME")

    ad_meta = {"batch", "line_no", "mol_id", "submitted_smiles",
               "returned_ik", "submitted_ik"}
    ad_value_cols = [c for c in ad.columns if c not in ad_meta | ADMET_DROP_COLS]
    ad = to_numeric_where_possible(ad, skip=ad_meta | {"raw_smiles", "smiles"})
    ad_flat, ad_drift = collapse(ad, ad_value_cols, "admetlab_", "ADMETlab")
    print(f"      SwissADME {len(sa_flat)} 个唯一结构 / ADMETlab {len(ad_flat)} 个")
    for tool, dr in (("SwissADME", sa_drift), ("ADMETlab", ad_drift)):
        if len(dr):
            n_bad = int((dr.n_cols_disagree > 0).sum())
            print(f"      {tool} 锚点 {len(dr)} 个，重复间有分歧的 {n_bad} 个"
                  f"（最大绝对差 {dr.max_abs_spread.max():.3g}）")

    print("[5/7] 贴回 1,274 行 …")
    sa_join = pd.read_csv(args.swissadme_dir / "Step3_02_Result_Join_Map.csv")
    ad_join = pd.read_csv(args.admetlab_dir / "Step3_03_Result_Join_Map.csv")
    sa_part = attach(base, sa_join, sa_flat, "swissadme_", excluded)
    ad_part = attach(base, ad_join, ad_flat, "admetlab_", excluded)
    out = pd.concat([base, sa_part, ad_part], axis=1)

    print("      CNS MPO（拐点取自 cn100008c.pdf，先自检六条曲线）…")
    problems_mpo = cns_mpo.selfcheck()
    if problems_mpo:
        for x in problems_mpo:
            print(f"      ⚠ {x}")
        sys.exit("[FATAL] CNS MPO 曲线自检不通过，拒绝出分")
    out = add_cns_mpo(out)

    print("[6/7] 自检 …")
    problems = []
    if len(out) != n_in:
        problems.append(f"行数变了：进 {n_in} 出 {len(out)}")
    if out.mol_id.duplicated().any():
        problems.append("mol_id 出现重复")
    n_ctrl = int(out.get("gka_is_positive_control", pd.Series(dtype=bool)).fillna(False).sum())
    if n_ctrl != N_GKA_POSITIVE_CONTROLS:
        problems.append(f"GKA 身份对照 {n_ctrl} 个，应为 {N_GKA_POSITIVE_CONTROLS} 个")
    for pre, name in (("swissadme_", "SwissADME"), ("admetlab_", "ADMETlab")):
        miss = out[~out[f"{pre}ok"]]
        unexplained = miss[miss[f"{pre}missing_reason"] == ""]
        if len(unexplained):
            problems.append(f"{name} 有 {len(unexplained)} 行缺结果且没写原因")
        not_manual = [m for m in miss.mol_id if m not in excluded]
        if not_manual:
            problems.append(f"{name} 有非人工排除的缺口：{not_manual[:10]}")
    # 用工具自己算的 MW 交叉校验本地 RDKit 的 MW，差得多说明贴错了行。
    # ⚠ 两个工具的 MW 口径不同：SwissADME 给平均分子量，ADMETlab 给单同位素质量
    #   （实测 admetlab_mw 与 mw_exact 全表一致；含 Cl/Br 的分子与平均分子量差可达 3 Da）。
    #   比错基准会得出「73 行回填错位」的假警报。
    for pre, ref in (("swissadme_", "mw"), ("admetlab_", "mw_exact")):
        c = f"{pre}mw"
        if c in out.columns:
            d = (out[c] - out[ref]).abs()
            bad = int((d > 1.0).sum())
            if bad:
                problems.append(f"{c} 与本地 {ref} 相差 >1 Da 的有 {bad} 行（回填可能错位）")
            else:
                print(f"      ✓ {c} 与本地 {ref} 全部一致（最大差 {d.max():.3g} Da）")
    if problems:
        for x in problems:
            print(f"  ⚠ {x}")
    else:
        print("  ✓ 全部通过")

    print("[7/7] 写出 …")
    args.outdir.mkdir(parents=True, exist_ok=True)
    p_main = args.outdir / "Step3_04_Integrated_Brain_Penetration_Results.csv"
    out.to_csv(p_main, index=False)
    fails = pd.concat([sa_fail, ad_fail], ignore_index=True) if len(sa_fail) or len(ad_fail) \
        else pd.DataFrame(columns=["tool", "batch", "line_no", "mol_id", "reason"])
    drift = pd.concat([sa_drift, ad_drift], ignore_index=True)
    prov = pd.concat([sa_prov, ad_prov], ignore_index=True)
    fails.to_csv(args.outdir / "Step3_04_Verification_Failures.csv", index=False)
    drift.to_csv(args.outdir / "Step3_04_Anchor_Drift.csv", index=False)
    prov.to_csv(args.outdir / "Step3_04_Batch_Provenance.csv", index=False)
    print(f"      {p_main.name}  {out.shape[0]} 行 × {out.shape[1]} 列")

    sa_allnd = [f"swissadme_{slug(c)}" for c in sa_value_cols
                if f"swissadme_{slug(c)}" in out.columns
                and out.loc[out.swissadme_ok, f"swissadme_{slug(c)}"].isna().all()]
    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "rdkit_table": str(args.rdkit),
            "swissadme_dir": str(args.swissadme_dir),
            "admetlab_dir": str(args.admetlab_dir),
        },
        "tools": {
            "SwissADME": {"url": "https://www.swissadme.ch/", "version": "网页版未标版本号",
                          "batches": N_SA_BATCHES,
                          "note": "iLOGP 与 5 个 CYP 抑制预测整列 n/d（ChemAxon 支持终止的后果）"},
            "ADMETlab": {"url": "https://admetlab3.scbdd.com/", "version": "3.0",
                         "batches": N_AD_BATCHES,
                         "note": "BBB 列是 0–1 概率，不是 Yes/No 标签"},
        },
        "n_rows_in": int(n_in), "n_rows_out": int(len(out)),
        "n_cols_out": int(out.shape[1]),
        "swissadme": {
            "n_returned": int(sa_prov.n_returned.sum()),
            "n_accepted": int(len(sa)), "n_rejected": int(len(sa_fail)),
            "n_unique_structures": int(len(sa_flat)),
            "n_rows_covered": int(out.swissadme_ok.sum()),
            "all_nd_columns": sa_allnd,
        },
        "admetlab": {
            "n_returned": int(ad_prov.n_returned.sum()),
            "n_accepted": int(len(ad)), "n_rejected": int(len(ad_fail)),
            "n_unique_structures": int(len(ad_flat)),
            "n_rows_covered": int(out.admetlab_ok.sum()),
        },
        "manual_exclusions": excluded,
        "cnsmpo": {
            "computed": True,
            "source": "Wager TT et al. ACS Chem Neurosci 2010;1(6):435-449, "
                      "doi:10.1021/cn100008c —— Table 1 + Figure 4（拐点）、"
                      "Methods eq 1/2（分段线性 + 等权求和）；本地 PDF cn100008c.pdf",
            "curve_selfcheck": "对照原文 Table 4 算例与 Table 3 的三个候选，"
                               "六项 T0 与总分全部复现（不通过则脚本直接退出）",
            "inflection_points": {k: v for k, v in cns_mpo.INFLECTIONS.items()},
            "input_columns": CNSMPO_INPUTS,
            "variants": {k: v for k, v in CNSMPO_VARIANTS.items()},
            "n_rows_inputs_complete": int(out.cnsmpo_inputs_complete.sum()),
            "n_rows_scored": int(out.cnsmpo_score.notna().sum()),
            "note": "原文 p.446 明写 the algorithm is not intended to be used purely as "
                    "a predictor of CNS penetration —— 它是成药性对齐工具，不是入脑预测器。"
                    "阈值（如常被引用的 ≥4）属 Step3_05，本步不套。",
        },
        "selfcheck_problems": problems,
    }
    (args.outdir / "Step3_04_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n完成：{len(out)} 行 × {out.shape[1]} 列，"
          f"SwissADME 覆盖 {int(out.swissadme_ok.sum())} 行 / "
          f"ADMETlab 覆盖 {int(out.admetlab_ok.sum())} 行，"
          f"用时 {(datetime.now() - started).total_seconds():.1f}s")
    print("⚠ 本步产物是数据，不是结论：阈值、排序、流程验收都属 Step3_05。")


if __name__ == "__main__":
    main()
