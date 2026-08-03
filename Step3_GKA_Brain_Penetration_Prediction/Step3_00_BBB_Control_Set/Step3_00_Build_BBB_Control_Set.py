#!/usr/bin/env python
"""
Step3_00：构建入脑（BBB）对照集

两个来源，角色不同：

  1. B3DB regression  → 正负分类对照（100 BBB+ / 56 BBB-）
     判据 logBB，用于「分得开 / 分不开」的二分类检验。
  2. Fridén 2009      → 全部定量脑暴露参考（K(p,uu,brain)）
     判据 Kp,uu，数值连续，用于看预测值能否重现真实的暴露梯度。

选取规则（来自设计约定）：
  BBB+ : regression group B，logBB >= -0.5，每骨架 1 个，取 100
  BBB- : regression group B，logBB <= -1.1，每骨架 1 个，全取（56 个）
  两侧排序均优先靠近 GKA 候选的理化空间

⚠ 三点已量出来的约束，改参数前先看 README：
  - BBB- 侧 group B 池只有 76 个分子 / 56 个骨架，两侧数量不对称是有意的：
    宁可少，也不为凑 100 而放弃骨架去冗余或掺入单一来源的 group A
  - BBB+ 池中位 MW 292 / TPSA 42，GKA 候选是 463 / 105，两者基本不重叠，
    「靠近 GKA 空间」只能作次级排序，不能作硬筛
  - group 的 A/B/C/D 含义 B3DB 官方文档未逐个定义，此处按 reference 数量的实测特征理解

用法：
  micromamba run -n GKA_in_Brain python Step3_00_Build_BBB_Control_Set.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------- 配置

B3DB_REGRESSION = HERE / "B3DB_regression.tsv"
FRIDEN_RAW = HERE / "Friden2009_CHEMBL1798466_raw.tsv"
GKA_CANDIDATES = (
    HERE.parents[1] / "Step1_Find_GKA_from_ChEMBL"
    / "Step1_GKA_Candidates_with_Properties.csv"
)

SHA256 = {
    "B3DB_regression.tsv":
        "1be4e33ab1fa1d99897541b6e0a9a00cd6061cf970a36dc838e640407735eba1",
    "Friden2009_CHEMBL1798466_raw.tsv":
        "8926fad7dceb502e1c2d63b7d86a078c5887de06c4cbe1e78ec1cd078e6fc464",
}

# 选取规则
POS_GROUP = "B"
POS_LOGBB_MIN = -0.5
NEG_GROUPS = ["B", "A"]          # 顺序即优先级：先取尽 group B，不足再用 group A 补
NEG_LOGBB_MAX = -1.1

# ⚠ 不设固定目标数：合格的都要，宁可多留可靠分子。
#   固定 N 会逼着在「凑数」与「保多样性」之间二选一，两次都得不偿失。
TARGET_POS = None
TARGET_NEG = None

# 骨架去冗余是**软性**的：同一 Murcko 骨架最多保留 SCAFFOLD_CAP 个代表。
# 不做「一骨架一分子」的硬去重——那会把证据充分的同系物一并砍掉。
SCAFFOLD_CAP = 3

# 簇内排序的权重。前两项决定「留谁」，GKA 邻近度只作小幅倾斜：
#   logbb_margin  离 logBB = -1 分类边界越远，类别归属越确定
#   n_reference   独立文献来源越多，证据越可靠（group B 多来源，group A 恒为 1）
#   gka_distance  离 GKA 候选理化空间越近越好（负权重）
W_MARGIN, W_NREF, W_GKA = 1.0, 1.0, -0.5

# 无环分子的 Murcko 骨架是空串。True = 全部当同一簇；False = 各自独立成簇。
# 阴性侧只差 1 个，但阳性池有 70 个无环分子，True 时它们最多只能入选 cap 个。
ACYCLIC_AS_ONE_CLUSTER = True

# 用于计算「与 GKA 候选的理化距离」的描述符
PROP_COLS = ["mw", "tpsa", "clogp", "hbd", "hba", "rtb"]
GKA_PROP_MAP = {                 # GKA 主表里的对应列
    "mw": "full_mwt", "tpsa": "psa", "clogp": "alogp",
    "hbd": "hbd", "hba": "hba", "rtb": "rtb",
}


# ---------------------------------------------------------------- 工具


def check_sha(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = SHA256.get(path.name)
    if expected and digest != expected:
        print(f"  ⚠ {path.name} SHA256 与记录不符，版本可能已变更")
        print(f"     记录：{expected}")
        print(f"     实际：{digest}")
    return digest


def add_structure_columns(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """解析结构、算 Murcko 骨架与 RDKit 描述符。解析失败的保留并标注。"""
    df = df.copy()
    mols = [Chem.MolFromSmiles(s) if isinstance(s, str) else None
            for s in df[smiles_col]]
    df["parse_ok"] = [m is not None for m in mols]

    def _s(m):
        if m is None:
            return None
        try:
            return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        except Exception:
            return None

    df["murcko_scaffold"] = [_s(m) for m in mols]
    df["inchikey"] = [Chem.MolToInchiKey(m) if m is not None else None for m in mols]
    for name, fn in [("mw", Descriptors.MolWt), ("tpsa", Descriptors.TPSA),
                     ("clogp", Crippen.MolLogP), ("hbd", Descriptors.NumHDonors),
                     ("hba", Descriptors.NumHAcceptors),
                     ("rtb", Descriptors.NumRotatableBonds)]:
        df[name] = [fn(m) if m is not None else np.nan for m in mols]
    # 无环分子的 Murcko 骨架是空串，不能和其它分子归成一簇
    df["is_acyclic"] = df.murcko_scaffold == ""
    return df


def gka_proximity(df: pd.DataFrame, gka: pd.DataFrame) -> pd.Series:
    """到 GKA 候选理化中心的标准化欧氏距离；越小越像 GKA 空间。

    ⚠ 只作次级排序用。BBB+ 分子天然小而低极性，与 GKA 空间基本不重叠，
    把它当硬筛会把候选池砍到凑不够数。
    """
    ref = gka[[GKA_PROP_MAP[c] for c in PROP_COLS]].astype(float)
    mu, sd = ref.mean().values, ref.std().values
    sd = np.where(sd == 0, 1.0, sd)
    x = df[PROP_COLS].astype(float).values
    return pd.Series(np.sqrt((((x - mu) / sd) ** 2).mean(axis=1)), index=df.index)


def cluster_key(pool: pd.DataFrame) -> pd.Series:
    """骨架簇键。无环分子的 Murcko 骨架是空串，是否合并成一簇由开关决定。"""
    if ACYCLIC_AS_ONE_CLUSTER:
        return pool.murcko_scaffold
    return pd.Series(
        np.where(pool.is_acyclic, "ACYCLIC::" + pool.index.astype(str),
                 pool.murcko_scaffold),
        index=pool.index,
    )


def quality_score(pool: pd.DataFrame) -> pd.Series:
    """簇内排序用的综合分，越大越优先保留。

    三项各自在**本池内**标准化后加权求和，避免量纲不同的一项独霸排序。
    """
    def z(v: pd.Series) -> pd.Series:
        sd = v.std()
        return (v - v.mean()) / (sd if sd and sd > 0 else 1.0)

    return (
        W_MARGIN * z(pool.logbb_margin)
        + W_NREF * z(np.log1p(pool.n_reference))
        + W_GKA * z(pool.gka_distance)
    )


def pick_with_scaffold_cap(pool: pd.DataFrame, n: int | None) -> pd.DataFrame:
    """同一骨架最多留 SCAFFOLD_CAP 个，簇内按综合分取优。

    n=None 表示不设目标数——合格的都要。给定 n 时只截取综合分最高的 n 个。
    """
    pool = pool.copy()
    pool["_cluster"] = cluster_key(pool)
    pool["quality_score"] = quality_score(pool)
    pool = pool.sort_values("quality_score", ascending=False)
    picked = pool.groupby("_cluster", group_keys=False).head(SCAFFOLD_CAP)
    if n is not None:
        picked = picked.head(n)
    out = picked.copy()
    out["scaffold_cap_used"] = SCAFFOLD_CAP
    return out.drop(columns="_cluster")


# ---------------------------------------------------------------- 选取


def select_b3db(reg: pd.DataFrame, gka: pd.DataFrame,
                exclude_inchikeys: set[str] | None = None) -> pd.DataFrame:
    """从 B3DB regression 选 BBB+ / BBB-，全表返回并打 selected 标记（不删行）。"""
    reg = reg.copy()
    reg["gka_distance"] = gka_proximity(reg, gka)
    # 独立文献来源数：reference 列形如 "R18|R26|R27|"，末尾有空段要去掉
    reg["n_reference"] = (
        reg.reference.fillna("").str.strip("|").str.split("|")
        .map(lambda xs: len([x for x in xs if x]))
    )
    # 离 logBB = -1 分类边界的余量，越大类别归属越确定
    reg["logbb_margin"] = (reg.logBB - (-1.0)).abs()
    reg["selected"] = False
    reg["control_class"] = ""
    reg["selection_note"] = ""
    reg["exclude_reason"] = ""

    # 与 Fridén 定量参考集去重：同一分子不在两个集合里各出现一次。
    # 从 B3DB 侧剔除而不是 Fridén 侧——Fridén 只有 42 个且带定量 Kp,uu，更金贵；
    # 在**选取之前**剔除，这样同骨架簇会自动补上次优的那个，不会留下空位。
    reg["dropped_dup_with_friden"] = False
    ok = reg.parse_ok
    if exclude_inchikeys:
        dup = reg.inchikey.isin(exclude_inchikeys)
        reg.loc[dup, "dropped_dup_with_friden"] = True
        ok = ok & ~dup

    # ---- BBB+ ----
    pos_pool = reg[ok & (reg.group == POS_GROUP) & (reg.logBB >= POS_LOGBB_MIN)]
    pos = pick_with_scaffold_cap(pos_pool, TARGET_POS)
    reg.loc[pos.index, ["selected", "control_class"]] = True, "control_positive"
    reg.loc[pos.index, "selection_note"] = (
        f"group {POS_GROUP}, logBB>={POS_LOGBB_MIN}, scaffold_cap={SCAFFOLD_CAP}"
    )

    # ---- BBB- ----
    # ⚠ 骨架上限必须在**合并池**上只应用一次。
    #   按 group 分轮各自设限时，同一骨架会在每轮各拿满 cap 个，上限形同虚设。
    #   group B 优先无需单独一轮：quality_score 含 n_reference，
    #   而 group A 恒为单一来源（n_reference = 1），簇内自然排在 group B 之后。
    neg_pool = reg[ok & reg.group.isin(NEG_GROUPS) & (reg.logBB <= NEG_LOGBB_MAX)]
    neg = pick_with_scaffold_cap(neg_pool, TARGET_NEG)
    reg.loc[neg.index, ["selected", "control_class"]] = True, "control_negative"
    reg.loc[neg.index, "selection_note"] = (
        "group " + reg.loc[neg.index, "group"].astype(str)
        + f", logBB<={NEG_LOGBB_MAX}, scaffold_cap={SCAFFOLD_CAP}"
    )

    # ---- 未入选的原因 ----
    unsel = ~reg.selected
    reg.loc[unsel & ~reg.parse_ok, "exclude_reason"] = "SMILES 无法解析"
    reg.loc[unsel & reg.dropped_dup_with_friden, "exclude_reason"] = (
        "与 Fridén 定量参考集重复，已在该侧保留"
    )
    reg.loc[unsel & ok & (reg.exclude_reason == "") & reg.logBB.between(NEG_LOGBB_MAX, POS_LOGBB_MIN,
                                           inclusive="neither"),
            "exclude_reason"] = f"logBB 落在 ({NEG_LOGBB_MAX}, {POS_LOGBB_MIN}) 边界带"
    reg.loc[unsel & ok & (reg.exclude_reason == "") & ~reg.group.isin(NEG_GROUPS),
            "exclude_reason"] = "group 不在选取范围"
    reg.loc[unsel & ok & (reg.exclude_reason == ""),
            "exclude_reason"] = "同骨架已有代表 / 超出目标数"
    return reg


def load_friden(path: Path) -> pd.DataFrame:
    """Fridén 2009 定量参考集，全部保留，不做筛选。"""
    raw = pd.read_csv(path, sep="\t")
    df = raw.rename(columns={
        "standard_value": "kpuu_brain",
        "standard_relation": "kpuu_relation",
        "canonical_smiles": "smiles",
    })[["molecule_chembl_id", "compound_name", "kpuu_brain", "kpuu_relation",
        "smiles", "assay_chembl_id", "doc_chembl_id", "pubmed_id", "doi",
        "journal", "year", "assay_organism"]]
    return df.sort_values(["kpuu_brain", "molecule_chembl_id"],
                          ascending=[False, True],
                          na_position="last").reset_index(drop=True)


# ---------------------------------------------------------------- 自检


def selfcheck(b3db_sel: pd.DataFrame, friden: pd.DataFrame) -> list[str]:
    p = []
    npos = (b3db_sel.control_class == "control_positive").sum()
    nneg = (b3db_sel.control_class == "control_negative").sum()
    if TARGET_POS is not None and npos != TARGET_POS:
        p.append(f"BBB+ 选到 {npos} 个，目标 {TARGET_POS}")
    if TARGET_NEG is not None and nneg != TARGET_NEG:
        p.append(f"BBB- 选到 {nneg} 个，目标 {TARGET_NEG}")
    if npos < 30 or nneg < 30:
        p.append(f"某一类样本过少（BBB+ {npos} / BBB- {nneg}），二分类检验会不稳")
    # 骨架软上限不得被突破
    for cls in ("control_positive", "control_negative"):
        sub = b3db_sel[b3db_sel.control_class == cls]
        over = sub.murcko_scaffold.value_counts()
        bad = over[over > SCAFFOLD_CAP]
        if len(bad):
            p.append(f"{cls} 有 {len(bad)} 个骨架超过上限 {SCAFFOLD_CAP}："
                     f"最大 {bad.iloc[0]} 个")

    pos = b3db_sel[b3db_sel.control_class == "control_positive"]
    neg = b3db_sel[b3db_sel.control_class == "control_negative"]
    if len(pos) and pos.logBB.min() < POS_LOGBB_MIN:
        p.append(f"BBB+ 中出现 logBB < {POS_LOGBB_MIN}：{pos.logBB.min()}")
    if len(neg) and neg.logBB.max() > NEG_LOGBB_MAX:
        p.append(f"BBB- 中出现 logBB > {NEG_LOGBB_MAX}：{neg.logBB.max()}")

    dup = b3db_sel[b3db_sel.inchikey.notna()].inchikey.duplicated().sum()
    if dup:
        p.append(f"B3DB 入选集内有 {dup} 个重复 InChIKey")

    # 去重后两个集合不应再有交集
    ov = b3db_sel[b3db_sel.inchikey.isin(set(friden.inchikey.dropna()))]
    if len(ov):
        p.append(f"去重失效：{len(ov)} 个分子仍同时出现在两个集合中："
                 + ", ".join(str(x) for x in ov.compound_name.head(8)))
    return p


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--b3db", type=Path, default=B3DB_REGRESSION)
    ap.add_argument("--friden", type=Path, default=FRIDEN_RAW)
    ap.add_argument("--gka", type=Path, default=GKA_CANDIDATES)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    started = datetime.now()
    for f in (args.b3db, args.friden, args.gka):
        if not f.exists():
            sys.exit(f"[FATAL] 找不到输入：{f}")

    print("[1/6] 读入并校验 …")
    check_sha(args.b3db)
    check_sha(args.friden)
    reg = pd.read_csv(args.b3db, sep="\t")
    gka = pd.read_csv(args.gka)
    print(f"      B3DB regression {len(reg)} 行；GKA 候选 {len(gka)} 行")

    print("[2/6] 解析结构、算骨架与描述符 …")
    reg = add_structure_columns(reg, "SMILES")
    print(f"      SMILES 可解析 {reg.parse_ok.sum()}/{len(reg)}；"
          f"骨架 {reg.murcko_scaffold.nunique()} 种")

    print("[3/6] 读入 Fridén 定量参考集 …")
    friden = load_friden(args.friden)
    friden = add_structure_columns(friden, "smiles")
    print(f"      {len(friden)} 个分子；有 Kp,uu 数值 {friden.kpuu_brain.notna().sum()}")

    print("[4/6] 选取 B3DB 分类对照（已排除与 Fridén 重复的分子）…")
    reg = select_b3db(reg, gka, set(friden.inchikey.dropna()))
    sel = reg[reg.selected].copy()
    print(f"      与 Fridén 重复而未选入 B3DB 侧的: "
          f"{int(reg.dropped_dup_with_friden.sum())}")
    print(f"      BBB+ {(sel.control_class=='control_positive').sum()}  "
          f"BBB- {(sel.control_class=='control_negative').sum()}")

    print("[5/6] 自检 …")
    problems = selfcheck(sel, friden)
    if problems:
        for x in problems:
            print(f"  ⚠ {x}")
    else:
        print("  ✓ 全部通过")

    print("[6/6] 写出 …")
    args.outdir.mkdir(parents=True, exist_ok=True)

    sel = sel.assign(control_set="b3db_classification", measure_basis="logBB")
    sel = sel.rename(columns={"SMILES": "smiles", "logBB": "logbb",
                              "BBB+/BBB-": "b3db_label", "group": "b3db_group",
                              "compound_name": "compound_name"})
    fr = friden.assign(control_set="friden_quantitative",
                       measure_basis="K(p,uu,brain)",
                       selected=True)

    keep_common = ["control_set", "control_class", "compound_name", "smiles",
                   "inchikey", "murcko_scaffold", "measure_basis",
                   *PROP_COLS, "gka_distance"]
    fr["control_class"] = ""
    fr["gka_distance"] = gka_proximity(fr, gka)
    combined = pd.concat(
        [sel.reindex(columns=keep_common + ["logbb", "b3db_label", "b3db_group",
                                            "selection_note"]),
         fr.reindex(columns=keep_common + ["kpuu_brain", "molecule_chembl_id",
                                           "pubmed_id", "doi"])],
        ignore_index=True,
    )
    combined.insert(0, "control_id",
                    [f"BBBCTRL_{i:03d}" for i in range(1, len(combined) + 1)])

    p_main = args.outdir / "Step3_00_BBB_Control_Set.csv"
    p_pool = args.outdir / "Step3_00_B3DB_Selection_Pool.csv"
    combined.to_csv(p_main, index=False)
    reg.to_csv(p_pool, index=False)
    print(f"      {p_main.name}  ({combined.shape[0]} × {combined.shape[1]})")
    print(f"      {p_pool.name}  ({reg.shape[0]} × {reg.shape[1]})  全池，加列不删行")

    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {p.name: check_sha(p) for p in (args.b3db, args.friden)},
        "gka_reference": str(args.gka),
        "rules": {
            "positive": f"regression group {POS_GROUP}, logBB >= {POS_LOGBB_MIN}",
            "negative": f"regression group {'->'.join(NEG_GROUPS)}, "
                        f"logBB <= {NEG_LOGBB_MAX}",
            "target_positive": TARGET_POS,
            "target_negative": TARGET_NEG,
            "scaffold_cap": SCAFFOLD_CAP,
            "weights": {"logbb_margin": W_MARGIN, "n_reference": W_NREF,
                        "gka_distance": W_GKA},
            "acyclic_as_one_cluster": ACYCLIC_AS_ONE_CLUSTER,
        },
        "counts": {
            "b3db_positive": int((sel.control_class == "control_positive").sum()),
            "b3db_negative": int((sel.control_class == "control_negative").sum()),
            "friden_quantitative": int(len(friden)),
            "total": int(len(combined)),
        },
        "selfcheck_problems": problems,
    }
    (args.outdir / "Step3_00_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n合计 {len(combined)} 个对照 "
          f"= {summary['counts']['b3db_positive']} BBB+ "
          f"+ {summary['counts']['b3db_negative']} BBB- "
          f"+ {summary['counts']['friden_quantitative']} Fridén")


if __name__ == "__main__":
    main()
