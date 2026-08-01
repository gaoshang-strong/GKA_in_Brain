#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step2_02_GKA_Patent_Tiering_and_Compound_Extraction
===================================================

把 Step2_01 的 GCK 相关专利分层，从选定层抽权利要求里的化合物，
给每个化合物算噪声指标 —— 产出 SureChEMBL 侧**独立**的 GKA 候选池。

**纯规则**：不调 API、不用 LLM。
**独立于 ChEMBL**：Step1 的分子不参与任何筛选，只在报告里做事后比较。

专利分层
--------
    L1  标题写 glucokinase activator      ← 唯一带方向信号的层
    L2  （标题提 GCK 或 权要提 GCK）且权要有化合物
    L3  其余                              ← 记录但不进后续

L2 要求「权要里有化合物」是为了排掉纯方法/用途权利要求（没有化学主体）。

⚠ 专利的层级不能继承给化合物
----------------------------
与 CLAUDE.md 里 Step1_04 那条同构：「assay 有目的性，但 assay 里的分子未必都符合」。
实测 `EP-4725482-A1`（GLUCOKINASE ACTIVATOR FOR COGNITIVE DISORDERS）权要里
只有 3 个化合物，**其中 2 个是葡萄糖**（开链式 + 环式）——权要写「激活葡萄糖激酶」
必然提到底物，抽取管道就把它注册成化合物了。

所以化合物层要独立判别。区分力最强的是 `n_global`（该结构在全库出现的专利数）：

    水  `O`      10,983,177 篇        特异性 0.00002
    多格列艾汀            187 篇        特异性 0.18

差 4 个数量级。指标**全部算成列、不直接删行**，改阈值不用重跑抽取。

输出
----
    Step2_02_GKA_Patent_Tiers.csv     专利分层表（含 L3）
    Step2_02_GKA_Compound_Pool.csv    主产物，一行一个化合物
    Step2_02_GKA_Patent_Tiering_and_Compound_Extraction.md   报告
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import duckdb
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_SNAP = REPO / "SureChEMBL" / "SureChEMBL_2026-07-17"
DEFAULT_STEP2_01 = (HERE.parent / "Step2_01_GCK_Related_Patent_Retrieval"
                    / "Step2_01_GCK_Related_Patents.csv")
DEFAULT_STEP1 = REPO / "Step1_Find_GKA_from_ChEMBL" / "Step1_GKA_Candidates_with_Properties.csv"

CLAIMS_FIELD_ID = 2          # fields.parquet: 2 = clms 权利要求

# L2 的第三支：说明书里反复讲 GCK 的化合物专利。
#
# ⚠ 初版 L2 只认 `by_title OR hit_clms>0`，漏掉了**典型的化合物专利**——
#   标题是纯化学名、权利要求里一次都不提靶点，靶点只写在说明书。实测被漏的：
#       Therapeutic agents                       权要提 GCK 0 次，权要 256 个化合物
#       2-(3,5-DISUBSTITUTEDPHENYL)PYRIMIDIN…    权要提 GCK 0 次，说明书提 37 次
#       NOVEL 2-PYRIDINECARBOXAMIDE DERIVATIVES  权要提 GCK 0 次，说明书提 32 次
#       Heteroaryl benzamide … as GLK activators 说明书提 7 次
#   这些 L3 专利的权要里装着 80 个 Step1 已知分子。
#
# 阈值 3 是量出来的：说明书提 1-2 次的有 23,787 篇（"顺带列举靶点"），
# >=3 次的只有 2,364 篇，信噪比陡变。
DESC_MENTION_MIN = 3

# 噪声过滤阈值 —— **每一条都用阳性对照标定过，误杀 0 个才留下**。
#
# ⚠ 初版拍了 5 个阈值没验对照，结果 12 个已知 GKA 只有 3 个通过。被误杀的两条：
#
#   specificity >= 0.01  误杀 2 个（AZD-1656 0.0063、Piraglitin 0.0072）
#       ——「越重要的药被后续专利引用越多，全库出现数越大，特异性反而越低」，
#         这个指标对成名已久的药有系统性偏见。
#   n_in_sel >= 2        误杀 5 个（MK-0941、PF-04991532、Neriglitin、LY-2608204、Piraglitin）
#       ——理由本身就错：一个药通常**只在它自己那一篇（族）专利的权要里**被具体画出来，
#         别的专利提它是当现有技术、出现在说明书而不是权要。
#         「只出现在 1 篇权要里」恰恰是原研化合物的特征，不是噪声的特征。
#
# 这两条已降级为标注列（`specificity` / `n_in_sel` 照常输出，不参与 keep 判定）。
# 本步骤是候选池，**宁可宽也不能把已知的药筛掉**，收窄留给后面有方向证据的步骤。
THRESH = {
    "n_global_max": 10_000,   # 砍水/溶剂/试剂/元素符号。对照最大值 316，留 30 倍余量
    "mw_min": 250.0,          # 砍元素、气体、小片段、葡萄糖(180)。对照最小 378.5
    "mw_max": 700.0,          # 砍聚合物。对照最大 559.8
}
# 只作标注、不参与 keep 判定的指标（保留在 CSV 里供下游自己收窄）
ANNOTATION_ONLY = ["specificity", "n_in_sel"]

PATENT_COLUMNS = [
    "patent_id", "patent_number", "country", "publication_date", "family_id",
    "title", "tier", "tier_reason", "title_says_activator", "by_title",
    "hit_clms", "hit_desc", "n_compounds_clms", "has_biomedical_annotation",
    "anchor_sources",
]

CMP_COLUMNS = [
    "compound_id", "inchi_key", "smiles", "mol_weight",
    # --- 噪声指标 ---
    "n_global", "n_in_sel", "n_in_L1", "n_in_L2", "specificity",
    "rdkit_valid", "has_carbon", "n_heavy_atoms",
    # --- 判定 ---
    "keep", "exclude_reason",
    # --- 溯源 ---
    "n_families", "example_patents",
    # --- 事后比较（不参与筛选）---
    "cmp_in_step1", "cmp_step1_chembl_ids",
]

T0 = time.time()


def log(m: str) -> None:
    print(f"[{time.time() - T0:6.1f}s] {m}", file=sys.stderr, flush=True)


def fmt(n) -> str:
    return "—" if n is None else f"{int(n):,}"


def pct(a, b) -> str:
    return "—" if not b else f"{100.0 * a / b:.1f}%"


def connect(snap: Path, threads: int, mem: str):
    con = duckdb.connect()
    con.execute(f"SET threads={threads}; SET memory_limit='{mem}';")
    con.execute("SET enable_progress_bar=false;")
    for t in ("patents", "compounds", "patent_compound_map"):
        p = snap / f"{t}.parquet"
        if not p.is_file():
            sys.exit(f"错误：找不到 {p}")
        con.execute(f"CREATE OR REPLACE VIEW {t} AS SELECT * FROM '{p}'")
    return con


def load_tiers(con, path: Path) -> dict:
    """读 Step2_01 主表并分层。"""
    if not path.is_file():
        sys.exit(f"错误：找不到 {path}，请先运行 Step2_01。")
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    con.execute("""CREATE OR REPLACE TABLE s1 (
        patent_id BIGINT, patent_number VARCHAR, country VARCHAR,
        publication_date VARCHAR, family_id BIGINT, title VARCHAR,
        title_says_activator VARCHAR, by_title VARCHAR, hit_clms INT,
        hit_desc INT, n_compounds_clms INT, has_biomedical_annotation VARCHAR,
        anchor_sources VARCHAR)""")
    con.executemany("INSERT INTO s1 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    [(int(r["patent_id"]), r["patent_number"], r["country"],
                      r["publication_date"], int(r["family_id"] or 0), r["title"],
                      r["title_says_activator"], r["by_title"],
                      int(r["hit_clms"] or 0), int(r["hit_desc"] or 0),
                      int(r["n_compounds_clms"] or 0),
                      r["has_biomedical_annotation"], r["anchor_sources"])
                     for r in rows])
    d = DESC_MENTION_MIN
    con.execute(f"""
        CREATE OR REPLACE TABLE tiers AS
        SELECT *,
          CASE WHEN title_says_activator = 'TRUE' THEN 'L1'
               WHEN (by_title = 'TRUE' OR hit_clms > 0 OR hit_desc >= {d})
                    AND n_compounds_clms > 0 THEN 'L2'
               ELSE 'L3' END AS tier,
          CASE WHEN title_says_activator = 'TRUE'
                 THEN '标题写 glucokinase activator（含方向信号）'
               WHEN (by_title = 'TRUE' OR hit_clms > 0) AND n_compounds_clms > 0
                 THEN '标题或权要提到 GCK，且权要有化合物'
               WHEN hit_desc >= {d} AND n_compounds_clms > 0
                 THEN '说明书反复讲 GCK（>= {d} 次）且权要有化合物 —— 典型化合物专利'
               WHEN (by_title = 'TRUE' OR hit_clms > 0 OR hit_desc >= {d})
                      AND n_compounds_clms = 0
                 THEN '提到 GCK 但权要没有化合物（纯方法/用途权利要求）'
               ELSE '仅顺带提及（说明书 < {d} 次）' END AS tier_reason
        FROM s1
    """)
    stat = {r[0]: (r[1], r[2]) for r in con.execute(
        "SELECT tier, COUNT(*), COUNT(DISTINCT family_id) FILTER (WHERE family_id>0) "
        "FROM tiers GROUP BY 1").fetchall()}
    for t in ("L1", "L2", "L3"):
        if t in stat:
            log(f"{t}: {stat[t][0]:,} 篇 / {stat[t][1]:,} 同族")
    return stat


def extract_compounds(con) -> int:
    """从 L1+L2 的权利要求抽化合物，并算噪声指标。"""
    con.execute(f"""
        CREATE OR REPLACE TABLE sel AS
        SELECT patent_id, family_id, tier FROM tiers WHERE tier IN ('L1', 'L2')
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE cand AS
        SELECT m.compound_id AS cid,
               COUNT(DISTINCT m.patent_id)                                        AS n_in_sel,
               COUNT(DISTINCT CASE WHEN s.tier='L1' THEN m.patent_id END)          AS n_in_L1,
               COUNT(DISTINCT CASE WHEN s.tier='L2' THEN m.patent_id END)          AS n_in_L2,
               COUNT(DISTINCT s.family_id) FILTER (WHERE s.family_id > 0)          AS n_families,
               list_sort(list_distinct(list(s.patent_id)))[1:3]                    AS ex_pid
        FROM patent_compound_map m
        JOIN sel s ON s.patent_id = m.patent_id
        WHERE m.field_id = {CLAIMS_FIELD_ID}
        GROUP BY 1
    """)
    n = con.execute("SELECT COUNT(*) FROM cand").fetchone()[0]
    log(f"L1+L2 权要化合物 {n:,} 个")

    # 全库出现的专利数 —— 噪声判别的核心指标
    con.execute("""
        CREATE OR REPLACE TABLE gfreq AS
        SELECT m.compound_id AS cid, COUNT(DISTINCT m.patent_id) AS n_global
        FROM patent_compound_map m JOIN cand ON cand.cid = m.compound_id
        GROUP BY 1
    """)
    log("全库出现专利数算完")
    return n


def enrich_and_judge(con) -> list:
    rows = con.execute("""
        SELECT c.cid, cp.inchi_key, cp.smiles, cp.mol_weight,
               g.n_global, c.n_in_sel, c.n_in_L1, c.n_in_L2, c.n_families, c.ex_pid
        FROM cand c
        JOIN gfreq g ON g.cid = c.cid
        JOIN compounds cp ON cp.id = c.cid
    """).fetchall()
    log(f"取回 {len(rows):,} 行，开始 RDKit 校验")

    out = []
    for cid, ik, smi, mw, ng, nsel, n1, n2, nfam, ex in rows:
        m = Chem.MolFromSmiles(smi) if smi else None
        valid = m is not None
        heavy = m.GetNumHeavyAtoms() if valid else 0
        carbon = bool(valid and any(a.GetSymbol() == "C" for a in m.GetAtoms()))
        spec = (nsel / ng) if ng else 0.0

        reasons = []
        if not valid:
            reasons.append("RDKit 无法解析")
        else:
            if heavy <= 1:
                reasons.append("单原子/空结构")
            if not carbon:
                reasons.append("不含碳")
        if ng > THRESH["n_global_max"]:
            reasons.append(f"全库出现 {ng:,} 篇 > {THRESH['n_global_max']:,}（通用物质）")
        if mw is None or not (THRESH["mw_min"] <= mw <= THRESH["mw_max"]):
            reasons.append(f"分子量 {mw:.1f} 不在 {THRESH['mw_min']:.0f}–{THRESH['mw_max']:.0f}"
                           if mw is not None else "分子量缺失")
        # specificity 与 n_in_sel 只作标注，**不参与 keep 判定**（见 THRESH 上方注释）

        out.append({
            "compound_id": cid, "inchi_key": ik, "smiles": smi,
            "mol_weight": round(mw, 3) if mw is not None else "",
            "n_global": ng, "n_in_sel": nsel, "n_in_L1": n1, "n_in_L2": n2,
            "specificity": round(spec, 6),
            "rdkit_valid": "TRUE" if valid else "FALSE",
            "has_carbon": "TRUE" if carbon else "FALSE",
            "n_heavy_atoms": heavy, "n_families": nfam,
            "example_patents": json.dumps(list(ex or []), ensure_ascii=False),
            "keep": "FALSE" if reasons else "TRUE",
            "exclude_reason": "；".join(reasons),
        })
    log(f"判定完成，keep = {sum(1 for r in out if r['keep'] == 'TRUE'):,}")
    return out


def compare_step1(con, cmps: list, path: Path) -> dict:
    """事后比较（不参与筛选）：与 ChEMBL 侧的重叠。"""
    if not path.is_file():
        log(f"警告：找不到 {path}，跳过事后比较")
        return {}
    with path.open(encoding="utf-8") as f:
        s1 = list(csv.DictReader(f))
    by_ik = {}
    for r in s1:
        ik = r.get("standard_inchi_key")
        if ik:
            by_ik.setdefault(ik, []).append(r["molecule_chembl_id"])
    for c in cmps:
        ids = by_ik.get(c["inchi_key"], [])
        c["cmp_in_step1"] = "TRUE" if ids else "FALSE"
        c["cmp_step1_chembl_ids"] = json.dumps(sorted(ids), ensure_ascii=False) if ids else ""
    kept = [c for c in cmps if c["keep"] == "TRUE"]
    ov = sum(1 for c in kept if c["cmp_in_step1"] == "TRUE")
    return {"n_step1": len(by_ik), "n_kept": len(kept), "n_overlap": ov,
            "n_new": len(kept) - ov}


def selfcheck_controls(cmps: list, step1_path: Path) -> dict:
    """阳性对照自检：已知 GKA 必须活着通过过滤。

    **这是硬性检查，不是可选项。** 初版没做，结果 12 个已知 GKA 只有 3 个通过
    （见 THRESH 上方注释）。任何阈值改动都必须先过这一关。
    """
    if not step1_path.is_file():
        return {}
    with step1_path.open(encoding="utf-8") as f:
        ctl = {r["standard_inchi_key"]: r["control_name"]
               for r in csv.DictReader(f)
               if r.get("is_positive_control") == "TRUE" and r.get("standard_inchi_key")}
    by_ik = {c["inchi_key"]: c for c in cmps if c["inchi_key"]}
    rows, killed = [], []
    for ik, name in sorted(ctl.items(), key=lambda kv: kv[1]):
        c = by_ik.get(ik)
        if c is None:
            rows.append({"name": name, "in_pool": False, "keep": None, "why": ""})
            continue
        ok = c["keep"] == "TRUE"
        rows.append({"name": name, "in_pool": True, "keep": ok,
                     "why": c["exclude_reason"], "n_global": c["n_global"],
                     "specificity": c["specificity"], "n_in_sel": c["n_in_sel"]})
        if not ok:
            killed.append((name, c["exclude_reason"]))
    in_pool = [r for r in rows if r["in_pool"]]
    n_ok = sum(1 for r in in_pool if r["keep"])
    log(f"阳性对照自检：池中 {len(in_pool)}/{len(ctl)}，通过过滤 {n_ok}/{len(in_pool)}")
    if killed:
        log("  ⚠⚠ 有已知 GKA 被过滤掉了，阈值必须放宽：")
        for n, w in killed:
            log(f"     {n}: {w}")
    return {"rows": rows, "n_ctl": len(ctl), "n_in_pool": len(in_pool),
            "n_pass": n_ok, "killed": killed}


def write_csv(rows: list, cols: list, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})


def write_report(con, cmps, tier_stat, cmp_stat, sc, path, s201, s1) -> None:
    L = []
    q = lambda s: con.execute(s).fetchall()          # noqa: E731
    one = lambda s: con.execute(s).fetchone()[0]     # noqa: E731
    kept = [c for c in cmps if c["keep"] == "TRUE"]

    L.append("# Step2_02 GKA 专利分层与化合物抽取")
    L.append("")
    L.append("- 数据源：**SureChEMBL 2.0，`2026-07-17` 全量快照**")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 输入：`{s201.name}`（Step2_01，{one('SELECT COUNT(*) FROM tiers'):,} 篇专利）")
    L.append(f"- 抽出化合物：**{len(cmps):,}** 个 → 通过噪声过滤 **{len(kept):,}** 个")
    L.append("")
    L.append("> **纯规则**，不调 API、不用 LLM。**不用 ChEMBL 做任何筛选**——"
             "Step1 的分子只在第六节做事后比较。")
    L.append("")

    L.append("## 一、专利分层")
    L.append("")
    rows = q("""SELECT tier, tier_reason, COUNT(*),
                       COUNT(DISTINCT family_id) FILTER (WHERE family_id>0)
                FROM tiers GROUP BY 1,2 ORDER BY 1, 3 DESC""")
    L.append("| 层 | 判据 | 专利 | 同族 |")
    L.append("| --- | --- | ---: | ---: |")
    for t, why, c, f in rows:
        L.append(f"| **{t}** | {why} | {fmt(c)} | {fmt(f)} |")
    L.append("")
    n_sel = sum(v[0] for k, v in tier_stat.items() if k in ("L1", "L2"))
    f_sel = one("SELECT COUNT(DISTINCT family_id) FROM tiers "
                "WHERE tier IN ('L1','L2') AND family_id > 0")
    f_all = one("SELECT COUNT(DISTINCT family_id) FROM tiers WHERE family_id > 0")
    L.append(f"**进入后续的是 L1 + L2：{fmt(n_sel)} 篇 / {fmt(f_sel)} 个同族**"
             f"（从 {fmt(f_all)} 个同族收窄到 {fmt(f_sel)}，砍掉 {pct(f_all - f_sel, f_all)}）。")
    L.append("")
    L.append("L1 是唯一带**方向信号**的层——标题直接写了是激活剂。"
             "L2 要求「权要里有化合物」，排掉纯方法/用途权利要求。"
             "L3 不进后续，但记录在 `Step2_02_GKA_Patent_Tiers.csv` 里，随时可回取。")
    L.append("")

    L.append("## 二、⚠ 专利的层级不能继承给化合物")
    L.append("")
    L.append("与 CLAUDE.md 里 Step1_04 那条同构——"
             "**「assay 有目的性，但 assay 里的分子未必都符合那个目的」**，"
             "换到专利就是：**一篇 GKA 专利的权利要求里，不是每个化合物都是 GKA**。")
    L.append("")
    L.append("实测 `EP-4725482-A1`（GLUCOKINASE ACTIVATOR FOR COGNITIVE DISORDERS）"
             "权要里只有 3 个化合物：")
    L.append("")
    L.append("| MW | SMILES | 是什么 |")
    L.append("| ---: | --- | --- |")
    L.append("| 180.156 | `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO` | **葡萄糖**（开链式） |")
    L.append("| 180.156 | `OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O` | **葡萄糖**（环式） |")
    L.append("| 462.927 | `CC(C)C[C@@H](C(=O)Nc1ccn(C[C@@H](O)CO)n1)N1CC(Oc2ccccc2Cl)…` | 多格列艾汀 |")
    L.append("")
    L.append("**3 个里 2 个是葡萄糖**——权要写「激活葡萄糖激酶」必然提到底物。"
             "所以化合物层必须独立判别。")
    L.append("")

    L.append("## 三、噪声指标的区分力")
    L.append("")
    L.append("`n_global`（该结构在**全库**出现的专利数）是区分力最强的指标。")
    L.append("")
    dirty = sorted(cmps, key=lambda c: -c["n_global"])[:6]
    clean = sorted([c for c in kept if c["n_global"] < 300],
                   key=lambda c: -c["n_in_L1"])[:6]
    L.append("| | `n_global` | 特异性 | MW | SMILES |")
    L.append("| --- | ---: | ---: | ---: | --- |")
    for c in dirty:
        L.append(f"| 最脏 | {fmt(c['n_global'])} | {c['specificity']:.5f} | "
                 f"{c['mol_weight']} | `{c['smiles'][:30]}` |")
    for c in clean:
        L.append(f"| 最干净 | {fmt(c['n_global'])} | {c['specificity']:.4f} | "
                 f"{c['mol_weight']} | `{c['smiles'][:30]}` |")
    L.append("")
    L.append("**特异性差 4 个数量级。** 分布：")
    L.append("")
    buckets = Counter()
    for c in cmps:
        g = c["n_global"]
        b = ("1" if g == 1 else "2–9" if g < 10 else "10–99" if g < 100
             else "100–999" if g < 1_000 else "1k–10k" if g < 10_000 else ">10k")
        buckets[b] += 1
    L.append("| `n_global` 区间 | 化合物数 |")
    L.append("| --- | ---: |")
    for b in ("1", "2–9", "10–99", "100–999", "1k–10k", ">10k"):
        if buckets.get(b):
            L.append(f"| {b} | {fmt(buckets[b])} |")
    L.append("")

    L.append("## 四、过滤结果")
    L.append("")
    L.append("阈值（可调，指标已存列，改阈值不用重跑抽取）：")
    L.append("")
    L.append("```")
    for k, v in THRESH.items():
        L.append(f"{k:<18} = {v}")
    L.append("```")
    L.append("")
    rc = Counter()
    for c in cmps:
        for r in filter(None, c["exclude_reason"].split("；")):
            rc[r.split("（")[0].split(" ")[0] if "分子量" not in r else "分子量不在窗口内"] += 1
    L.append("| 排除理由 | 化合物数 |")
    L.append("| --- | ---: |")
    for k, v in rc.most_common():
        L.append(f"| {k} | {fmt(v)} |")
    L.append("")
    L.append(f"**通过全部过滤：{fmt(len(kept))} / {fmt(len(cmps))}**"
             f"（{pct(len(kept), len(cmps))}）。"
             "被排除的行**保留在 CSV 里**，`keep = FALSE` 且写明 `exclude_reason`，可复核可反悔。")
    L.append("")
    n_l1 = sum(1 for c in kept if c["n_in_L1"] > 0)
    L.append(f"通过过滤的化合物里，**{fmt(n_l1)}** 个出现在 L1（标题写 activator）的专利中。")
    L.append("")

    if sc:
        L.append("## 五、阳性对照自检 ⭐")
        L.append("")
        L.append("**已知 GKA 必须活着通过过滤。这是硬性检查，任何阈值改动都要先过这一关。**")
        L.append("")
        L.append("| 对照 | 在池中 | 通过过滤 | `n_global` | 特异性 | `n_in_sel` |")
        L.append("| --- | :---: | :---: | ---: | ---: | ---: |")
        for r in sc["rows"]:
            if not r["in_pool"]:
                L.append(f"| {r['name']} | ❌ 不在 | — | — | — | — |")
            else:
                L.append(f"| {r['name']} | ✅ | {'✅' if r['keep'] else '**❌ ' + r['why'] + '**'} | "
                         f"{fmt(r['n_global'])} | {r['specificity']} | {r['n_in_sel']} |")
        L.append("")
        L.append(f"**结论：池中 {sc['n_in_pool']}/{sc['n_ctl']} 个对照，"
                 f"通过过滤 {sc['n_pass']}/{sc['n_in_pool']}。**")
        L.append("")
        if sc["killed"]:
            L.append("⚠⚠ **有已知 GKA 被过滤掉，阈值需要放宽。**")
        else:
            L.append("不在池中的是盐型（SureChEMBL 不单独注册盐）与专利未进 L1/L2 的，"
                     "属分层召回问题，与过滤规则无关。")
        L.append("")
        L.append("### 初版的教训")
        L.append("")
        L.append("初版拍了 5 个阈值**没验对照**，结果 12 个已知 GKA 只有 3 个通过。"
                 "被误杀的两条规则：")
        L.append("")
        L.append("| 规则 | 误杀 | 为什么错 |")
        L.append("| --- | ---: | --- |")
        L.append("| `specificity >= 0.01` | 2 | **越重要的药被后续专利引用越多**，"
                 "全库出现数越大、特异性反而越低（AZD-1656 0.0063、Piraglitin 0.0072）。"
                 "这个指标对成名已久的药有系统性偏见 |")
        L.append("| `n_in_sel >= 2` | 5 | 理由本身就错：一个药通常**只在它自己那一族专利的"
                 "权要里**被具体画出来，别处提它是当现有技术、出现在说明书。"
                 "**「只出现在 1 篇权要里」恰恰是原研化合物的特征** |")
        L.append("")
        L.append("两条已降级为标注列。这与 CLAUDE.md 里「**阈值要用阳性对照标定，"
                 "不能按整体分布拍**」是同一条——Step1_05 守住了，这里第一次没守住。")
        L.append("")

    L.append("## 六、⚠ 这套规则解决不了的")
    L.append("")
    L.append("| 问题 | 说明 |")
    L.append("| --- | --- |")
    L.append("| **特异性高 ≠ 是 GKA** | 合成中间体、砌块只出现在同一批专利里，"
             "特异性同样高。分子量下限挡掉大部分，**挡不干净** |")
    L.append("| **判不了「是不是权要主张的主体」** | 通式枚举物、中间体、参比化合物"
             "在 `patent_compound_map` 里没有区别，只能读权利要求原文 |")
    L.append("| **判不了方向** | 只有 L1 的标题信号，而且那是**专利级**不是化合物级 |")
    L.append("| **马库什枚举程度不均** | 有的专利枚举上万结构，有的只画通式，"
             "按化合物计数的统计会被大专利带偏 |")
    L.append("")

    if cmp_stat:
        L.append("## 七、事后比较：与 ChEMBL 侧的重叠")
        L.append("")
        L.append("> **以下不参与任何筛选。** SureChEMBL 侧是独立检索出来的，"
                 "这里只回答「两库各自找到了什么」。")
        L.append("")
        L.append("| | 数量 |")
        L.append("| --- | ---: |")
        L.append(f"| ChEMBL 侧（Step1 整合表） | {fmt(cmp_stat['n_step1'])} |")
        L.append(f"| SureChEMBL 侧（本表 `keep = TRUE`） | {fmt(cmp_stat['n_kept'])} |")
        L.append(f"| **两侧重叠** | {fmt(cmp_stat['n_overlap'])} |")
        L.append(f"| **SureChEMBL 独有（ChEMBL 里没有）** | {fmt(cmp_stat['n_new'])} |")
        L.append("")
        L.append(f"重叠率只有 {pct(cmp_stat['n_overlap'], cmp_stat['n_kept'])}——"
                 "**专利侧绝大多数化学实体在 ChEMBL 里查不到**，"
                 "这正是做专利挖掘的理由。但也要清醒：这批没有任何活性数据，"
                 "「是不是 GKA」还没被证实。")
        L.append("")

    L.append("## 八、下一步该做什么")
    L.append("")
    L.append("| 想解决 | 需要什么 |")
    L.append("| --- | --- |")
    L.append("| 化合物级的方向与主体判定 | 读权利要求原文（SureChEMBL API `/document/{id}/contents`） |")
    L.append("| 与 ChEMBL 侧合并成统一候选池 | 按 InChIKey 对齐，标明来源与证据强度 |")
    L.append("| 脑暴露筛选 | 算物化性质（本步骤只有 `mol_weight`） |")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="专利分层 + 权要化合物抽取（纯规则）。")
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAP)
    ap.add_argument("--step2-01", type=Path, default=DEFAULT_STEP2_01)
    ap.add_argument("--step1-csv", type=Path, default=DEFAULT_STEP1,
                    help="仅用于事后比较，不参与筛选")
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="12GB")
    args = ap.parse_args()

    con = connect(args.snapshot, args.threads, args.memory)
    log(f"快照 {args.snapshot}")
    tier_stat = load_tiers(con, args.step2_01)
    extract_compounds(con)
    cmps = enrich_and_judge(con)
    cmp_stat = compare_step1(con, cmps, args.step1_csv)
    sc = selfcheck_controls(cmps, args.step1_csv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    tier_csv = args.outdir / "Step2_02_GKA_Patent_Tiers.csv"
    pool_csv = args.outdir / "Step2_02_GKA_Compound_Pool.csv"
    md = args.outdir / "Step2_02_GKA_Patent_Tiering_and_Compound_Extraction.md"

    trows = con.execute("SELECT * FROM tiers ORDER BY tier, patent_number").fetchall()
    tnames = [d[0] for d in con.description]
    write_csv([dict(zip(tnames, r)) for r in trows], PATENT_COLUMNS, tier_csv)

    cmps.sort(key=lambda c: (c["keep"] != "TRUE", -c["n_in_L1"], -c["n_in_sel"],
                             c["compound_id"]))
    write_csv(cmps, CMP_COLUMNS, pool_csv)
    write_report(con, cmps, tier_stat, cmp_stat, sc, md, args.step2_01, args.step1_csv)

    kept = sum(1 for c in cmps if c["keep"] == "TRUE")
    print(f"\n专利分层表：{tier_csv}  （{len(trows):,} 篇）")
    print(f"化合物池：  {pool_csv}  （{len(cmps):,} 个，keep {kept:,}）")
    if cmp_stat:
        print(f"  与 ChEMBL 侧重叠 {cmp_stat['n_overlap']:,}，"
              f"SureChEMBL 独有 {cmp_stat['n_new']:,}")
    print(f"报告：      {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
