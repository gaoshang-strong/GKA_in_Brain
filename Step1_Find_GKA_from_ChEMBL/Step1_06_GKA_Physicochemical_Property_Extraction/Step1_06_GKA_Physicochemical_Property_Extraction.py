#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_06_GKA_Physicochemical_Property_Extraction
================================================

把 Step1_05 后续实验清单里的分子，理化性质全部从 ChEMBL 拉出来。

输入
----
    ../Step1_05_GCK_Activator_Ranking_and_Candidate_Selection/
        Step1_05_Followup_Candidates.csv          782 个分子

路径
----
    候选 chembl_id → molecule_dictionary.molregno
                   → compound_structures    SMILES / InChI / InChIKey
                   → compound_properties    MW / ALogP / PSA / HBD / HBA / RTB / QED / Ro5
                   → molecule_hierarchy     盐型 → 母体
                   → molecule_synonyms      名称、研发代号

**本步骤只取数，不做筛选、不做判断。** 性质窗口（能不能进脑）留给 Step2。

ChEMBL 37 缺哪些性质
--------------------
这一版的 `compound_properties` 只有 15 列，**没有** `cx_logp` / `cx_logd` /
`cx_most_apka` / `cx_most_bpka` / `molecular_species`。缺的恰好是判断脑暴露最关键的：

- **`cx_logd`（pH 7.4 分配系数）**：GKA 里有相当一部分是羧酸，在生理 pH 下带负电，
  logD 会比 logP 低几个数量级。**只看 `alogp` 会高估这批分子的膜通透性。**
- **`molecular_species`（酸/碱/中性）**：酸性化合物入脑普遍差，这是 CNS 项目的第一刀。

这两项 Step2 必须自己算（RDKit 能算 logP 与酸性基团，pKa/logD 需另找工具）。
本步骤把 `alogp` 如实取出并标注它**不是 logD**，不代替。

注解字段基本是空的
------------------
`chirality` / `prodrug` / `first_in_class` / `inorganic_flag` 在 777/782 上是 `-1`
（未标注），`oral` / `parenteral` / `therapeutic_flag` 等全是 0。
按项目约定**空值也是事实，如实记录不省略字段**——这批分子绝大多数是文献化合物，
没进过开发流程，本来就不会有这些注解。

输出
----
    Step1_06_GKA_Physicochemical_Properties.csv        一行一个分子（主产物）
    Step1_06_GKA_Physicochemical_Property_Extraction.md 报告

用法
----
    python3 Step1_06_GKA_Physicochemical_Property_Extraction.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_DB = REPO / "ChEMBL" / "ChEMBL_37" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"
DEFAULT_IN = (HERE.parent / "Step1_05_GCK_Activator_Ranking_and_Candidate_Selection"
              / "Step1_05_Followup_Candidates.csv")

# 从 Step1_05 带过来的链路列，保证这张表能直接接着用
CHAIN_COLUMNS = [
    "priority", "is_positive_control", "control_name",
    "potency_band", "pactivity_median", "potency_nm_min", "direction",
    "murcko_scaffold", "scaffold_cluster_size", "is_scaffold_representative",
]

PROP_COLUMNS = [
    "mw_freebase", "full_mwt", "alogp", "hba", "hbd", "psa", "rtb",
    "ro3_pass", "num_ro5_violations", "aromatic_rings", "heavy_atoms",
    "qed_weighted", "full_molformula", "np_likeness_score",
]

OUT_COLUMNS = (
    # --- 身份 ---
    ["molecule_chembl_id", "molregno", "molecule_pref_name", "molecule_type",
     "structure_type", "max_phase", "first_approval", "availability_type",
     "chirality", "prodrug", "natural_product", "therapeutic_flag",
     "first_in_class", "inorganic_flag", "chemical_probe", "orphan",
     "usan_stem", "usan_year"]
    # --- 母体/盐型 ---
    + ["parent_molregno", "parent_chembl_id", "is_parent"]
    # --- 结构 ---
    + ["canonical_smiles", "standard_inchi_key", "standard_inchi"]
    # --- 理化性质 ---
    + PROP_COLUMNS
    # --- 同义词 ---
    + ["n_synonyms", "synonyms"]
    # --- 链路 ---
    + CHAIN_COLUMNS
)

MOL_SQL = """
SELECT md.chembl_id, md.molregno, md.pref_name, md.molecule_type,
       md.structure_type, md.max_phase, md.first_approval, md.availability_type,
       md.chirality, md.prodrug, md.natural_product, md.therapeutic_flag,
       md.first_in_class, md.inorganic_flag, md.chemical_probe, md.orphan,
       md.usan_stem, md.usan_year,
       mh.parent_molregno, pmd.chembl_id AS parent_chembl_id,
       cs.canonical_smiles, cs.standard_inchi_key, cs.standard_inchi,
       cp.mw_freebase, cp.full_mwt, cp.alogp, cp.hba, cp.hbd, cp.psa, cp.rtb,
       cp.ro3_pass, cp.num_ro5_violations, cp.aromatic_rings, cp.heavy_atoms,
       cp.qed_weighted, cp.full_molformula, cp.np_likeness_score,
       cp.molregno AS has_props
FROM molecule_dictionary md
LEFT JOIN molecule_hierarchy   mh  ON mh.molregno = md.molregno
LEFT JOIN molecule_dictionary  pmd ON pmd.molregno = mh.parent_molregno
LEFT JOIN compound_structures  cs  ON cs.molregno = md.molregno
LEFT JOIN compound_properties  cp  ON cp.molregno = md.molregno
WHERE md.chembl_id IN ({placeholders})
"""

SYN_SQL = """
SELECT md.chembl_id, ms.syn_type, ms.synonyms
FROM molecule_synonyms ms
JOIN molecule_dictionary md ON md.molregno = ms.molregno
WHERE md.chembl_id IN ({placeholders})
ORDER BY md.chembl_id, ms.syn_type, ms.synonyms
"""


def connect_readonly(db_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.text_factory = lambda b: b.decode("utf-8", "replace")
    return con


def chembl_version(con: sqlite3.Connection) -> str:
    row = con.execute(
        "SELECT chembl_release, creation_date FROM chembl_release "
        "ORDER BY chembl_release_id DESC LIMIT 1"
    ).fetchone()
    return f"{row[0]}（{str(row[1])[:10]}）" if row else "未知"


def fetch(con: sqlite3.Connection, cands: list) -> tuple:
    ids = [c["molecule_chembl_id"] for c in cands]
    ph = ",".join("?" * len(ids))

    syns = defaultdict(list)
    for r in con.execute(SYN_SQL.format(placeholders=ph), ids):
        syns[r["chembl_id"]].append({"type": r["syn_type"], "name": r["synonyms"]})

    by_id = {}
    for r in con.execute(MOL_SQL.format(placeholders=ph), ids):
        d = dict(r)
        d["molecule_chembl_id"] = d.pop("chembl_id")
        d["molecule_pref_name"] = d.pop("pref_name")
        d["is_parent"] = ("TRUE" if (d["parent_molregno"] is None
                                     or d["parent_molregno"] == d["molregno"])
                          else "FALSE")
        s = syns.get(d["molecule_chembl_id"], [])
        d["n_synonyms"] = len(s)
        d["synonyms"] = json.dumps(s, ensure_ascii=False) if s else ""
        by_id[d["molecule_chembl_id"]] = d

    rows, missing = [], []
    for c in cands:
        mid = c["molecule_chembl_id"]
        d = by_id.get(mid)
        if d is None:
            missing.append(mid)
            continue
        for k in CHAIN_COLUMNS:
            d[k] = c.get(k, "")
        rows.append(d)
    return rows, missing


def write_csv(rows: list, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in OUT_COLUMNS})


def num(rows: list, field: str) -> list:
    out = []
    for r in rows:
        v = r.get(field)
        if v is None or v == "":
            continue
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            pass
    return out


def write_report(rows, missing, path, in_csv, db, version) -> None:
    L = []
    n = len(rows)
    L.append("# Step1_06 GKA 候选分子理化性质提取")
    L.append("")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db}`")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 输入：`{in_csv.name}`（Step1_05 后续实验清单）")
    L.append(f"- 提取分子：**{n:,}** 个"
             + (f"，**{len(missing)}** 个在 ChEMBL 里查不到" if missing else "，全部命中"))
    L.append("")
    L.append("> **本步骤只取数，不做筛选、不做判断。** 性质窗口（能不能进脑）留给 Step2。")
    L.append("")

    # --- 覆盖 ---
    L.append("## 一、字段覆盖")
    L.append("")
    L.append("| 字段 | 有值 | 缺失 | 说明 |")
    L.append("| --- | ---: | ---: | --- |")
    desc = {
        "canonical_smiles": "结构",
        "standard_inchi_key": "结构哈希，跨库对齐用",
        "mw_freebase": "游离碱分子量",
        "full_mwt": "含盐分子量",
        "alogp": "计算 logP，**不是 logD**",
        "psa": "极性表面积（TPSA）",
        "hbd": "氢键供体",
        "hba": "氢键受体",
        "rtb": "可旋转键",
        "aromatic_rings": "芳环数",
        "heavy_atoms": "重原子数",
        "qed_weighted": "类药性 QED (0–1)",
        "num_ro5_violations": "Lipinski 五规则违规数",
        "ro3_pass": "是否通过 Rule of 3（片段筛选用）",
        "np_likeness_score": "天然产物相似度",
        "full_molformula": "分子式",
    }
    for f, d in desc.items():
        have = sum(1 for r in rows if r.get(f) not in (None, ""))
        L.append(f"| `{f}` | {have:,} | {n - have:,} | {d} |")
    L.append("")

    # --- 缺什么 ---
    L.append("## 二、ChEMBL 37 缺的性质（Step2 必须自己算）")
    L.append("")
    L.append("这一版的 `compound_properties` 只有 15 列，**没有**下面这些。"
             "缺的恰好是判断脑暴露最关键的：")
    L.append("")
    L.append("| 缺失字段 | 是什么 | 为什么这里要命 |")
    L.append("| --- | --- | --- |")
    L.append("| `cx_logd` | pH 7.4 下的分配系数 | GKA 里有相当一部分是羧酸，"
             "生理 pH 下带负电，logD 比 logP 低几个数量级。**只看 `alogp` 会高估膜通透性** |")
    L.append("| `molecular_species` | 酸 / 碱 / 中性 | 酸性化合物入脑普遍差，"
             "这是 CNS 项目的第一刀 |")
    L.append("| `cx_most_apka` / `cx_most_bpka` | 最强酸/碱解离常数 | 没有 pKa 就算不出 logD |")
    L.append("")
    acid = sum(1 for r in rows
               if r.get("canonical_smiles")
               and ("C(=O)O" in r["canonical_smiles"] or "C(O)=O" in r["canonical_smiles"]))
    L.append(f"粗查 SMILES，**{acid:,} / {n:,}** 个分子含羧酸基团"
             "（只是字符串匹配，Step2 应改用 SMARTS 正式判定）。"
             "本步骤把 `alogp` 如实取出并标注它**不是 logD**，不代替。")
    L.append("")

    # --- 分布 ---
    L.append("## 三、性质分布")
    L.append("")
    L.append("**描述性统计，不构成筛选。** 阈值列只是让你先看到量级。")
    L.append("")
    L.append("| 性质 | 中位 | 最小 | 最大 | 常用 CNS 窗口 | 落在窗口内 |")
    L.append("| --- | ---: | ---: | ---: | --- | ---: |")
    windows = [
        ("mw_freebase", "MW", "≤ 450", lambda v: v <= 450),
        ("alogp", "ALogP", "1 – 5", lambda v: 1 <= v <= 5),
        ("psa", "PSA", "≤ 90", lambda v: v <= 90),
        ("hbd", "HBD", "≤ 2", lambda v: v <= 2),
        ("hba", "HBA", "≤ 7", lambda v: v <= 7),
        ("rtb", "RotB", "≤ 8", lambda v: v <= 8),
        ("aromatic_rings", "芳环数", "≤ 3", lambda v: v <= 3),
        ("qed_weighted", "QED", "≥ 0.5", lambda v: v >= 0.5),
    ]
    for f, label, win, ok in windows:
        v = num(rows, f)
        if not v:
            L.append(f"| {label} (`{f}`) | — | — | — | {win} | — |")
            continue
        L.append(f"| {label} (`{f}`) | {statistics.median(v):.2f} | {min(v):.2f} | "
                 f"{max(v):.2f} | {win} | {sum(1 for x in v if ok(x)):,} |")
    L.append("")
    ro5 = Counter(str(r.get("num_ro5_violations", "")) for r in rows)
    L.append("| Lipinski 违规数 | 分子数 |")
    L.append("| ---: | ---: |")
    for k in sorted(ro5, key=lambda x: (x == "", x)):
        L.append(f"| {k or '(空)'} | {ro5[k]:,} |")
    L.append("")
    all_ok = sum(1 for r in rows if all(
        (lambda v: v is not None and ok(v))(
            float(r[f]) if r.get(f) not in (None, "") else None)
        for f, _, _, ok in windows[:6]))
    L.append(f"六项一起看（MW / ALogP / PSA / HBD / HBA / RotB 全部落在窗口内）："
             f"**{all_ok:,} / {n:,}**。")
    L.append("")
    L.append("> 这批分子是冲着**肝和胰腺**做的，不是冲着脑做的——"
             "PSA 中位数已经超过常用 CNS 上限，是意料之中的结果。"
             "**这不是筛选结论**，Step2 会用正式的判据重做。")
    L.append("")

    # --- 注解 ---
    L.append("## 四、注解字段（大多为空，如实记录）")
    L.append("")
    L.append("`-1` 是 ChEMBL 的「未标注」，不是 0。这批分子绝大多数是文献化合物，"
             "没进过开发流程，本来就不会有这些注解。")
    L.append("")
    L.append("| 字段 | 取值分布 |")
    L.append("| --- | --- |")
    for f in ("max_phase", "chirality", "prodrug", "natural_product",
              "first_in_class", "inorganic_flag", "therapeutic_flag",
              "availability_type", "chemical_probe", "orphan"):
        c = Counter(str(r.get(f)) if r.get(f) not in (None, "") else "(空)" for r in rows)
        L.append(f"| `{f}` | " + "、".join(f"`{k}`×{v:,}" for k, v in c.most_common(5)) + " |")
    L.append("")

    # --- 盐型与同义词 ---
    L.append("## 五、盐型归并与同义词")
    L.append("")
    n_salt = sum(1 for r in rows if r.get("is_parent") == "FALSE")
    L.append(f"`molecule_hierarchy` 里 **{n - n_salt:,}** 个分子本身就是母体，"
             f"**{n_salt:,}** 个是盐型/衍生记录（`is_parent = FALSE`，"
             "`parent_chembl_id` 给出母体）。")
    if n_salt:
        L.append("")
        L.append("| 分子 | 母体 |")
        L.append("| --- | --- |")
        for r in rows:
            if r.get("is_parent") == "FALSE":
                L.append(f"| `{r['molecule_chembl_id']}` | "
                         f"`{r.get('parent_chembl_id') or '—'}` |")
    L.append("")
    n_syn = sum(1 for r in rows if r.get("n_synonyms"))
    L.append(f"`molecule_synonyms` 只覆盖 **{n_syn:,} / {n:,}** 个分子——"
             "有名字的基本就是进过临床的那几个，其余是文献化合物只有 ChEMBL ID。")
    if n_syn:
        L.append("")
        L.append("| 分子 | 名称 / 研发代号 |")
        L.append("| --- | --- |")
        for r in rows:
            if r.get("n_synonyms"):
                names = "、".join(f"{s['name']}（{s['type']}）"
                                 for s in json.loads(r["synonyms"]))
                L.append(f"| `{r['molecule_chembl_id']}` | {names[:120]} |")
    L.append("")

    if missing:
        L.append("## 六、未命中的分子")
        L.append("")
        L.append("、".join(f"`{m}`" for m in missing))
        L.append("")

    L.append("## 附：按 priority 的性质概览")
    L.append("")
    L.append("| priority | 分子数 | MW 中位 | ALogP 中位 | PSA 中位 | 六项全落窗口 |")
    L.append("| --- | ---: | ---: | ---: | ---: | ---: |")
    for pri in ("P1", "P2", "P3", "P4"):
        g = [r for r in rows if r.get("priority") == pri]
        if not g:
            continue
        ok_n = sum(1 for r in g if all(
            (lambda v: v is not None and ok(v))(
                float(r[f]) if r.get(f) not in (None, "") else None)
            for f, _, _, ok in windows[:6]))
        L.append(f"| {pri} | {len(g):,} | "
                 f"{statistics.median(num(g, 'mw_freebase')):.1f} | "
                 f"{statistics.median(num(g, 'alogp')):.2f} | "
                 f"{statistics.median(num(g, 'psa')):.1f} | {ok_n:,} |")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(
        description="把 Step1_05 候选清单的分子理化性质从 ChEMBL 全部取出。")
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--in-csv", type=Path, default=DEFAULT_IN)
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    if not args.db.is_file():
        print(f"错误：找不到数据库 {args.db}", file=sys.stderr)
        return 1
    if not args.in_csv.is_file():
        print(f"错误：找不到 {args.in_csv}，请先运行 Step1_05。", file=sys.stderr)
        return 1

    with args.in_csv.open(encoding="utf-8") as f:
        cands = list(csv.DictReader(f))

    con = connect_readonly(args.db)
    version = chembl_version(con)
    print(f"数据库：{args.db}")
    print(f"ChEMBL 版本：{version}")
    print(f"输入：{args.in_csv.name}（{len(cands):,} 个候选）\n")

    rows, missing = fetch(con, cands)
    print(f"命中 {len(rows):,} / {len(cands):,}"
          + (f"，未命中 {len(missing)}：{missing[:5]}" if missing else "，全部命中"))
    have_props = sum(1 for r in rows if r.get("has_props") is not None)
    print(f"  有 compound_properties：{have_props:,}")
    print(f"  有 SMILES：            {sum(1 for r in rows if r.get('canonical_smiles')):,}")
    print(f"  非母体（盐型等）：      {sum(1 for r in rows if r.get('is_parent') == 'FALSE'):,}")
    print(f"  有同义词：              {sum(1 for r in rows if r.get('n_synonyms')):,}")

    rows.sort(key=lambda r: (r.get("priority") or "ZZ", r["molecule_chembl_id"]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = args.outdir / "Step1_06_GKA_Physicochemical_Properties.csv"
    out_md = args.outdir / "Step1_06_GKA_Physicochemical_Property_Extraction.md"
    write_csv(rows, out_csv)
    write_report(rows, missing, out_md, args.in_csv, args.db, version)

    print(f"\n性质主表：{out_csv}")
    print(f"报告：    {out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
