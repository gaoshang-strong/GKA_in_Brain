#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_01_GCK_Target_Anchoring
=============================

以 UniProt accession 为锚点，定位 ChEMBL 中所有对应该蛋白组件的 target 记录。
默认锚点为 P35557（人葡萄糖激酶 GCK / ChEMBL 中的 pref_name 为 Hexokinase-4）。

为什么用 accession 而不是名字锚定
--------------------------------
GCK 在 ChEMBL 里的 pref_name 是 "Hexokinase-4"，按 "Glucokinase" 搜名字会漏掉主靶点，
而按 "Glucokinase" 搜到的 CHEMBL1075152 其实是 GKRP（调节蛋白，另一个基因 GCKR）。
UniProt accession 是 ChEMBL 与蛋白世界之间唯一稳定的桥梁，因此锚定必须走 accession。

查询链路
--------
    accession → component_sequences → component_id
              → target_components   → tid
              → target_dictionary   → chembl_id

一个蛋白组件可以挂在多个 target 上（单一蛋白靶点、蛋白复合物、PPI 靶点…），
所以输出是一张多行的映射表，不是一行。

输出
----
    Step1_01_GCK_Target_Anchoring.csv   映射表（主产物）
    Step1_01_GCK_Target_Anchoring.md    带出处信息的可读报告

用法
----
    python3 Step1_01_GCK_Target_Anchoring.py
    python3 Step1_01_GCK_Target_Anchoring.py --accession P35557 --db /path/to/chembl_37.db
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                      # /ShangGaoAIProjects/GKA_in_Brain
DEFAULT_DB = REPO / "ChEMBL" / "ChEMBL_37" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"

# 输出列：README 要求的 7 列在最前，其余为下游步骤有用的补充信息
COLUMNS = [
    # --- README 要求 ---
    "accession",
    "component_id",
    "tid",
    "target_chembl_id",
    "pref_name",
    "organism",
    "target_type",
    # --- 组件侧补充 ---
    "component_type",
    "component_description",
    "component_organism",
    "component_tax_id",
    "component_db_source",
    "component_db_version",
    # --- 映射关系补充 ---
    "homologue",
    # --- 靶点侧补充 ---
    "target_tax_id",
    "species_group_flag",
    # --- 数据量，供下一步判断该靶点值不值得深挖 ---
    "n_assays",
    "n_assays_conf_ge8",
    "n_activities",
    "n_activities_with_pchembl",
]

ANCHOR_SQL = """
SELECT
    cs.accession                AS accession,
    cs.component_id             AS component_id,
    td.tid                      AS tid,
    td.chembl_id                AS target_chembl_id,
    td.pref_name                AS pref_name,
    td.organism                 AS organism,
    td.target_type              AS target_type,
    cs.component_type           AS component_type,
    cs.description              AS component_description,
    cs.organism                 AS component_organism,
    cs.tax_id                   AS component_tax_id,
    cs.db_source                AS component_db_source,
    cs.db_version               AS component_db_version,
    tc.homologue                AS homologue,
    td.tax_id                   AS target_tax_id,
    td.species_group_flag       AS species_group_flag
FROM component_sequences cs
JOIN target_components  tc ON tc.component_id = cs.component_id
JOIN target_dictionary  td ON td.tid          = tc.tid
WHERE cs.accession = ?
ORDER BY td.target_type, td.chembl_id
"""

# 每个 target 的数据量。分开查是有意的：把锚定逻辑和计数逻辑解耦，
# 前者必须精确，后者只是参考信息。
COUNT_SQL = """
SELECT
    (SELECT COUNT(*) FROM assays WHERE tid = :tid),
    (SELECT COUNT(*) FROM assays WHERE tid = :tid AND confidence_score >= 8),
    (SELECT COUNT(*) FROM activities act JOIN assays a ON act.assay_id = a.assay_id
      WHERE a.tid = :tid),
    (SELECT COUNT(*) FROM activities act JOIN assays a ON act.assay_id = a.assay_id
      WHERE a.tid = :tid AND act.pchembl_value IS NOT NULL)
"""


def connect_readonly(db_path: Path) -> sqlite3.Connection:
    """只读打开，杜绝脚本意外写入这个 30 GB 的库。"""
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.text_factory = lambda b: b.decode("utf-8", "replace")
    return con


def chembl_version(con: sqlite3.Connection) -> str:
    try:
        row = con.execute(
            "SELECT chembl_release, creation_date FROM chembl_release "
            "ORDER BY chembl_release_id DESC LIMIT 1"
        ).fetchone()
        if row:
            return f"{row[0]}（{str(row[1])[:10]}）"
    except sqlite3.Error:
        pass
    return "未知"


def anchor(con: sqlite3.Connection, accession: str) -> list[dict]:
    rows = [dict(r) for r in con.execute(ANCHOR_SQL, (accession,))]
    for row in rows:
        counts = con.execute(COUNT_SQL, {"tid": row["tid"]}).fetchone()
        (row["n_assays"], row["n_assays_conf_ge8"],
         row["n_activities"], row["n_activities_with_pchembl"]) = counts
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in COLUMNS})


def write_report(rows: list[dict], path: Path, accession: str,
                 db_path: Path, version: str) -> None:
    L: list[str] = []
    L.append(f"# Step1_01 GCK 靶点锚定结果")
    L.append("")
    L.append(f"- 锚点 UniProt accession：**{accession}**")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db_path}`")
    L.append(f"- 运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- 命中 target 数：**{len(rows)}**")
    L.append("")

    if rows:
        comp = rows[0]
        L.append("## 蛋白组件")
        L.append("")
        L.append(f"- `component_id` = **{comp['component_id']}**"
                 f"（{comp['component_type']}）")
        L.append(f"- 描述：{comp['component_description']}")
        L.append(f"- 物种：{comp['component_organism']}（tax_id {comp['component_tax_id']}）")
        L.append(f"- 序列来源：{comp['component_db_source']} {comp['component_db_version']}")
        L.append("")
        L.append("> 注：这里的版本取自 `component_sequences.db_version`（逐条记录）。"
                 "ChEMBL 的 `version` 表另有一条全局声明，两者在 ChEMBL 37 中并不一致"
                 "（全局声明为 Swiss-Prot 2025_03，而绝大多数组件记录标注为 2026_01）。"
                 "以逐条记录的值为准。")
        L.append("")

    L.append("## 命中的 ChEMBL target")
    L.append("")
    L.append("| target_chembl_id | tid | pref_name | organism | target_type | homologue | "
             "assays | assays(conf≥8) | activities | 有 pChEMBL |")
    L.append("| --- | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |")
    for r in rows:
        L.append(
            f"| `{r['target_chembl_id']}` | {r['tid']} | {r['pref_name']} | {r['organism']} | "
            f"{r['target_type']} | {r['homologue']} | {r['n_assays']:,} | "
            f"{r['n_assays_conf_ge8']:,} | {r['n_activities']:,} | "
            f"{r['n_activities_with_pchembl']:,} |"
        )
    L.append("")
    L.append("> `homologue = 0` 表示该组件是这个 target 的直接组成部分，而非同源映射。")
    L.append("")
    L.append("## 说明")
    L.append("")
    L.append("同一个蛋白组件会挂在多个 target 上：既有把它单独作为作用对象的 "
             "`SINGLE PROTEIN` 靶点，也有把它与其他蛋白的相互作用作为作用对象的 "
             "`PROTEIN-PROTEIN INTERACTION` 靶点。下游取活性数据时要**显式决定保留哪些**，"
             "不能默认只有一个。")
    L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(
        description="以 UniProt accession 锚定 ChEMBL 中的 target 记录。")
    p.add_argument("--accession", default="P35557",
                   help="UniProt accession，默认 P35557（人 GCK）")
    p.add_argument("--db", type=Path, default=DEFAULT_DB,
                   help=f"ChEMBL SQLite 路径，默认 {DEFAULT_DB}")
    p.add_argument("--outdir", type=Path, default=HERE,
                   help="输出目录，默认与本脚本同目录")
    args = p.parse_args()

    if not args.db.is_file():
        print(f"错误：找不到数据库 {args.db}", file=sys.stderr)
        print("      ChEMBL 数据库不入 git，请先按 ChEMBL/ChEMBL_37/README_ChEMBL_37.md 下载。",
              file=sys.stderr)
        return 1

    con = connect_readonly(args.db)
    version = chembl_version(con)
    print(f"数据库：{args.db}")
    print(f"ChEMBL 版本：{version}")
    print(f"锚点：{args.accession}\n")

    rows = anchor(con, args.accession)
    if not rows:
        print(f"错误：accession {args.accession} 在 component_sequences 中无匹配，"
              "或它没有关联到任何 target。", file=sys.stderr)
        return 2

    # 控制台输出核心几列，完整结果看 CSV
    w = [max(len(str(r[c])) for r in rows + [{c: c}]) for c in
         ("target_chembl_id", "pref_name", "organism", "target_type")]
    hdr = ("target_chembl_id", "pref_name", "organism", "target_type")
    print("  ".join(h.ljust(x) for h, x in zip(hdr, w)) + "  activities")
    print("  ".join("-" * x for x in w) + "  ----------")
    for r in rows:
        print("  ".join(str(r[c]).ljust(x) for c, x in zip(hdr, w))
              + f"  {r['n_activities']:,}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "Step1_01_GCK_Target_Anchoring.csv"
    md_path = args.outdir / "Step1_01_GCK_Target_Anchoring.md"
    write_csv(rows, csv_path)
    write_report(rows, md_path, args.accession, args.db, version)

    print(f"\n命中 {len(rows)} 个 target")
    print(f"映射表：{csv_path}")
    print(f"报告：  {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
