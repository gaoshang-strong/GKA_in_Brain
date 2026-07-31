#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_02_GCK_Assay_Mapping
==========================

由 Step1_01 锚定到的 GCK 靶点出发，提取全部相关 assay，汇总每个 assay 的
实验类型、描述、物种、靶点可信度、来源及活性数据规模。

查询路径
--------
    target_dictionary.tid → assays.tid → assay_id / chembl_id / description

本步骤**只做清点，不做判定**：不推断某个 assay 测的是激活还是抑制。
输出表提供做这个判断所需的全部原始信息（描述原文、实测的 standard_type
分布、BAO 实验格式、文献出处），判断由人来做。

输入
----
默认读取上一步的输出 `../Step1_01_GCK_Target_Anchoring/Step1_01_GCK_Target_Anchoring.csv`
取其中的 tid；也可用 --tid 直接指定。

输出
----
    Step1_02_GCK_Assay_Mapping.csv   assay 清单（主产物，一行一个 assay）
    Step1_02_GCK_Assay_Mapping.md    汇总报告（分布统计 + 完整清单）

用法
----
    python3 Step1_02_GCK_Assay_Mapping.py
    python3 Step1_02_GCK_Assay_Mapping.py --tid 20095
    python3 Step1_02_GCK_Assay_Mapping.py --db /path/to/chembl_37.db
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_DB = REPO / "ChEMBL" / "ChEMBL_37" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"
DEFAULT_ANCHOR = HERE.parent / "Step1_01_GCK_Target_Anchoring" / "Step1_01_GCK_Target_Anchoring.csv"

COLUMNS = [
    # --- 靶点 ---
    "tid", "target_chembl_id", "target_pref_name", "target_type",
    # --- assay 身份与描述 ---
    "assay_id", "assay_chembl_id", "assay_description",
    # --- 实验类型 ---
    "assay_type", "assay_type_desc", "assay_test_type", "assay_category",
    "bao_format", "bao_label",
    # --- 物种与实验体系 ---
    "assay_organism", "assay_tax_id", "assay_strain",
    "assay_tissue", "assay_cell_type", "assay_subcellular_fraction",
    # --- 突变体 ---
    "variant_id", "variant_mutation", "variant_accession", "variant_isoform",
    # --- 靶点可信度 ---
    "confidence_score", "confidence_desc", "relationship_type", "relationship_desc",
    # --- 来源与出处 ---
    "src_id", "src_short_name", "src_assay_id",
    "doc_chembl_id", "doc_year", "doc_journal", "pubmed_id", "doi", "doc_title",
    # --- 聚合字段（JSON 列表，保持一个 assay 一行）---
    "assay_parameters", "n_assay_parameters",
    "assay_classifications", "n_assay_classifications",
    # --- 活性数据规模 ---
    "n_activities", "n_activities_with_pchembl", "n_distinct_compounds",
    "standard_types", "activity_comments", "n_activity_comments",
]

ASSAY_SQL = """
SELECT
    a.assay_id                      AS assay_id,
    a.chembl_id                     AS assay_chembl_id,
    a.description                   AS assay_description,
    a.assay_type                    AS assay_type,
    at.assay_desc                   AS assay_type_desc,
    a.assay_test_type               AS assay_test_type,
    a.assay_category                AS assay_category,
    a.bao_format                    AS bao_format,
    bao.label                       AS bao_label,
    a.assay_organism                AS assay_organism,
    a.assay_tax_id                  AS assay_tax_id,
    a.assay_strain                  AS assay_strain,
    a.assay_tissue                  AS assay_tissue,
    a.assay_cell_type               AS assay_cell_type,
    a.assay_subcellular_fraction    AS assay_subcellular_fraction,
    a.variant_id                    AS variant_id,
    vs.mutation                     AS variant_mutation,
    vs.accession                    AS variant_accession,
    vs.isoform                      AS variant_isoform,
    a.confidence_score              AS confidence_score,
    csl.description                 AS confidence_desc,
    a.relationship_type             AS relationship_type,
    rt.relationship_desc            AS relationship_desc,
    a.src_id                        AS src_id,
    s.src_short_name                AS src_short_name,
    a.src_assay_id                  AS src_assay_id,
    d.chembl_id                     AS doc_chembl_id,
    d.year                          AS doc_year,
    d.journal                       AS doc_journal,
    d.pubmed_id                     AS pubmed_id,
    d.doi                           AS doi,
    d.title                         AS doc_title,
    a.tid                           AS tid,
    td.chembl_id                    AS target_chembl_id,
    td.pref_name                    AS target_pref_name,
    td.target_type                  AS target_type
FROM assays a
JOIN target_dictionary td            ON td.tid = a.tid
LEFT JOIN assay_type at              ON at.assay_type = a.assay_type
LEFT JOIN bioassay_ontology bao      ON bao.bao_id = a.bao_format
LEFT JOIN confidence_score_lookup csl ON csl.confidence_score = a.confidence_score
LEFT JOIN relationship_type rt       ON rt.relationship_type = a.relationship_type
LEFT JOIN source s                   ON s.src_id = a.src_id
LEFT JOIN docs d                     ON d.doc_id = a.doc_id
LEFT JOIN variant_sequences vs       ON vs.variant_id = a.variant_id
WHERE a.tid IN ({placeholders})
ORDER BY a.tid, a.assay_id
"""

# 每个 assay 上实测到的活性规模。standard_type 分布是事实统计，
# 不是方向判定 —— 但它是人工判断该 assay 是否用于测 GKA 的最强线索之一。
ACT_SQL = """
SELECT assay_id,
       COUNT(*)                                              AS n_act,
       SUM(CASE WHEN pchembl_value IS NOT NULL THEN 1 ELSE 0 END) AS n_pchembl,
       COUNT(DISTINCT molregno)                              AS n_cmpd
FROM activities
WHERE assay_id IN ({placeholders})
GROUP BY assay_id
"""

STD_SQL = """
SELECT assay_id, COALESCE(standard_type, '(空)') AS st, COUNT(*) AS n
FROM activities
WHERE assay_id IN ({placeholders})
GROUP BY assay_id, st
"""

# 一个 assay 可以有多条实验参数/多个分类。为了保持「一个 assay 一行」，
# 这两张表在 Python 侧聚合成 JSON 列表写入单元格。
PARAM_SQL = """
SELECT assay_id, type, relation, value, units, text_value,
       standard_type, standard_relation, standard_value, standard_units,
       standard_text_value, comments
FROM assay_parameters
WHERE assay_id IN ({placeholders})
ORDER BY assay_id, assay_param_id
"""

# activity_comment 是自由文本备注（"Not Active"、"Dose-dependent effect" 等），
# 对判断该 assay 测的是什么有参考价值，按取值去重后聚合。
COMMENT_SQL = """
SELECT assay_id, activity_comment AS c, COUNT(*) AS n
FROM activities
WHERE assay_id IN ({placeholders}) AND activity_comment IS NOT NULL
GROUP BY assay_id, activity_comment
ORDER BY assay_id, n DESC
"""

CLASS_SQL = """
SELECT acm.assay_id, ac.assay_class_id, ac.l1, ac.l2, ac.l3, ac.class_type, ac.source
FROM assay_class_map acm
JOIN assay_classification ac ON ac.assay_class_id = acm.assay_class_id
WHERE acm.assay_id IN ({placeholders})
ORDER BY acm.assay_id, ac.assay_class_id
"""


def _compact(row: sqlite3.Row, keys: list[str]) -> dict:
    """丢掉值为空的键，避免 JSON 里塞满 null 影响可读性。"""
    return {k: row[k] for k in keys if row[k] not in (None, "")}


def _json(items: list[dict]) -> str:
    """空列表写成空字符串，便于在 CSV / pandas 里判空。"""
    return json.dumps(items, ensure_ascii=False) if items else ""


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


def load_tids(anchor_csv: Path) -> list[int]:
    """从 Step1_01 的输出里取 tid，保持两步之间的衔接是显式的。"""
    with anchor_csv.open(encoding="utf-8") as f:
        return sorted({int(r["tid"]) for r in csv.DictReader(f) if r.get("tid")})


def fetch(con: sqlite3.Connection, tids: list[int]) -> list[dict]:
    ph = ",".join("?" * len(tids))
    rows = [dict(r) for r in con.execute(ASSAY_SQL.format(placeholders=ph), tids)]
    if not rows:
        return rows

    ids = [r["assay_id"] for r in rows]
    ph2 = ",".join("?" * len(ids))

    counts = {r["assay_id"]: r for r in con.execute(ACT_SQL.format(placeholders=ph2), ids)}
    std: dict[int, Counter] = {}
    for r in con.execute(STD_SQL.format(placeholders=ph2), ids):
        std.setdefault(r["assay_id"], Counter())[r["st"]] = r["n"]

    params: dict[int, list] = {}
    pkeys = ["type", "relation", "value", "units", "text_value",
             "standard_type", "standard_relation", "standard_value",
             "standard_units", "standard_text_value", "comments"]
    for r in con.execute(PARAM_SQL.format(placeholders=ph2), ids):
        params.setdefault(r["assay_id"], []).append(_compact(r, pkeys))

    classes: dict[int, list] = {}
    ckeys = ["assay_class_id", "l1", "l2", "l3", "class_type", "source"]
    for r in con.execute(CLASS_SQL.format(placeholders=ph2), ids):
        classes.setdefault(r["assay_id"], []).append(_compact(r, ckeys))

    comments: dict = {}
    for r in con.execute(COMMENT_SQL.format(placeholders=ph2), ids):
        comments.setdefault(r["assay_id"], []).append({"comment": r["c"], "n": r["n"]})

    for row in rows:
        aid = row["assay_id"]
        row["assay_parameters"] = _json(params.get(aid, []))
        row["n_assay_parameters"] = len(params.get(aid, []))
        row["assay_classifications"] = _json(classes.get(aid, []))
        row["n_assay_classifications"] = len(classes.get(aid, []))
        row["activity_comments"] = _json(comments.get(aid, []))
        row["n_activity_comments"] = len(comments.get(aid, []))
        c = counts.get(row["assay_id"])
        row["n_activities"] = c["n_act"] if c else 0
        row["n_activities_with_pchembl"] = c["n_pchembl"] if c else 0
        row["n_distinct_compounds"] = c["n_cmpd"] if c else 0
        row["standard_types"] = "; ".join(
            f"{k}:{v}" for k, v in std.get(row["assay_id"], Counter()).most_common()
        )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in COLUMNS})


def _dist(rows: list[dict], key: str, blank: str = "(空)") -> list[tuple]:
    return Counter(r[key] if r[key] not in (None, "") else blank for r in rows).most_common()


def write_report(rows: list[dict], path: Path, tids: list[int],
                 db_path: Path, version: str) -> None:
    L: list[str] = []
    tot_act = sum(r["n_activities"] for r in rows)

    L.append("# Step1_02 GCK assay 清单")
    L.append("")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db_path}`")
    L.append(f"- 运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- 输入 tid：{', '.join(str(t) for t in tids)}")
    L.append(f"- assay 总数：**{len(rows)}**，其上活性数据点合计 **{tot_act:,}**")
    L.append("")
    L.append("> 本步骤只清点、不判定方向。下面的分布统计与完整清单用于人工判断"
             "哪些 assay 真正用于测量 GKA。")
    L.append("")

    L.append("## 按靶点")
    L.append("")
    L.append("| target | pref_name | target_type | assay 数 | 活性数 |")
    L.append("| --- | --- | --- | ---: | ---: |")
    by_t: dict[str, list] = {}
    for r in rows:
        by_t.setdefault(r["target_chembl_id"], []).append(r)
    for cid, rs in sorted(by_t.items(), key=lambda kv: -len(kv[1])):
        L.append(f"| `{cid}` | {rs[0]['target_pref_name']} | {rs[0]['target_type']} | "
                 f"{len(rs)} | {sum(x['n_activities'] for x in rs):,} |")
    L.append("")

    for title, key, note in [
        ("实验类型 assay_type", "assay_type",
         "注意：绝大多数 GCK assay 被归为 B(Binding)，**该字段在这里几乎没有区分力**，"
         "不能用它判断实验测的是激活还是抑制。"),
        ("实验格式 bao_format", "bao_label",
         "BAO 本体对实验形式的描述，比 assay_type 细一些。"),
        ("物种 assay_organism", "assay_organism", ""),
        ("靶点可信度 confidence_score", "confidence_score",
         "官方含义：9 = 直接指认到单一蛋白，8 = 同源单一蛋白，5 = 可能对应多个蛋白。"
         "它描述的是靶点归因的粒度与方式，**不是数据质量分**；"
         "每行对应的官方原文见 CSV 的 `confidence_desc` 列。"),
        ("数据来源 src_short_name", "src_short_name", ""),
    ]:
        L.append(f"## {title}")
        L.append("")
        if note:
            L.append(f"> {note}")
            L.append("")
        L.append("| 取值 | assay 数 | 活性数 |")
        L.append("| --- | ---: | ---: |")
        for val, n in _dist(rows, key):
            act = sum(r["n_activities"] for r in rows
                      if (r[key] if r[key] not in (None, "") else "(空)") == val)
            L.append(f"| {val} | {n} | {act:,} |")
        L.append("")

    # standard_type 全局分布 —— 判断 GKA 实验最有用的线索
    st_all: Counter = Counter()
    for r in rows:
        for part in filter(None, r["standard_types"].split("; ")):
            k, v = part.rsplit(":", 1)
            st_all[k] += int(v)
    L.append("## 全部 assay 上实测的 standard_type")
    L.append("")
    L.append("> 这是**事实统计**，不是方向判定。但它是判断某个 assay 是否在测 GKA 的最强线索："
             "`EC50`/`Emax`/`%max`/`S0.5`/`Vmax` 指向激活与酶动力学表征，"
             "`IC50`/`Inhibition` 指向抑制，`Kd` 指向结合。")
    L.append("")
    L.append("| standard_type | 活性数 |")
    L.append("| --- | ---: |")
    for k, v in st_all.most_common():
        L.append(f"| `{k}` | {v:,} |")
    L.append("")

    # 三个聚合/可选字段的覆盖情况 —— 空值本身也是需要如实记录的事实
    n_var = sum(1 for r in rows if r["variant_id"])
    n_par = sum(1 for r in rows if r["n_assay_parameters"])
    n_cls = sum(1 for r in rows if r["n_assay_classifications"])
    L.append("## 突变体 / 实验参数 / 实验分类的覆盖情况")
    L.append("")
    L.append("| 字段 | 有值的 assay 数 | 说明 |")
    L.append("| --- | ---: | --- |")
    L.append(f"| `variant_id` | {n_var} / {len(rows)} | "
             "非空表示该实验用的是突变体蛋白（如激酶耐药突变）；"
             "GCK 的天然激活突变若被单独建为 variant 会出现在这里 |")
    L.append(f"| `assay_parameters` | {n_par} / {len(rows)} | "
             "实验参数，已聚合为 JSON 列表 |")
    L.append(f"| `assay_classifications` | {n_cls} / {len(rows)} | "
             "实验分类，已聚合为 JSON 列表；该表主要覆盖体内/治疗领域实验 |")
    L.append("")
    if n_par:
        ptypes = Counter()
        for r in rows:
            for item in json.loads(r["assay_parameters"] or "[]"):
                ptypes[item.get("type")] += 1
        L.append("实际出现的参数类型：" +
                 "、".join(f"`{k}`（{v}）" for k, v in ptypes.most_common()))
        L.append("")
        L.append("> 注意：这里的参数**不是**葡萄糖浓度、孵育时间这类实验条件，"
                 "而是 EFO/MONDO 疾病本体标注。GCK 这批 assay 的实验条件"
                 "（葡萄糖浓度等）只写在 `assay_description` 自由文本里，"
                 "没有被结构化到 `assay_parameters`。要按葡萄糖浓度分层，"
                 "必须自己解析描述文本。")
        L.append("")

    L.append("## 文献年份分布")
    L.append("")
    L.append("| 年份 | assay 数 |")
    L.append("| ---: | ---: |")
    for y, n in sorted(_dist(rows, "doc_year"), key=lambda kv: str(kv[0])):
        L.append(f"| {y} | {n} |")
    L.append("")

    L.append("## 完整 assay 清单")
    L.append("")
    L.append("按活性数据量降序。完整字段见同目录 CSV。")
    L.append("")
    L.append("| assay_chembl_id | type | conf | 活性数 | 化合物数 | standard_types | 描述 |")
    L.append("| --- | --- | ---: | ---: | ---: | --- | --- |")
    for r in sorted(rows, key=lambda x: -x["n_activities"]):
        desc = (r["assay_description"] or "").replace("|", "\\|")
        L.append(
            f"| `{r['assay_chembl_id']}` | {r['assay_type']} | {r['confidence_score']} | "
            f"{r['n_activities']:,} | {r['n_distinct_compounds']:,} | "
            f"{r['standard_types'] or '—'} | {desc} |"
        )
    L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(description="提取 GCK 靶点关联的全部 assay 并汇总。")
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--anchor-csv", type=Path, default=DEFAULT_ANCHOR,
                   help=f"Step1_01 的输出，默认 {DEFAULT_ANCHOR}")
    p.add_argument("--tid", type=int, nargs="+", default=None,
                   help="直接指定 tid，给出后忽略 --anchor-csv")
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    if not args.db.is_file():
        print(f"错误：找不到数据库 {args.db}", file=sys.stderr)
        return 1

    if args.tid:
        tids = sorted(set(args.tid))
        origin = "命令行 --tid"
    else:
        if not args.anchor_csv.is_file():
            print(f"错误：找不到 Step1_01 的输出 {args.anchor_csv}，"
                  "请先运行 Step1_01 或用 --tid 指定。", file=sys.stderr)
            return 1
        tids = load_tids(args.anchor_csv)
        origin = str(args.anchor_csv)

    con = connect_readonly(args.db)
    version = chembl_version(con)
    print(f"数据库：{args.db}")
    print(f"ChEMBL 版本：{version}")
    print(f"tid 来源：{origin}")
    print(f"tid：{tids}\n")

    rows = fetch(con, tids)
    if not rows:
        print("错误：这些 tid 下没有找到任何 assay。", file=sys.stderr)
        return 2

    for cid, n in Counter(r["target_chembl_id"] for r in rows).most_common():
        act = sum(r["n_activities"] for r in rows if r["target_chembl_id"] == cid)
        print(f"  {cid}: {n} 个 assay，{act:,} 条活性")
    print(f"\nassay 合计 {len(rows)} 个，活性合计 "
          f"{sum(r['n_activities'] for r in rows):,} 条")

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "Step1_02_GCK_Assay_Mapping.csv"
    md_path = args.outdir / "Step1_02_GCK_Assay_Mapping.md"
    write_csv(rows, csv_path)
    write_report(rows, md_path, tids, args.db, version)
    print(f"\n清单：{csv_path}")
    print(f"报告：{md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
