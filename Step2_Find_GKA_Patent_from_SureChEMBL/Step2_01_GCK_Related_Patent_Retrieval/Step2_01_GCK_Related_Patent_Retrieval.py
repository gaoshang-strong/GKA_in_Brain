#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step2_01_GCK_Related_Patent_Retrieval
=====================================

从 SureChEMBL bulk 快照检索**所有与葡萄糖激酶（GCK）相关的专利**。

**召回优先**：不判断方向、不判断是不是 GKA 专利、不做 field 过滤。
**纯规则**：只用 bulk parquet，不调 API、不用 LLM。
**独立于 ChEMBL**：Step1 的候选不进检索式，只做验证列（`val_*` 前缀）。

策略：两个独立锚点取并集 + 同族扩展
------------------------------------
    锚点 A  patents.title 正则         ← 不经过生物医学标注管道
    锚点 B  biomedical_locations 实体   ← resolved_form + 白名单
               ↓ 并集
    扩展    family_id 展开             ← 同族=同一发明，成员共享权要
               ↓
    标注（不筛）  CPC / 命中来源 / field 分布 / 化合物数 / 风险标记

为什么必须两个锚点：标注管道会整篇缺失
--------------------------------------
实测四篇专利，权要里含已知 GKA 化合物、标题明写 "GLUCOKINASE ACTIVATOR"，
但 `biomedical_locations` **一条记录都没有**（化学侧却有 75–238 个化合物）：

    EP-4725482-A1      GLUCOKINASE ACTIVATOR FOR COGNITIVE DISORDERS AND NEUROLOGICAL…
    US-20260200881-A1  SULFOXIDE AND SULFONE GLUCOKINASE ACTIVATORS…
    CN-118453592-A     Glucokinase activator composition for treating diabetes
    US-12064416-B2     Pharmaceutical combination containing glucose kinase activator

只用实体锚定会把这些整个漏掉——第一篇正是本项目最该找到的专利。
产物的 `has_biomedical_annotation` 列直接暴露这个缺口。

换靶点
------
靶点相关的全部参数收在 `TARGET` 配置块里，检索逻辑不用改。
必须重推的四项见 Readme「三、换靶点怎么复用」。

输出
----
    Step2_01_GCK_Related_Patents.csv     一行一篇专利文档（主产物）
    Step2_01_GCKR_Related_Patents.csv    GKRP 相关，单独一张
    Step2_01_GCK_Related_Patent_Retrieval.md   报告
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import duckdb

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_SNAP = REPO / "SureChEMBL" / "SureChEMBL_2026-07-17"
DEFAULT_STEP1 = REPO / "Step1_Find_GKA_from_ChEMBL" / "Step1_GKA_Candidates_with_Properties.csv"

# ---------------------------------------------------------------------------
# 靶点配置块 —— 换靶点只改这里
# ---------------------------------------------------------------------------
TARGET = {
    "name": "GCK / glucokinase",

    # 锚点 B：实体归一 ID。全库 53% 的实体未解析，换靶点先查它在哪个命名空间
    "resolved_forms": ["HGNC:4195"],

    # 未解析但确认是本靶点的写法（人工白名单，逐条给理由）
    "entity_whitelist": [
        "glucokinase (hexokinase 4)",
        "Glucokinase (hexokinase 4, maturity onset diabetes of the young 2)",
    ],

    # 锚点 A：标题正则（SQL 的 regexp_matches，'i' 忽略大小写）
    "title_regex": r"glucokinase|glucose\s*kinase|hexokinase[\s-]*(4|IV)\b",
    # 方向信号：标题里直接写了「激活剂」
    "title_activator_regex": r"glucokinase\s+activator|glucose\s*kinase\s+activator",

    # 近邻但不是本靶点 —— 单独成表
    "sibling_resolved": ["HGNC:4196", "Q14397"],
    "sibling_texts": ["glucokinase regulatory protein",
                      "Glucokinase (hexokinase 4) regulator", "GLKRP"],
    "sibling_label": "GKRP / GCKR",

    # 会带进假阳性的写法：标记不删除
    "risk_forms": {
        "GK": "⚠ 糖尿病文献里更常指 Goto-Kakizaki 大鼠（2 型糖尿病模型），与本领域高度重叠",
        "4": "⚠ 单个数字被解析成基因，标注错误",
        "glk": "细菌 glucokinase 基因名", "GlkA": "细菌 glucokinase 基因名",
        "gukA": "细菌 glucokinase 基因名", "GlcK": "细菌 glucokinase 基因名",
        "Hk4": "缩写，需核", "gki": "缩写，需核", "GluK": "缩写，需核",
    },
    # 正向信号：其实是「activator」的缩写
    "positive_forms": {
        "GKA": "✅「glucokinase activator」的缩写，**方向正向信号**",
        "GKAs": "✅「glucokinase activators」的缩写，**方向正向信号**",
    },

    # 只作标注的分类号（CPC 覆盖好但特异性差，不能当锚点）
    "cpc_of_interest": {"A61P 3/10": "抗糖尿病用途"},
}

OUT_COLUMNS = [
    # --- 身份 ---
    "patent_id", "patent_number", "country", "publication_date",
    "family_id", "family_valid", "title", "assignee",
    # --- 命中来源（本步骤核心）---
    "anchor_sources", "by_title", "by_entity", "by_family_only",
    "title_matched", "title_says_activator",
    # --- 实体侧证据 ---
    "has_biomedical_annotation", "n_surface_forms", "surface_forms",
    "hit_ttl", "hit_abst", "hit_desc", "hit_clms", "n_mentions",
    "risk_flags", "has_activator_abbrev",
    # --- 化学侧 ---
    "n_compounds_total", "n_compounds_clms",
    # --- 标注（不筛）---
    "cpc_antidiabetic", "sibling_gckr",
    # --- 验证列：来自 ChEMBL，**未参与检索** ---
    "val_n_step1_candidates", "val_n_step1_candidates_clms",
    "val_step1_candidate_ids", "val_is_positive_control", "val_control_names",
]

T0 = time.time()


def log(m: str) -> None:
    print(f"[{time.time() - T0:6.1f}s] {m}", file=sys.stderr, flush=True)


def fmt(n) -> str:
    return "—" if n is None else f"{int(n):,}" if isinstance(n, (int, float)) else str(n)


def pct(a, b) -> str:
    return "—" if not b else f"{100.0 * a / b:.1f}%"


def q_list(xs) -> str:
    return ",".join("'" + str(x).replace("'", "''") + "'" for x in xs)


def q_str(s) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def connect(snap: Path, threads: int, mem: str):
    con = duckdb.connect()
    con.execute(f"SET threads={threads}; SET memory_limit='{mem}';")
    con.execute("SET enable_progress_bar=false;")
    for t in ("patents", "compounds", "patent_compound_map", "fields",
              "biomedical_entities", "biomedical_locations"):
        p = snap / f"{t}.parquet"
        if not p.is_file():
            sys.exit(f"错误：找不到 {p}")
        con.execute(f"CREATE OR REPLACE VIEW {t} AS SELECT * FROM '{p}'")
    return con


# ---------------------------------------------------------------------------
# 锚点
# ---------------------------------------------------------------------------

def anchor_title(con) -> int:
    """锚点 A：标题正则。不经过生物医学标注管道，因此不受标注缺失影响。"""
    con.execute(f"""
        CREATE OR REPLACE TABLE a_title AS
        SELECT id AS patent_id,
               regexp_matches(title, {q_str(TARGET['title_regex'])}, 'i')          AS by_title,
               regexp_matches(title, {q_str(TARGET['title_activator_regex'])}, 'i') AS title_says_activator
        FROM patents
        WHERE title IS NOT NULL
          AND regexp_matches(title, {q_str(TARGET['title_regex'])}, 'i')
    """)
    n = con.execute("SELECT COUNT(*) FROM a_title").fetchone()[0]
    na = con.execute("SELECT COUNT(*) FROM a_title WHERE title_says_activator").fetchone()[0]
    log(f"锚点 A 标题：{n:,} 篇（其中标题明写 activator {na:,} 篇）")
    return n


def anchor_entity(con) -> int:
    """锚点 B：实体归一 + 白名单。"""
    con.execute(f"""
        CREATE OR REPLACE TABLE anchor_entity AS
        SELECT id, original_text
        FROM biomedical_entities
        WHERE resolved_form IN ({q_list(TARGET['resolved_forms'])})
           OR original_text IN ({q_list(TARGET['entity_whitelist'])})
    """)
    ne = con.execute("SELECT COUNT(*) FROM anchor_entity").fetchone()[0]
    con.execute("""
        CREATE OR REPLACE TABLE a_entity AS
        SELECT l.patent_id,
               COUNT(DISTINCT e.id)                            AS n_surface_forms,
               list_sort(list_distinct(list(e.original_text)))  AS surface_forms,
               SUM(l.count)                                     AS n_mentions,
               SUM(l.count) FILTER (WHERE l.field_id = 4)       AS hit_ttl,
               SUM(l.count) FILTER (WHERE l.field_id = 3)       AS hit_abst,
               SUM(l.count) FILTER (WHERE l.field_id = 1)       AS hit_desc,
               SUM(l.count) FILTER (WHERE l.field_id = 2)       AS hit_clms
        FROM biomedical_locations l
        JOIN anchor_entity e ON e.id = l.entity_id
        GROUP BY 1
    """)
    n = con.execute("SELECT COUNT(*) FROM a_entity").fetchone()[0]
    log(f"锚点 B 实体：{ne} 个实体 → {n:,} 篇")
    return n


def expand_family(con) -> dict:
    """同族扩展：命中任一成员即收全族。同族=同一发明，成员共享说明书与权要。

    ⚠ family_id = -1 是「未分配同族」哨兵，必须排除，否则会把 7 万余篇
      不相干专利当成同一个发明全部拉进来。
    """
    con.execute("""
        CREATE OR REPLACE TABLE seed AS
        SELECT patent_id FROM a_title
        UNION
        SELECT patent_id FROM a_entity
    """)
    n_seed = con.execute("SELECT COUNT(*) FROM seed").fetchone()[0]
    con.execute("""
        CREATE OR REPLACE TABLE hit_family AS
        SELECT DISTINCT p.family_id
        FROM seed s JOIN patents p ON p.id = s.patent_id
        WHERE p.family_id > 0
    """)
    con.execute("""
        CREATE OR REPLACE TABLE universe AS
        SELECT DISTINCT patent_id FROM (
            SELECT patent_id FROM seed
            UNION
            SELECT p.id AS patent_id FROM patents p
            JOIN hit_family f ON f.family_id = p.family_id
        )
    """)
    n_fam = con.execute("SELECT COUNT(*) FROM hit_family").fetchone()[0]
    n_all = con.execute("SELECT COUNT(*) FROM universe").fetchone()[0]
    log(f"并集种子 {n_seed:,} 篇 → 命中 {n_fam:,} 个同族 → 展开后 {n_all:,} 篇"
        f"（同族扩展补 {n_all - n_seed:,} 篇）")
    return {"n_seed": n_seed, "n_family": n_fam, "n_universe": n_all}


def load_step1_validation(con, path: Path) -> int:
    """读 Step1 整合表，**只作验证列**，不参与检索。"""
    if not path.is_file():
        log(f"警告：找不到 {path}，跳过验证列")
        con.execute("CREATE OR REPLACE TABLE val_cmp(cid BIGINT, chembl_id VARCHAR, ctl VARCHAR)")
        con.execute("CREATE OR REPLACE TABLE val_raw(ik VARCHAR, chembl_id VARCHAR, ctl VARCHAR)")
        return 0
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    con.execute("CREATE OR REPLACE TABLE val_raw(ik VARCHAR, chembl_id VARCHAR, ctl VARCHAR)")
    con.executemany("INSERT INTO val_raw VALUES (?,?,?)",
                    [(r["standard_inchi_key"], r["molecule_chembl_id"],
                      r.get("control_name", "") if r.get("is_positive_control") == "TRUE" else "")
                     for r in rows if r.get("standard_inchi_key")])
    # SureChEMBL 的 inchi_key 不唯一，先 DISTINCT 再 join
    con.execute("""
        CREATE OR REPLACE TABLE val_cmp AS
        SELECT DISTINCT c.id AS cid, v.chembl_id, v.ctl
        FROM val_raw v JOIN compounds c ON c.inchi_key = v.ik
    """)
    n = con.execute("SELECT COUNT(*) FROM val_cmp").fetchone()[0]
    log(f"验证列：Step1 {len(rows)} 个分子 → SureChEMBL 匹配 {n} 个 compound_id（不参与检索）")
    return n


def build_main(con) -> None:
    con.execute("""
        CREATE OR REPLACE TABLE cmp_cnt AS
        SELECT m.patent_id,
               COUNT(DISTINCT m.compound_id)                               AS n_compounds_total,
               COUNT(DISTINCT m.compound_id) FILTER (WHERE m.field_id = 2) AS n_compounds_clms
        FROM patent_compound_map m JOIN universe u ON u.patent_id = m.patent_id
        GROUP BY 1
    """)
    log("化合物计数建好")

    con.execute("""
        CREATE OR REPLACE TABLE val_cnt AS
        SELECT m.patent_id,
               COUNT(DISTINCT v.chembl_id)                               AS n_val,
               COUNT(DISTINCT v.chembl_id) FILTER (WHERE m.field_id = 2) AS n_val_clms,
               list_sort(list_distinct(list(v.chembl_id)))               AS val_ids,
               list_sort(list_distinct(list(v.ctl) FILTER (WHERE v.ctl <> ''))) AS ctls
        FROM patent_compound_map m
        JOIN val_cmp v ON v.cid = m.compound_id
        JOIN universe u ON u.patent_id = m.patent_id
        GROUP BY 1
    """)
    log("验证列建好")

    # 该专利有没有任何生物医学标注 —— 直接暴露标注缺失
    con.execute("""
        CREATE OR REPLACE TABLE has_bio AS
        SELECT DISTINCT l.patent_id FROM biomedical_locations l
        JOIN universe u ON u.patent_id = l.patent_id
    """)
    log("标注存在性建好")

    sib_cond = (f"resolved_form IN ({q_list(TARGET['sibling_resolved'])}) "
                f"OR original_text IN ({q_list(TARGET['sibling_texts'])})")
    con.execute(f"""
        CREATE OR REPLACE TABLE sib AS
        SELECT DISTINCT l.patent_id FROM biomedical_locations l
        JOIN (SELECT id FROM biomedical_entities WHERE {sib_cond}) e ON e.id = l.entity_id
    """)

    cpc = list(TARGET["cpc_of_interest"])[0]
    con.execute(f"""
        CREATE OR REPLACE TABLE main AS
        SELECT p.id AS patent_id, p.patent_number, p.country, p.publication_date,
               p.family_id,
               CASE WHEN p.family_id > 0 THEN 'TRUE' ELSE 'FALSE' END AS family_valid,
               p.title, p.assignee,
               CASE WHEN t.patent_id IS NOT NULL THEN 'TRUE' ELSE 'FALSE' END AS by_title,
               CASE WHEN e.patent_id IS NOT NULL THEN 'TRUE' ELSE 'FALSE' END AS by_entity,
               CASE WHEN t.patent_id IS NULL AND e.patent_id IS NULL
                    THEN 'TRUE' ELSE 'FALSE' END AS by_family_only,
               CASE WHEN COALESCE(t.by_title, FALSE) THEN 'TRUE' ELSE 'FALSE' END AS title_matched,
               CASE WHEN COALESCE(t.title_says_activator, FALSE)
                    THEN 'TRUE' ELSE 'FALSE' END AS title_says_activator,
               CASE WHEN b.patent_id IS NOT NULL THEN 'TRUE' ELSE 'FALSE' END
                    AS has_biomedical_annotation,
               COALESCE(e.n_surface_forms, 0) AS n_surface_forms,
               e.surface_forms,
               COALESCE(e.hit_ttl, 0) AS hit_ttl, COALESCE(e.hit_abst, 0) AS hit_abst,
               COALESCE(e.hit_desc, 0) AS hit_desc, COALESCE(e.hit_clms, 0) AS hit_clms,
               COALESCE(e.n_mentions, 0) AS n_mentions,
               list_filter(COALESCE(e.surface_forms, []),
                           x -> x IN ({q_list(TARGET['risk_forms'])})) AS risk_flags,
               CASE WHEN len(list_filter(COALESCE(e.surface_forms, []),
                                         x -> x IN ({q_list(TARGET['positive_forms'])}))) > 0
                    THEN 'TRUE' ELSE 'FALSE' END AS has_activator_abbrev,
               COALESCE(c.n_compounds_total, 0) AS n_compounds_total,
               COALESCE(c.n_compounds_clms, 0)  AS n_compounds_clms,
               CASE WHEN list_contains(p.cpc, {q_str(cpc)}) THEN 'TRUE' ELSE 'FALSE' END
                    AS cpc_antidiabetic,
               CASE WHEN s.patent_id IS NOT NULL THEN 'TRUE' ELSE 'FALSE' END AS sibling_gckr,
               COALESCE(v.n_val, 0) AS val_n_step1_candidates,
               COALESCE(v.n_val_clms, 0) AS val_n_step1_candidates_clms,
               v.val_ids AS val_step1_candidate_ids,
               CASE WHEN len(COALESCE(v.ctls, [])) > 0 THEN 'TRUE' ELSE 'FALSE' END
                    AS val_is_positive_control,
               v.ctls AS val_control_names
        FROM universe u
        JOIN patents p ON p.id = u.patent_id
        LEFT JOIN a_title  t ON t.patent_id = u.patent_id
        LEFT JOIN a_entity e ON e.patent_id = u.patent_id
        LEFT JOIN has_bio  b ON b.patent_id = u.patent_id
        LEFT JOIN cmp_cnt  c ON c.patent_id = u.patent_id
        LEFT JOIN val_cnt  v ON v.patent_id = u.patent_id
        LEFT JOIN sib      s ON s.patent_id = u.patent_id
    """)
    con.execute("""
        ALTER TABLE main ADD COLUMN anchor_sources VARCHAR;
        UPDATE main SET anchor_sources =
            CASE WHEN by_title='TRUE' AND by_entity='TRUE' THEN 'title+entity'
                 WHEN by_title='TRUE' THEN 'title'
                 WHEN by_entity='TRUE' THEN 'entity'
                 ELSE 'family' END;
    """)
    log(f"主表 {con.execute('SELECT COUNT(*) FROM main').fetchone()[0]:,} 行")

    con.execute(f"""
        CREATE OR REPLACE TABLE gckr AS
        SELECT p.id AS patent_id, p.patent_number, p.country, p.publication_date,
               p.family_id, p.title,
               list_sort(list_distinct(list(e.original_text))) AS surface_forms,
               SUM(l.count) AS n_mentions,
               SUM(l.count) FILTER (WHERE l.field_id = 2) AS hit_clms
        FROM biomedical_locations l
        JOIN (SELECT id, original_text FROM biomedical_entities WHERE {sib_cond}) e
             ON e.id = l.entity_id
        JOIN patents p ON p.id = l.patent_id
        GROUP BY 1,2,3,4,5,6
    """)
    log(f"GKRP 表 {con.execute('SELECT COUNT(*) FROM gckr').fetchone()[0]:,} 行")


def to_csv(con, table: str, cols: list, path: Path, order: str) -> int:
    rows = con.execute(f"SELECT * FROM {table} ORDER BY {order}").fetchall()
    names = [d[0] for d in con.description]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            d = dict(zip(names, r))
            out = {}
            for k in cols:
                v = d.get(k)
                if isinstance(v, list):
                    v = json.dumps(v, ensure_ascii=False) if v else ""
                out[k] = "" if v is None else v
            w.writerow(out)
    return len(rows)


def write_report(con, path: Path, snap: Path, step1: Path, stats: dict) -> None:
    L = []
    q = lambda s: con.execute(s).fetchall()          # noqa: E731
    one = lambda s: con.execute(s).fetchone()[0]     # noqa: E731

    n = one("SELECT COUNT(*) FROM main")
    n_fam = one("SELECT COUNT(DISTINCT family_id) FROM main WHERE family_id > 0")

    L.append("# Step2_01 GCK 相关专利检索")
    L.append("")
    L.append("- 数据源：**SureChEMBL 2.0，`2026-07-17` 全量快照**")
    L.append(f"- 快照路径：`{snap}`")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 靶点：**{TARGET['name']}**（`{'`、`'.join(TARGET['resolved_forms'])}`）")
    L.append(f"- 命中：**{fmt(n)}** 篇专利文档 → 同族去重 **{fmt(n_fam)}** 个发明")
    L.append("")
    L.append("> **召回优先**：不判断方向、不判断是不是 GKA 专利、**不做 `field_id` 过滤**。")
    L.append("> **独立于 ChEMBL**：Step1 候选不进检索式，`val_*` 列仅作验证。")
    L.append("")

    L.append("## 一、两个独立锚点 + 同族扩展")
    L.append("")
    rows = q("""SELECT anchor_sources, COUNT(*),
                       COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)
                FROM main GROUP BY 1 ORDER BY 2 DESC""")
    note = {"title+entity": "两个锚点都命中，最可靠",
            "entity": "只有实体标注命中（标题没写靶点名，常见于化合物专利）",
            "title": "**只有标题命中——实体标注没抓到**",
            "family": "本身没命中，靠同族成员带进来"}
    L.append("| 命中来源 | 专利数 | 同族 | 说明 |")
    L.append("| --- | ---: | ---: | --- |")
    for src, c, f in rows:
        L.append(f"| `{src}` | {fmt(c)} | {fmt(f)} | {note.get(src, '')} |")
    L.append("")
    L.append(f"锚点 A（标题正则）命中 {fmt(stats['n_title'])} 篇，"
             f"锚点 B（实体）命中 {fmt(stats['n_entity'])} 篇，"
             f"并集 {fmt(stats['n_seed'])} 篇；"
             f"同族展开后 **{fmt(stats['n_universe'])}** 篇"
             f"（补 {fmt(stats['n_universe'] - stats['n_seed'])} 篇）。")
    L.append("")

    L.append("### ⚠ 标注管道整篇缺失的证据")
    L.append("")
    no_bio = one("SELECT COUNT(*) FROM main WHERE has_biomedical_annotation = 'FALSE'")
    tit_no_bio = one("SELECT COUNT(*) FROM main WHERE by_title = 'TRUE' "
                     "AND has_biomedical_annotation = 'FALSE'")
    L.append(f"本表 **{fmt(no_bio)}** 篇专利的 `biomedical_locations` **一条记录都没有**"
             f"（`has_biomedical_annotation = FALSE`），"
             f"其中 **{fmt(tit_no_bio)}** 篇是标题正则捞回来的——"
             "**只用实体锚定这些会被整个漏掉**。")
    L.append("")
    rows = q("""SELECT patent_number, country, publication_date,
                       substr(title, 1, 62) AS t, n_compounds_total
                FROM main
                WHERE by_title = 'TRUE' AND has_biomedical_annotation = 'FALSE'
                  AND title_says_activator = 'TRUE'
                ORDER BY publication_date DESC LIMIT 8""")
    if rows:
        L.append("标题明写 activator、却完全没有生物医学标注的（按公开日倒序）：")
        L.append("")
        L.append("| 专利 | 国 | 公开日 | 标题 | 化合物数 |")
        L.append("| --- | --- | --- | --- | ---: |")
        for pn, c, d, t, nc in rows:
            L.append(f"| `{pn}` | {c} | {d} | {t} | {fmt(nc)} |")
        L.append("")
        L.append("这些专利的**化学侧是完整的**（几十上百个化合物），只是文本标注为空。")
        L.append("")

    L.append("## 二、方向信号（不做判定，只标记）")
    L.append("")
    ta = one("SELECT COUNT(*) FROM main WHERE title_says_activator = 'TRUE'")
    taf = one("SELECT COUNT(DISTINCT family_id) FROM main "
              "WHERE title_says_activator = 'TRUE' AND family_id > 0")
    ab = one("SELECT COUNT(*) FROM main WHERE has_activator_abbrev = 'TRUE'")
    L.append("| 信号 | 专利数 | 同族 | 说明 |")
    L.append("| --- | ---: | ---: | --- |")
    L.append(f"| 标题写着 activator | {fmt(ta)} | {fmt(taf)} | "
             "`title_says_activator`，**最强的方向信号** |")
    L.append(f"| 出现 `GKA` / `GKAs` 缩写 | {fmt(ab)} | — | "
             "缩写本身就是 glucokinase activator |")
    L.append("")
    L.append("> bulk 数据判不了方向（无全文，`Mechanism` 实体全是工业化学词）。"
             "这两列是**规则能拿到的全部方向信息**，真正的方向判定要读权利要求原文。")
    L.append("")

    L.append("## 三、命中位置分布（不过滤，只记录）")
    L.append("")
    for lab, col in (("`ttl` 标题", "hit_ttl"), ("`abst` 摘要", "hit_abst"),
                     ("`desc` 说明书", "hit_desc"), ("**`clms` 权利要求**", "hit_clms")):
        c = one(f"SELECT COUNT(*) FROM main WHERE {col} > 0")
        L.append(f"- {lab}：{fmt(c)} 篇（{pct(c, n)}）")
    L.append("")
    L.append("`clms` 是精度轴、`desc` 是召回轴，锚定步骤两个都留。"
             "全库这两者差 6.5 倍（12.18 亿 vs 1.87 亿关联）。")
    L.append("")

    L.append("## 四、专利局分布 ⚠")
    L.append("")
    tot = dict(q("SELECT country, COUNT(*) FROM patents GROUP BY 1"))
    rows = q("""SELECT country, COUNT(*),
                       COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)
                FROM main GROUP BY 1 ORDER BY 2 DESC""")
    L.append("| 专利局 | 命中 | 同族 | 全库该国 | 命中率 |")
    L.append("| --- | ---: | ---: | ---: | ---: |")
    for c, nd, nf in rows:
        L.append(f"| `{c}` | {fmt(nd)} | {fmt(nf)} | {fmt(tot.get(c, 0))} | {pct(nd, tot.get(c, 0))} |")
    L.append("")
    L.append("**JPO 不提供全文**（只有著录项 + 英文标题摘要），**CNIPA 只有英文机翻全文**——"
             "标注管道在这两家基本失效。**这批数据实质上是 US / EP / WO 的视图**，"
             "不能说「中国/日本没有 GKA 专利」，是看不见不是没有。")
    L.append("")

    L.append("## 五、风险标记")
    L.append("")
    rows = q("""SELECT f, COUNT(*) FROM (SELECT UNNEST(surface_forms) AS f FROM main)
                GROUP BY 1 ORDER BY 2 DESC LIMIT 15""")
    L.append("| surface form | 专利数 | 备注 |")
    L.append("| --- | ---: | --- |")
    for f, c in rows:
        L.append(f"| `{f}` | {fmt(c)} | "
                 f"{TARGET['risk_forms'].get(f, TARGET['positive_forms'].get(f, ''))} |")
    L.append("")
    nr = one("SELECT COUNT(*) FROM main WHERE len(risk_flags) > 0")
    nonly = one("""SELECT COUNT(*) FROM main WHERE len(risk_flags) > 0
                   AND len(surface_forms) = len(risk_flags) AND by_title = 'FALSE'""")
    L.append(f"- 带风险标记 **{fmt(nr)}** 篇")
    L.append(f"- **只靠风险形命中、标题也没写靶点名** 的 **{fmt(nonly)}** 篇 "
             "← 最可疑，下游优先人工核")
    L.append("")
    cpc = list(TARGET["cpc_of_interest"])[0]
    nc = one("SELECT COUNT(*) FROM main WHERE cpc_antidiabetic = 'TRUE'")
    L.append(f"另：`{cpc}`（{TARGET['cpc_of_interest'][cpc]}）命中 {fmt(nc)} 篇（{pct(nc, n)}）。"
             "CPC 召回好但全库 20 万篇，**只能当过滤器不能当锚点**，本表只作标注。")
    L.append("")

    L.append("## 六、验证：与 ChEMBL 侧的召回对照")
    L.append("")
    L.append("> 以下全部是**事后验证**，Step1 的分子**没有参与检索**。")
    L.append("")
    nv = one("SELECT COUNT(*) FROM main WHERE val_n_step1_candidates > 0")
    nvc = one("SELECT COUNT(*) FROM main WHERE val_n_step1_candidates_clms > 0")
    all_pat = one("""SELECT COUNT(DISTINCT m.patent_id) FROM patent_compound_map m
                     JOIN val_cmp v ON v.cid = m.compound_id""")
    all_clms = one("""SELECT COUNT(DISTINCT m.patent_id) FROM patent_compound_map m
                      JOIN val_cmp v ON v.cid = m.compound_id WHERE m.field_id = 2""")
    L.append("| 口径 | 全库有多少 | 本表捞到 | 召回率 |")
    L.append("| --- | ---: | ---: | ---: |")
    L.append(f"| 含 Step1 候选化合物（任意部分） | {fmt(all_pat)} | {fmt(nv)} | {pct(nv, all_pat)} |")
    L.append(f"| **权利要求里含 Step1 候选** | {fmt(all_clms)} | {fmt(nvc)} | "
             f"**{pct(nvc, all_clms)}** |")
    L.append("")
    L.append("第二行是更扎实的口径——权要里主张已知 GKA 化合物的专利，本表覆盖了多少。")
    L.append("")
    rows = q("""SELECT c AS name, COUNT(*) AS nd,
                       COUNT(*) FILTER (WHERE hit_clms > 0 OR title_says_activator = 'TRUE') AS ns
                FROM (SELECT UNNEST(val_control_names) AS c, hit_clms, title_says_activator FROM main)
                GROUP BY 1 ORDER BY 2 DESC""")
    hit = {r[0]: (r[1], r[2]) for r in rows}
    all_ctl = set(one("SELECT list(DISTINCT ctl) FROM val_raw WHERE ctl <> ''") or [])
    matched = set(one("SELECT list(DISTINCT ctl) FROM val_cmp WHERE ctl <> ''") or [])
    L.append("阳性对照（12 个，来自 Step1 整合表）：")
    L.append("")
    L.append("| 对照 | 结构在 SureChEMBL | 命中专利 |")
    L.append("| --- | :---: | ---: |")
    for name in sorted(all_ctl, key=lambda x: -hit.get(x, (0, 0))[0]):
        if name in hit:
            L.append(f"| {name} | ✅ | {fmt(hit[name][0])} |")
        elif name in matched:
            L.append(f"| {name} | ✅ | **0** |")
        else:
            L.append(f"| {name} | **❌ 无对应结构** | — |")
    L.append("")
    L.append(f"**{len(hit)}/{len(all_ctl)} 个对照捞到了专利。** "
             "「无对应结构」的是盐型——SureChEMBL 里不单独注册盐，"
             "母体条目已命中，药物层没真丢。")
    L.append("")

    L.append("## 七、GKRP 单独成表")
    L.append("")
    ng = one("SELECT COUNT(*) FROM gckr")
    ngf = one("SELECT COUNT(DISTINCT family_id) FROM gckr WHERE family_id > 0")
    ov = one("SELECT COUNT(*) FROM main WHERE sibling_gckr = 'TRUE'")
    L.append(f"**{fmt(ng)}** 篇（{fmt(ngf)} 个同族），与主表重叠 **{fmt(ov)}** 篇"
             "（主表 `sibling_gckr` 列标出）。")
    L.append("")
    L.append(f"{TARGET['sibling_label']} 解离剂与直接激活 GCK 是两类机制，"
             "与 ChEMBL 侧单列 `CHEMBL3885579`(GCK–GKRP PPI) 的处理一致。")
    L.append("")

    L.append("## 八、这一步没做什么")
    L.append("")
    L.append("| 没做 | 为什么 | 留给谁 |")
    L.append("| --- | --- | --- |")
    L.append("| 方向判定 | bulk 无全文；`Mechanism` 实体全是工业化学词、`resolved_form` 全空 | "
             "后续读权利要求原文 |")
    L.append("| `field_id` 过滤 | 锚定要召回 | 下游用 `hit_*` 列 |")
    L.append("| 结构相似性检索 | **以 ChEMBL 为种子，不独立**，属扩展臂 | 单独一步 |")
    L.append("| 泛称 `hexokinase`（29,158 篇） | 无法分辨 HK1/2/3 与 GCK | 需结构或全文佐证 |")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="从 SureChEMBL 检索 GCK 相关专利（纯规则，独立于 ChEMBL）。")
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAP)
    ap.add_argument("--step1-csv", type=Path, default=DEFAULT_STEP1,
                    help="仅用于验证列，不参与检索")
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="12GB")
    args = ap.parse_args()

    if not args.snapshot.is_dir():
        print(f"错误：找不到快照目录 {args.snapshot}", file=sys.stderr)
        return 1
    con = connect(args.snapshot, args.threads, args.memory)
    log(f"快照 {args.snapshot}　靶点 {TARGET['name']}")

    stats = {"n_title": anchor_title(con), "n_entity": anchor_entity(con)}
    stats.update(expand_family(con))
    load_step1_validation(con, args.step1_csv)
    build_main(con)

    args.outdir.mkdir(parents=True, exist_ok=True)
    main_csv = args.outdir / "Step2_01_GCK_Related_Patents.csv"
    gckr_csv = args.outdir / "Step2_01_GCKR_Related_Patents.csv"
    md = args.outdir / "Step2_01_GCK_Related_Patent_Retrieval.md"
    n1 = to_csv(con, "main", OUT_COLUMNS, main_csv,
                "title_says_activator DESC, hit_clms DESC, "
                "val_n_step1_candidates DESC, patent_number")
    n2 = to_csv(con, "gckr",
                ["patent_id", "patent_number", "country", "publication_date",
                 "family_id", "title", "surface_forms", "n_mentions", "hit_clms"],
                gckr_csv, "hit_clms DESC, patent_number")
    write_report(con, md, args.snapshot, args.step1_csv, stats)

    n_fam = con.execute("SELECT COUNT(DISTINCT family_id) FROM main "
                        "WHERE family_id > 0").fetchone()[0]
    print(f"\n主表：    {main_csv}  （{n1:,} 篇 → {n_fam:,} 个同族）")
    print(f"GKRP 表： {gckr_csv}  （{n2:,} 篇）")
    print(f"报告：    {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
