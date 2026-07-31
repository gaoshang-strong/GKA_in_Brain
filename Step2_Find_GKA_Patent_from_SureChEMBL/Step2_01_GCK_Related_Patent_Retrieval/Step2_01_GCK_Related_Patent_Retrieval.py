#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step2_01_GCK_Related_Patent_Retrieval
=====================================

从 SureChEMBL bulk 快照检索**所有与葡萄糖激酶（GCK）相关的专利**，作为 Step2 的锚定表。

**召回优先**：不判断方向（激活/抑制）、不判断是不是 GKA 专利、不做 field 过滤。
**纯规则检索**：只用 bulk parquet 的实体标注，不调 API、不用 LLM。

锚定：为什么按 resolved_form 而不是文本匹配
------------------------------------------
专利全文来自扫描件 OCR，`glucokinase` 在 SureChEMBL 里有 33 种字面写法——
`glucokmase`、`giucokinase`、`gl uc ok i na se`、`Glucokinas e`、`Hexokinase 4`……
按 `original_text` 做 `LIKE '%glucokinase%'` 会漏掉其中一大半。
`resolved_form = 'HGNC:4195'` 一次兜住。

但 `resolved_form` 之外还有两类必须手工处理的边界（都在库上实测过，见 `ANCHOR_*` 常量）：

1. **未解析但明确是 GCK** → 白名单补进来（如 `glucokinase (hexokinase 4)`；
   MODY2 就是 GCK，所以那条也算）。
2. **名字带 glucokinase / hexokinase 但不是人 GCK** → 明确排除，理由写进报告。
   包括 ADPGK、细菌 ATP/多聚磷酸 glucokinase、HK1/2/3、酮己糖激酶、假基因，
   以及**泛称 `hexokinase`（29,158 篇，规模与锚点集相当但无法分辨，单独统计）**。

GKRP（`HGNC:4196` / `Q14397` 等）**另出一张表**——与 ChEMBL 侧把
`CHEMBL3885579`(GCK–GKRP PPI) 单列一致，GKRP 解离剂与直接激活 GCK 是两类机制。

输出
----
    Step2_01_GCK_Related_Patents.csv    一行一篇专利文档（主产物）
    Step2_01_GCKR_Related_Patents.csv   GKRP 相关，单独一张
    Step2_01_GCK_Related_Patent_Retrieval.md   报告

用法
----
    python3 Step2_01_GCK_Related_Patent_Retrieval.py
    python3 Step2_01_GCK_Related_Patent_Retrieval.py --snapshot /path/to/SureChEMBL_2026-07-17
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_SNAP = REPO / "SureChEMBL" / "SureChEMBL_2026-07-17"
DEFAULT_STEP1 = (REPO / "Step1_Find_GKA_from_ChEMBL"
                 / "Step1_06_GKA_Physicochemical_Property_Extraction"
                 / "Step1_06_GKA_Physicochemical_Properties.csv")

# ---------------------------------------------------------------------------
# 锚定规则（每一条都在库上实测过，文档数记在 Readme 与报告里）
# ---------------------------------------------------------------------------

# 核心锚点：人 GCK 基因
ANCHOR_RESOLVED = ["HGNC:4195"]

# 未解析但明确是人 GCK 的写法，白名单补召回
ANCHOR_TEXT_WHITELIST = [
    "glucokinase (hexokinase 4)",
    "Glucokinase (hexokinase 4, maturity onset diabetes of the young 2)",
]

# GKRP：不是 GCK，单独成表
GCKR_RESOLVED = ["HGNC:4196", "Q14397"]
GCKR_TEXT = [
    "glucokinase regulatory protein",
    "Glucokinase (hexokinase 4) regulator",
]

# 名字带 glucokinase/hexokinase 但不是人 GCK —— 明确排除，理由写进报告
EXCLUDED = [
    ("ADP-dependent glucokinase 及变体", ["ADP-dependent glucokinase", "ADP dependent glucokinase",
                                          "ADP -dependent glucokinase", "ADP- dependent glucokinase"],
     ["HGNC:25250"], "ADPGK，另一个酶"),
    ("ATP-dependent glucokinase 及变体", ["ATP-dependent glucokinase", "ATP dependent glucokinase",
                                          "ATP- dependent glucokinase"],
     [], "细菌/古菌的酶"),
    ("polyphosphate glucokinase", ["polyphosphate glucokinase"], [], "细菌的酶"),
    ("glucokinase 1 / -1", ["glucokinase 1", "glucokinase-1"], ["Q9GTW9"], "非人物种"),
    ("glucokinase-associated dual specificity phosphatase",
     ["glucokinase-associated dual specificity phosphatase"], ["Q9JIM4"], "另一个蛋白"),
    ("glucokinase 假基因", ["glucokinase activity, related sequence 1",
                             "glucokinase activity, related sequence 2",
                             "Glucokinase-Like"], [], "假基因/类似序列"),
    ("己糖激酶家族其他成员", [], ["HGNC:4922", "HGNC:4923", "HGNC:4925", "HGNC:6315"],
     "HK1 / HK2 / HK3 / 酮己糖激酶，不是 GCK"),
]

# 泛称 hexokinase：规模与锚点集相当但无法分辨，单独统计不进主表
AMBIGUOUS_GENERIC = "hexokinase"

# 会带进假阳性的 surface form —— 标记不删除
RISK_FORMS = {
    "GK": "⚠ 糖尿病文献里更常指 Goto-Kakizaki 大鼠（2 型糖尿病模型），与本领域高度重叠",
    "4": "⚠ 单个数字被解析成基因，标注错误",
    "glk": "细菌 glucokinase 基因名",
    "GlkA": "细菌 glucokinase 基因名",
    "gukA": "细菌 glucokinase 基因名",
    "GlcK": "细菌 glucokinase 基因名",
    "Hk4": "缩写，需核",
    "gki": "缩写，需核",
    "GluK": "缩写，需核",
}
# 正向信号：这两个其实是 "glucokinase activator" 的缩写
POSITIVE_FORMS = {"GKA": "「glucokinase activator」的缩写，方向正向信号",
                  "GKAs": "同上"}

# 阳性对照（后两个是 Step1 之后补的，见 Readme）
EXTRA_CONTROLS = {
    "OYUDYQMFVRHPIY-UHFFFAOYSA-N": ("BMS-820132", "CHEMBL5072532"),
    "HMUMWSORCUWQJO-QAPCUYQASA-N": ("DORZAGLIATIN", "CHEMBL4297508"),
}

OUT_COLUMNS = [
    "patent_id", "patent_number", "country", "publication_date",
    "family_id", "family_valid", "title", "assignee",
    "n_surface_forms", "surface_forms", "matched_by",
    "hit_ttl", "hit_abst", "hit_desc", "hit_clms", "n_mentions",
    "risk_flags", "has_activator_abbrev",
    "n_compounds_total", "n_compounds_clms",
    "n_step1_candidates", "n_step1_candidates_clms", "step1_candidate_ids",
    "is_positive_control_patent", "positive_control_names",
]

T0 = time.time()


def log(m: str) -> None:
    print(f"[{time.time() - T0:6.1f}s] {m}", file=sys.stderr, flush=True)


def fmt(n) -> str:
    if n is None:
        return "—"
    if isinstance(n, float) and abs(n - round(n)) < 1e-9:
        n = int(round(n))
    return f"{n:,}" if isinstance(n, int) else str(n)


def pct(a, b) -> str:
    return "—" if not b else f"{100.0 * a / b:.1f}%"


def sql_list(xs) -> str:
    return ",".join("'" + str(x).replace("'", "''") + "'" for x in xs)


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


def build_anchor(con) -> dict:
    """建锚点实体表，并把每一条纳入/排除的理由记下来备写报告。"""
    con.execute(f"""
        CREATE OR REPLACE TABLE anchor_entity AS
        SELECT id, original_text,
               CASE WHEN resolved_form IN ({sql_list(ANCHOR_RESOLVED)})
                    THEN 'resolved' ELSE 'whitelist' END AS matched_by
        FROM biomedical_entities
        WHERE resolved_form IN ({sql_list(ANCHOR_RESOLVED)})
           OR original_text IN ({sql_list(ANCHOR_TEXT_WHITELIST)})
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE gckr_entity AS
        SELECT id, original_text FROM biomedical_entities
        WHERE resolved_form IN ({sql_list(GCKR_RESOLVED)})
           OR original_text IN ({sql_list(GCKR_TEXT)})
    """)
    n_anchor = con.execute("SELECT COUNT(*) FROM anchor_entity").fetchone()[0]
    n_gckr = con.execute("SELECT COUNT(*) FROM gckr_entity").fetchone()[0]
    log(f"锚点实体 {n_anchor} 个，GKRP 实体 {n_gckr} 个")
    return {"n_anchor": n_anchor, "n_gckr": n_gckr}


def load_step1(con, path: Path) -> int:
    """读 Step1_06 的 782 个候选，按 InChIKey 对齐到 SureChEMBL compound_id。"""
    if not path.is_file():
        log(f"警告：找不到 {path}，跳过与 Step1 的交叉")
        con.execute("CREATE OR REPLACE TABLE step1_cmp(cid BIGINT, chembl_id VARCHAR, "
                    "ctl_name VARCHAR)")
        return 0
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    con.execute("CREATE OR REPLACE TABLE step1_raw(ik VARCHAR, chembl_id VARCHAR, "
                "ctl_name VARCHAR)")
    data = [(r["standard_inchi_key"], r["molecule_chembl_id"],
             r["control_name"] if r.get("is_positive_control") == "TRUE" else "")
            for r in rows if r.get("standard_inchi_key")]
    for ik, (name, cid) in EXTRA_CONTROLS.items():
        data.append((ik, cid, name))
    con.executemany("INSERT INTO step1_raw VALUES (?,?,?)", data)
    # SureChEMBL 的 inchi_key 不唯一，先 DISTINCT 再 join，否则会放大行数
    con.execute("""
        CREATE OR REPLACE TABLE step1_cmp AS
        SELECT DISTINCT c.id AS cid, s.chembl_id, s.ctl_name
        FROM step1_raw s JOIN compounds c ON c.inchi_key = s.ik
    """)
    n = con.execute("SELECT COUNT(*) FROM step1_cmp").fetchone()[0]
    log(f"Step1 候选 {len(data)} 个（含 {len(EXTRA_CONTROLS)} 个补充对照），"
        f"在 SureChEMBL 里匹配到 {n} 个 compound_id")
    return n


def build_main(con) -> None:
    """主表：一行一篇专利文档。"""
    con.execute("""
        CREATE OR REPLACE TABLE hit AS
        SELECT l.patent_id,
               COUNT(DISTINCT e.id)                                   AS n_surface_forms,
               list_sort(list_distinct(list(e.original_text)))        AS surface_forms,
               list_sort(list_distinct(list(e.matched_by)))           AS matched_by,
               SUM(l.count)                                           AS n_mentions,
               SUM(l.count) FILTER (WHERE l.field_id = 4)             AS hit_ttl,
               SUM(l.count) FILTER (WHERE l.field_id = 3)             AS hit_abst,
               SUM(l.count) FILTER (WHERE l.field_id = 1)             AS hit_desc,
               SUM(l.count) FILTER (WHERE l.field_id = 2)             AS hit_clms
        FROM biomedical_locations l
        JOIN anchor_entity e ON e.id = l.entity_id
        GROUP BY 1
    """)
    log("锚点命中表建好")

    # 每篇专利的化合物数（全部 / 仅权利要求）
    con.execute("""
        CREATE OR REPLACE TABLE cmp_cnt AS
        SELECT m.patent_id,
               COUNT(DISTINCT m.compound_id)                                       AS n_compounds_total,
               COUNT(DISTINCT m.compound_id) FILTER (WHERE m.field_id = 2)         AS n_compounds_clms
        FROM patent_compound_map m
        JOIN (SELECT patent_id FROM hit) h ON h.patent_id = m.patent_id
        GROUP BY 1
    """)
    log("化合物计数建好")

    # 与 Step1 候选的交叉
    con.execute("""
        CREATE OR REPLACE TABLE s1_cnt AS
        SELECT m.patent_id,
               COUNT(DISTINCT s.chembl_id)                                   AS n_step1_candidates,
               COUNT(DISTINCT s.chembl_id) FILTER (WHERE m.field_id = 2)     AS n_step1_candidates_clms,
               list_sort(list_distinct(list(s.chembl_id)))                   AS step1_candidate_ids,
               list_sort(list_distinct(list(s.ctl_name) FILTER (WHERE s.ctl_name <> ''))) AS ctl_names
        FROM patent_compound_map m
        JOIN step1_cmp s ON s.cid = m.compound_id
        JOIN (SELECT patent_id FROM hit) h ON h.patent_id = m.patent_id
        GROUP BY 1
    """)
    log("Step1 交叉建好")

    con.execute(f"""
        CREATE OR REPLACE TABLE main AS
        SELECT p.id                                   AS patent_id,
               p.patent_number,
               p.country,
               p.publication_date,
               p.family_id,
               CASE WHEN p.family_id > 0 THEN 'TRUE' ELSE 'FALSE' END AS family_valid,
               p.title,
               p.assignee,
               h.n_surface_forms,
               h.surface_forms,
               h.matched_by,
               COALESCE(h.hit_ttl, 0)  AS hit_ttl,
               COALESCE(h.hit_abst, 0) AS hit_abst,
               COALESCE(h.hit_desc, 0) AS hit_desc,
               COALESCE(h.hit_clms, 0) AS hit_clms,
               h.n_mentions,
               list_filter(h.surface_forms, x -> x IN ({sql_list(RISK_FORMS)})) AS risk_flags,
               CASE WHEN len(list_filter(h.surface_forms,
                                         x -> x IN ({sql_list(POSITIVE_FORMS)}))) > 0
                    THEN 'TRUE' ELSE 'FALSE' END AS has_activator_abbrev,
               COALESCE(c.n_compounds_total, 0) AS n_compounds_total,
               COALESCE(c.n_compounds_clms, 0)  AS n_compounds_clms,
               COALESCE(s.n_step1_candidates, 0) AS n_step1_candidates,
               COALESCE(s.n_step1_candidates_clms, 0) AS n_step1_candidates_clms,
               s.step1_candidate_ids,
               CASE WHEN len(COALESCE(s.ctl_names, [])) > 0 THEN 'TRUE' ELSE 'FALSE' END
                    AS is_positive_control_patent,
               s.ctl_names AS positive_control_names
        FROM hit h
        JOIN patents p ON p.id = h.patent_id
        LEFT JOIN cmp_cnt c ON c.patent_id = h.patent_id
        LEFT JOIN s1_cnt  s ON s.patent_id = h.patent_id
    """)
    n = con.execute("SELECT COUNT(*) FROM main").fetchone()[0]
    log(f"主表 {n:,} 行")

    con.execute("""
        CREATE OR REPLACE TABLE gckr AS
        SELECT p.id AS patent_id, p.patent_number, p.country, p.publication_date,
               p.family_id, p.title,
               list_sort(list_distinct(list(e.original_text))) AS surface_forms,
               SUM(l.count) AS n_mentions,
               SUM(l.count) FILTER (WHERE l.field_id = 2) AS hit_clms
        FROM biomedical_locations l
        JOIN gckr_entity e ON e.id = l.entity_id
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


def write_report(con, path: Path, snap: Path, step1_csv: Path, stats: dict) -> None:
    L = []
    q = lambda s: con.execute(s).fetchall()          # noqa: E731
    one = lambda s: con.execute(s).fetchone()[0]     # noqa: E731

    n_docs = one("SELECT COUNT(*) FROM main")
    n_fam = one("SELECT COUNT(DISTINCT family_id) FROM main WHERE family_id > 0")
    n_nofam = one("SELECT COUNT(*) FROM main WHERE family_id IS NULL OR family_id <= 0")

    L.append("# Step2_01 GCK 相关专利检索")
    L.append("")
    L.append("- 数据源：**SureChEMBL 2.0，`2026-07-17` 全量快照**")
    L.append(f"- 快照路径：`{snap}`")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 交叉引用：`{step1_csv.name}`（Step1_06 的 782 个候选）")
    L.append(f"- 命中专利文档：**{fmt(n_docs)}** 篇 → 同族去重后 **{fmt(n_fam)}** 个发明")
    L.append("")
    L.append("> **召回优先**：本步骤不判断方向（激活/抑制）、不判断是不是 GKA 专利、"
             "**不做 `field_id` 过滤**。纯规则检索，未调 API、未用 LLM。")
    L.append("")

    L.append("## 一、锚定规则与实测规模")
    L.append("")
    L.append(f"锚点实体 **{stats['n_anchor']}** 个。核心是 `resolved_form = 'HGNC:4195'`，"
             "**不做文本匹配**——专利全文来自 OCR，`glucokinase` 有几十种破碎写法。")
    L.append("")
    rows = q("""SELECT matched_by, COUNT(*) FROM anchor_entity GROUP BY 1 ORDER BY 2 DESC""")
    L.append("| 纳入方式 | 实体数 | 说明 |")
    L.append("| --- | ---: | --- |")
    for k, c in rows:
        note = ("`resolved_form = 'HGNC:4195'`，覆盖全部 OCR 变体与缩写" if k == "resolved"
                else "未解析但明确是 GCK 的写法，手工白名单")
        L.append(f"| `{k}` | {fmt(c)} | {note} |")
    L.append("")
    L.append("白名单逐条（未解析，人工确认是 GCK）：")
    L.append("")
    for t in ANCHOR_TEXT_WHITELIST:
        L.append(f"- `{t}`" + ("　← MODY2 就是 GCK" if "maturity onset" in t else ""))
    L.append("")

    L.append("### 明确排除的（名字带 glucokinase / hexokinase 但不是人 GCK）")
    L.append("")
    L.append("| 类别 | 判据 | 涉及专利 | 排除理由 |")
    L.append("| --- | --- | ---: | --- |")
    for label, texts, resolved, why in EXCLUDED:
        cond = []
        if texts:
            cond.append(f"original_text IN ({sql_list(texts)})")
        if resolved:
            cond.append(f"resolved_form IN ({sql_list(resolved)})")
        n = one(f"""SELECT COUNT(DISTINCT l.patent_id) FROM biomedical_locations l
                    JOIN biomedical_entities e ON e.id = l.entity_id
                    WHERE {' OR '.join(cond)}""")
        key = (", ".join(f"`{x}`" for x in (resolved or texts[:2]))
               + ("…" if len(texts) > 2 and not resolved else ""))
        L.append(f"| {label} | {key} | {fmt(n)} | {why} |")
    n_generic = one(f"""SELECT COUNT(DISTINCT l.patent_id) FROM biomedical_locations l
                        JOIN biomedical_entities e ON e.id = l.entity_id
                        WHERE e.original_text = '{AMBIGUOUS_GENERIC}'
                          AND COALESCE(e.resolved_form,'') = ''""")
    L.append(f"| **泛称 `hexokinase`** | 未解析 | **{fmt(n_generic)}** | "
             "规模与锚点集相当，但绝大多数是 HK1/HK2（肿瘤代谢），**无法分辨**，不进主表 |")
    L.append("")

    L.append("## 二、规模与去重")
    L.append("")
    L.append(f"**{fmt(n_docs)} 篇专利文档 → {fmt(n_fam)} 个同族**"
             f"（平均每族 {n_docs / n_fam:.2f} 篇），"
             f"另有 {fmt(n_nofam)} 篇没有有效 `family_id`（哨兵 `-1` 或 `NULL`）。")
    L.append("")
    L.append("> 去重必须 `COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)`——"
             "`-1` 是「未分配同族」的哨兵，不排除会把它们错当成同一个发明。")
    L.append("")

    L.append("## 三、⚠ 专利局覆盖严重不均（必读）")
    L.append("")
    rows = q("""SELECT country, COUNT(*) AS n_docs,
                       COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0) AS n_fam
                FROM main GROUP BY 1 ORDER BY 2 DESC""")
    total_by_country = dict(q("SELECT country, COUNT(*) FROM patents GROUP BY 1"))
    L.append("| 专利局 | 命中文档 | 同族 | 全库该国专利 | 命中率 |")
    L.append("| --- | ---: | ---: | ---: | ---: |")
    for c, nd, nf in rows:
        tot = total_by_country.get(c, 0)
        L.append(f"| `{c}` | {fmt(nd)} | {fmt(nf)} | {fmt(tot)} | {pct(nd, tot)} |")
    L.append("")
    L.append("**JP 与 CN 的命中率低到不能用**。原因不在检索式，在数据源本身：")
    L.append("")
    L.append("- **JPO 不提供全文**，SureChEMBL 只拿到著录项 + 英文标题摘要")
    L.append("- **CNIPA 只有英文机器翻译全文**，实体标注管道在机翻文本上效果差")
    L.append("")
    L.append("> **这批数据实质上是 US / EP / WO 的视图。**"
             "不能据此谈「全球 GKA 专利版图」，也不能说「中国/日本没有 GKA 专利」——"
             "是看不见，不是没有。")
    L.append("")

    L.append("## 四、命中位置分布")
    L.append("")
    rows = [("`ttl` 标题", one("SELECT COUNT(*) FROM main WHERE hit_ttl > 0")),
            ("`abst` 摘要", one("SELECT COUNT(*) FROM main WHERE hit_abst > 0")),
            ("`desc` 说明书", one("SELECT COUNT(*) FROM main WHERE hit_desc > 0")),
            ("**`clms` 权利要求**", one("SELECT COUNT(*) FROM main WHERE hit_clms > 0"))]
    L.append("| 位置 | 命中专利数 | 占比 |")
    L.append("| --- | ---: | ---: |")
    for k, c in rows:
        L.append(f"| {k} | {fmt(c)} | {pct(c, n_docs)} |")
    L.append("")
    L.append("本步骤**不按位置过滤**，四个位置的命中次数都写进主表的 "
             "`hit_ttl` / `hit_abst` / `hit_desc` / `hit_clms` 列，下游精筛直接用。"
             "`clms` 是精度轴、`desc` 是召回轴，别在锚定步骤就二选一。")
    L.append("")

    L.append("## 五、surface form 分解与风险标记")
    L.append("")
    rows = q("""SELECT f AS form, COUNT(*) AS n FROM (
                  SELECT UNNEST(surface_forms) AS f FROM main)
                GROUP BY 1 ORDER BY 2 DESC LIMIT 20""")
    L.append("| surface form | 命中专利数 | 备注 |")
    L.append("| --- | ---: | --- |")
    for f, c in rows:
        note = RISK_FORMS.get(f, POSITIVE_FORMS.get(f, ""))
        L.append(f"| `{f}` | {fmt(c)} | {note} |")
    L.append("")
    n_risk = one("SELECT COUNT(*) FROM main WHERE len(risk_flags) > 0")
    n_pos = one("SELECT COUNT(*) FROM main WHERE has_activator_abbrev = 'TRUE'")
    n_only_risk = one("""SELECT COUNT(*) FROM main
                         WHERE len(risk_flags) > 0 AND len(surface_forms) = len(risk_flags)""")
    L.append(f"- 带风险标记的专利 **{fmt(n_risk)}** 篇（`risk_flags` 列）")
    L.append(f"- **只靠风险形命中、没有任何可靠写法** 的 **{fmt(n_only_risk)}** 篇 "
             "← **这批最可疑，下游应优先人工核或直接排除**")
    L.append(f"- 出现 `GKA` / `GKAs` 缩写的 **{fmt(n_pos)}** 篇（`has_activator_abbrev = TRUE`）"
             "——这两个缩写本身就是「glucokinase activator」，是**方向正向信号**，"
             "虽然 SureChEMBL 把它解析成了基因")
    L.append("")

    L.append("## 六、与 Step1（ChEMBL 侧）的交叉")
    L.append("")
    n_with = one("SELECT COUNT(*) FROM main WHERE n_step1_candidates > 0")
    n_with_clms = one("SELECT COUNT(*) FROM main WHERE n_step1_candidates_clms > 0")
    L.append(f"用 **InChIKey** 对齐（唯一的跨库桥梁）。"
             f"{fmt(stats['n_step1_cmp'])} 个 Step1 候选在 SureChEMBL 里找到了对应结构。")
    L.append("")
    L.append("| 项目 | 专利数 |")
    L.append("| --- | ---: |")
    L.append(f"| 含至少 1 个 Step1 候选化合物 | {fmt(n_with)} |")
    L.append(f"| 权利要求里含 Step1 候选化合物 | {fmt(n_with_clms)} |")
    L.append("")
    L.append("> ⚠ SureChEMBL 的 `inchi_key` 不唯一（30,990,818 行 / 29,874,136 唯一），"
             "join 前必须 `DISTINCT`，否则行数会被放大。本脚本已处理。")
    L.append("")

    L.append("## 七、阳性对照自检")
    L.append("")
    L.append("已知 GKA 的专利必须能被这套规则捞到。落不进说明规则有问题，不是数据有问题。")
    L.append("")
    rows = q("""SELECT c AS name,
                       COUNT(*) AS n_docs,
                       COUNT(*) FILTER (WHERE hit_clms > 0) AS n_clms
                FROM (SELECT UNNEST(positive_control_names) AS c, hit_clms FROM main)
                GROUP BY 1 ORDER BY 2 DESC""")
    found = {r[0] for r in rows}
    L.append("| 对照 | 命中专利数 | 其中权利要求也提到 GCK |")
    L.append("| --- | ---: | ---: |")
    for name, nd, nc in rows:
        L.append(f"| {name} | {fmt(nd)} | {fmt(nc)} |")
    all_ctl = set(one("SELECT list(DISTINCT ctl_name) FROM step1_cmp WHERE ctl_name <> ''") or [])
    for miss in sorted(all_ctl - found):
        L.append(f"| {miss} | **0** | 0 |")
    L.append("")
    n_ok, n_all = len(found), len(all_ctl)
    L.append(f"**自检结论：{n_ok}/{n_all} 个对照的专利被捞到。**")
    if n_ok < n_all:
        L.append("")
        L.append("未命中的对照说明：该化合物出现的专利里，GCK 没有被实体标注管道识别出来。"
                 "**这是召回缺口的直接证据**——纯实体锚定会漏掉一部分真 GKA 专利，"
                 "下游若要补，只能靠结构相似性或全文检索。")
    L.append("")

    L.append("## 八、GKRP 单独成表")
    L.append("")
    ng = one("SELECT COUNT(*) FROM gckr")
    ngf = one("SELECT COUNT(DISTINCT family_id) FROM gckr WHERE family_id > 0")
    ov = one("""SELECT COUNT(*) FROM (SELECT patent_id FROM gckr
                INTERSECT SELECT patent_id FROM main)""")
    L.append(f"**{fmt(ng)}** 篇（{fmt(ngf)} 个同族），与主表重叠 **{fmt(ov)}** 篇。")
    L.append("")
    L.append("GKRP（`HGNC:4196` / `Q14397` / `glucokinase regulatory protein` 等）"
             "**不并入主表**——GKRP 解离剂与直接激活 GCK 是两类机制，"
             "与 ChEMBL 侧把 `CHEMBL3885579`(GCK–GKRP PPI) 单列的处理一致。"
             "重叠的那部分两张表都有，按需取用。")
    L.append("")

    L.append("## 九、这一步没做什么")
    L.append("")
    L.append("| 没做 | 为什么 | 留给谁 |")
    L.append("| --- | --- | --- |")
    L.append("| 方向判定（激活 / 抑制） | bulk 数据无全文，`Mechanism` 实体类型是工业化学词，"
             "`resolved_form` 全空，判不了 | 后续步骤（读权利要求原文） |")
    L.append("| `field_id` 过滤 | 锚定步骤要召回，`clms` 与 `desc` 都留着 | 下游精筛，用主表的 hit_* 列 |")
    L.append("| 结构相似性检索 | 本步骤只做实体锚定 | 后续用 `fpsim2_fingerprints.h5` |")
    L.append("| 泛称 `hexokinase` 的那 %s 篇 | 无法分辨 HK1/2/3 与 GCK | 若要补召回，需结构或全文佐证 |"
             % fmt(n_generic))
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="从 SureChEMBL 检索 GCK 相关专利（纯规则）。")
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAP)
    ap.add_argument("--step1-csv", type=Path, default=DEFAULT_STEP1)
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="12GB")
    args = ap.parse_args()

    if not args.snapshot.is_dir():
        print(f"错误：找不到快照目录 {args.snapshot}", file=sys.stderr)
        return 1

    con = connect(args.snapshot, args.threads, args.memory)
    log(f"快照 {args.snapshot}")
    stats = build_anchor(con)
    stats["n_step1_cmp"] = load_step1(con, args.step1_csv)
    build_main(con)

    args.outdir.mkdir(parents=True, exist_ok=True)
    main_csv = args.outdir / "Step2_01_GCK_Related_Patents.csv"
    gckr_csv = args.outdir / "Step2_01_GCKR_Related_Patents.csv"
    md = args.outdir / "Step2_01_GCK_Related_Patent_Retrieval.md"

    n1 = to_csv(con, "main", OUT_COLUMNS, main_csv,
                "hit_clms DESC, n_step1_candidates DESC, patent_number")
    n2 = to_csv(con, "gckr",
                ["patent_id", "patent_number", "country", "publication_date",
                 "family_id", "title", "surface_forms", "n_mentions", "hit_clms"],
                gckr_csv, "hit_clms DESC, patent_number")
    write_report(con, md, args.snapshot, args.step1_csv, stats)

    n_fam = con.execute("SELECT COUNT(DISTINCT family_id) FROM main "
                        "WHERE family_id > 0").fetchone()[0]
    print(f"\n主表：      {main_csv}  （{n1:,} 篇 → {n_fam:,} 个同族）")
    print(f"GKRP 表：   {gckr_csv}  （{n2:,} 篇）")
    print(f"报告：      {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
