#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
surechembl_profile.py — 调查一份 SureChEMBL 2.0 bulk parquet 快照，生成 Markdown 概览报告。

设计目标
--------
与同项目的 `ChEMBL/chembl_profile.py` 对齐：报告的读者假定为
**有生物学 / 生信背景、但没有专利检索经验**的人。因此不只给数字，
还解释每个概念是什么、为什么这么设计、怎么用、哪里会踩坑。

用法
----
    python3 surechembl_profile.py SureChEMBL_2026-07-17/
    python3 surechembl_profile.py SureChEMBL_2026-07-17/ -o report.md
    python3 surechembl_profile.py SureChEMBL_2026-07-17/ --quick   # 只用 parquet 元数据，秒级
    python3 surechembl_profile.py SureChEMBL_2026-07-17/ --deep    # 加做分位数/完整性全扫（约 2–3 分钟）

只读打开，不会改动任何文件。用 duckdb 直接查 parquet——
`patent_compound_map` 有 15.4 亿行，**任何情况下都不要整表读进内存**。
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime
from pathlib import Path

try:
    import duckdb
except ImportError:
    sys.exit("需要 duckdb：micromamba install -n GKA_in_Brain -c conda-forge duckdb")

# ---------------------------------------------------------------------------
# 知识库：把专利库的行话翻译成生物学家能懂的话
# ---------------------------------------------------------------------------

TABLE_DESC = {
    "patents": "专利主表：一行一篇专利文档，含专利号、国别、公开日、同族、分类号、申请人、标题",
    "compounds": "化合物主表：一行一个唯一化学结构（SMILES / InChI / InChIKey / 分子量）",
    "patent_compound_map": "核心事实表：某化合物出现在某专利的某个部分。**全库最大，15.4 亿行**",
    "fields": "字典：专利被切成哪几个部分（标题/摘要/说明书/权利要求/图片/MOL 附件）",
    "biomedical_entities": "生物医学实体词表：文中出现的基因/蛋白/疾病/机制词，及其归一化 ID",
    "biomedical_types": "字典：实体的四种类型",
    "biomedical_locations": "实体位置表：某实体出现在某专利某部分，及出现次数",
    "fpsim2_fingerprints.h5": "FPSim2 指纹库（HDF5，非 parquet），可直接做全库相似性检索",
}

# 速览表用的短名
SHORT_NAME = {
    "patents": "专利数",
    "compounds": "化合物数",
    "patent_compound_map": "专利-化合物关联数",
    "biomedical_locations": "实体位置记录数",
    "biomedical_entities": "生物医学实体数",
}

FIELD_NOTE = {
    "ttl": "标题。信息密度最高但覆盖最少",
    "abst": "摘要",
    "desc": "说明书正文。**包含背景技术**，所以这里出现的化合物很多是他人的、被引用的",
    "clms": "**权利要求**。这是专利真正主张保护的范围——判断「这篇专利要保护什么化合物」只能看这里",
    "image": "从图片里识别出的结构（化学结构图 OSR）",
    "molattachment": "专利附带的 MOL 文件（序列表/结构文件）",
}

COUNTRY_NOTE = {
    "CN": "中国国家知识产权局（CNIPA）。SureChEMBL 2.0 新增，专利数最多但化学信息密度最低",
    "US": "美国专利商标局（USPTO）。化合物贡献最大",
    "EP": "欧洲专利局（EPO）",
    "JP": "日本特许厅（JPO）",
    "WO": "世界知识产权组织（WIPO）的 PCT 国际申请。**不是某国专利**，是进入各国前的国际阶段申请",
}

TYPE_NOTE = {
    "GeneOrProtein": "基因/蛋白名，归一到 UniProt、HGNC、Entrez Gene",
    "Disease": "疾病名，归一到 MeSH、Disease Ontology、Wikipedia",
    "Mechanism": "化学物质的作用或作用机制",
    "Physquant": "各类物理量",
}

CLASSIFICATION_NOTE = {
    "cpc": "Cooperative Patent Classification，EPO 与 USPTO 联合分类体系，最细（约 25 万个条目）",
    "ipc": "International Patent Classification，WIPO 的国际分类，粗一些",
    "ipcr": "IPC 的 reformed 版本（2006 年后的写法）",
    "ecla": "European Classification，EPO 的旧体系，**已被 CPC 取代**，新专利基本没有",
}

# 与 GKA 项目相关的分类号，报告里用来举例
CPC_OF_INTEREST = {
    "A61P3/10": "治疗糖尿病/高血糖的药物用途",
    "C07D": "杂环化合物（绝大多数小分子药物落在这里）",
    "A61K31": "含有机活性成分的医药制剂",
}

# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:6.1f}s] {msg}", file=sys.stderr, flush=True)


def fmt(n) -> str:
    if n is None:
        return "—"
    if isinstance(n, float):
        if n != n:            # NaN
            return "—"
        if abs(n - round(n)) < 1e-9:
            n = int(round(n))
        else:
            return f"{n:,.2f}"
    return f"{n:,}"


def pct(part, whole) -> str:
    if not whole:
        return "—"
    return f"{100.0 * part / whole:.1f}%"


def human_bytes(n: int) -> str:
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.2f} {unit}" if unit not in ("B", "KB") else f"{x:.0f} {unit}"
        x /= 1024
    return f"{n} B"


def bar(value: float, vmax: float, width: int = 26) -> str:
    if not vmax or value is None:
        return ""
    return "█" * max(1, int(round(width * value / vmax))) if value > 0 else ""


class Snap:
    """一份 SureChEMBL bulk 快照。所有查询走 duckdb，绝不整表加载。"""

    TABLES = ["patents", "compounds", "patent_compound_map", "fields",
              "biomedical_entities", "biomedical_types", "biomedical_locations"]

    def __init__(self, root: Path, threads: int, mem: str):
        self.root = root.resolve()
        self.con = duckdb.connect()
        self.con.execute(f"SET threads={threads}; SET memory_limit='{mem}';")
        # 进度条会把 stdout 刷爆，必须关
        self.con.execute("SET enable_progress_bar=false;")
        self.files = {}
        for t in self.TABLES:
            p = self.root / f"{t}.parquet"
            if p.is_file():
                self.files[t] = p
        self.fps = self.root / "fpsim2_fingerprints.h5"

    def has(self, t: str) -> bool:
        return t in self.files

    def p(self, t: str) -> str:
        """给 duckdb 用的 parquet 路径字面量。"""
        return "'" + str(self.files[t]) + "'"

    def q(self, sql: str, label: str = "") -> list:
        s = time.time()
        try:
            rows = self.con.sql(sql).fetchall()
        except Exception as e:                       # noqa: BLE001
            log(f"查询失败（{label}）：{e}")
            return []
        if label:
            log(f"{label}  {time.time() - s:.1f}s")
        return rows

    def one(self, sql: str, label: str = ""):
        r = self.q(sql, label)
        return r[0][0] if r and r[0] else None

    def meta(self, t: str) -> dict:
        """只读 parquet footer，不扫数据。"""
        r = self.q(f"SELECT num_rows, num_row_groups FROM parquet_file_metadata({self.p(t)})")
        n_rows, n_rg = (r[0] if r else (None, None))
        cols = self.q(f"SELECT name, type FROM parquet_schema({self.p(t)}) WHERE num_children IS NULL")
        return {"rows": n_rows, "row_groups": n_rg, "cols": cols,
                "bytes": self.files[t].stat().st_size}


class Report:
    def __init__(self):
        self.L = []

    def h(self, level: int, text: str) -> None:
        self.L.append("")
        self.L.append("#" * level + " " + text)
        self.L.append("")

    def p(self, text: str = "") -> None:
        self.L.append(text)

    def table(self, header: list, rows: list, align: list = None) -> None:
        self.L.append("| " + " | ".join(header) + " |")
        if align:
            self.L.append("| " + " | ".join(align) + " |")
        else:
            self.L.append("| " + " | ".join("---" for _ in header) + " |")
        for r in rows:
            self.L.append("| " + " | ".join("" if c is None else str(c) for c in r) + " |")
        self.L.append("")

    def code(self, text: str, lang: str = "") -> None:
        self.L.append("```" + lang)
        self.L.extend(text.rstrip("\n").split("\n"))
        self.L.append("```")
        self.L.append("")

    def write(self, path: Path) -> None:
        path.write_text("\n".join(self.L).rstrip() + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# 各章节
# ---------------------------------------------------------------------------

def sec_intro(r: Report, s: Snap, args) -> dict:
    r.p("# SureChEMBL 专利化学数据库结构与内容概览报告")
    r.p("")
    r.p(f"> 自动生成于 {datetime.now():%Y-%m-%d %H:%M}　|　快照目录：`{s.root}`")
    r.p(">")
    r.p("> 本报告面向**有生物学 / 生信背景、但没有专利检索经验**的读者。"
        "每一节先解释「这是什么」，再给出这份快照里的实际统计。")
    r.p("")

    r.h(2, "0. 三十秒速览")
    r.p("**SureChEMBL 是什么？** EMBL-EBI 维护的**专利化学数据库**。"
        "它把专利全文（含扫描件 OCR、结构图识别）里的化学结构和生物医学名词自动抽取出来，"
        "变成可查询的表。与 ChEMBL 的关键区别：")
    r.p("")
    r.table(
        ["", "ChEMBL", "SureChEMBL"],
        [["数据来源", "论文 + 部分专利", "**专利全文**"],
         ["加工方式", "**人工审编**", "**全自动抽取**（NER + OSR + OCR）"],
         ["核心事实", "化合物 × 靶点 × 活性数值", "化合物 × 专利 × **出现位置**"],
         ["有没有活性数据", "有（IC50/EC50…）", "**没有**。只知道「这个化合物出现在这篇专利里」"],
         ["规模", "~292 万化合物", "~3,099 万化合物"]])
    r.p("**最重要的一条**：SureChEMBL **不告诉你化合物有没有活性**，"
        "只告诉你它出现在哪篇专利的哪个部分。活性要回 ChEMBL 或读专利原文。")
    r.p("")
    r.p("**一句话的数据模型**：")
    r.code("某个化合物  出现在某篇专利  的某个部分\n"
           "(compounds)   (patents)       (fields)\n"
           "         └────────┬────────────────┘\n"
           "            patent_compound_map\n"
           "\n"
           "某个生物医学实体  出现在某篇专利  的某个部分  出现了几次\n"
           "(biomedical_entities) (patents)   (fields)    (count)\n"
           "         └────────────┬─────────────────────────┘\n"
           "                biomedical_locations")

    info = {}
    rows = []
    total_bytes = 0
    for t in Snap.TABLES:
        if not s.has(t):
            continue
        m = s.meta(t)
        info[t] = m
        total_bytes += m["bytes"]
    if s.fps.is_file():
        total_bytes += s.fps.stat().st_size

    rows.append(["快照目录", f"`{s.root.name}`"])
    rows.append(["文件总大小", human_bytes(total_bytes)])
    for t in ("patents", "compounds", "patent_compound_map",
              "biomedical_locations", "biomedical_entities"):
        if t in info:
            rows.append([f"{SHORT_NAME[t]}（`{t}`）", fmt(info[t]["rows"])])
    rows.append(["duckdb 版本", duckdb.__version__])
    r.table(["条目", "值"], rows)

    r.p("> ⚠ **这份快照是某一期的全量镜像，不是增量。** SureChEMBL 每两周发一版并覆盖 "
        "`latest/`，因此本项目固定在日期目录上。任何数字都对应这一版。")
    r.p("")
    return info


def sec_glossary(r: Report, s: Snap) -> None:
    r.h(2, "1. 先搞懂这些词")

    r.h(3, "1.1 专利文档（patent document）≠ 发明")
    r.p("同一项发明会在多个国家/多个阶段分别公开，每次公开都是**一篇独立的专利文档**，"
        "有各自的专利号。所以「专利数」远大于「发明数」。")
    r.p("")
    r.p("把同一发明的所有文档串起来的字段是 **`family_id`（专利同族）**。"
        "**做去重统计必须按 `family_id`，不能按 `patent_number`**，"
        "否则一个化合物会因为在 5 个国家申请而被数 5 次。")

    r.h(3, "1.2 `WO` 不是一个国家")
    r.p("`country` 里的 `WO` 是 WIPO 的 **PCT 国际申请**——申请人先交一份国际申请，"
        "之后再决定进入哪些国家。它和 `US`/`EP`/`JP`/`CN` 不是并列的地理关系，"
        "**同一发明常常既有 WO 文档又有各国文档**。")

    r.h(3, "1.3 专利被切成 6 个部分（`fields`）⭐")
    r.p("这是全库**最容易踩坑**的设计。一篇专利被切成标题/摘要/说明书/权利要求/图片/MOL 附件，"
        "`patent_compound_map.field_id` 指出化合物是在哪个部分被发现的。")
    r.p("")
    rows = []
    for fid, name in s.q(f"SELECT id, field_name FROM {s.p('fields')} ORDER BY id"):
        rows.append([fid, f"`{name}`", FIELD_NOTE.get(name, "")])
    r.table(["id", "field_name", "说明"], rows)
    r.p("**`desc`（说明书）里出现 ≠ 这篇专利要保护它。** 说明书含背景技术，"
        "会大段引用他人的化合物做对比。真正主张保护的范围只在 **`clms`（权利要求）**。")

    r.h(3, "1.4 实体解析（entity resolution）：`resolved_form` 才是锚点")
    r.p("`biomedical_entities` 里，`original_text` 是原文里出现的字面写法，"
        "`resolved_form` 是归一后的标识符。因为专利全文来自扫描件 OCR，"
        "**同一个概念会有几十种字面写法**（断字、错字、缩写）。")
    r.p("")
    r.p("**按 `original_text` 做字符串匹配一定会漏，必须按 `resolved_form` 锚定。**")
    r.p("")
    rows = []
    for tid, name, desc in s.q(
            f"SELECT id, type_name, description FROM {s.p('biomedical_types')} ORDER BY id"):
        n = s.one(f"SELECT COUNT(*) FROM {s.p('biomedical_entities')} WHERE type_id={tid}")
        rows.append([tid, f"`{name}`", TYPE_NOTE.get(name, desc or ""), fmt(n)])
    r.table(["id", "type_name", "说明", "本快照实体数"], rows,
            ["---:", "---", "---", "---:"])

    resolved = s.one(f"SELECT COUNT(*) FROM {s.p('biomedical_entities')} "
                     "WHERE resolved_form IS NOT NULL AND resolved_form <> ''")
    total = s.one(f"SELECT COUNT(*) FROM {s.p('biomedical_entities')}")
    r.p(f"**只有 {fmt(resolved)} / {fmt(total)}（{pct(resolved, total)}）的实体有 "
        "`resolved_form`**，其余是抽出来了但没能归一到任何标识符的词。"
        "未归一的实体**无法可靠地用于检索**——你不知道它到底指什么。")
    r.p("")
    r.p("⚠ `resolved_form` **不是单一命名空间**：同一列里混着 "
        "`HGNC:4195`（带前缀）、`Q14397`（裸 UniProt accession）、MeSH ID 等。"
        "解析前必须先判断命名空间。")

    r.h(3, "1.5 分类号：CPC / IPC / IPCR / ECLA")
    rows = [[f"`{k}`", v] for k, v in CLASSIFICATION_NOTE.items()]
    r.table(["字段", "说明"], rows)
    r.p("这四个字段在 `patents.parquet` 里都是 **list 类型**（一篇专利可以有多个分类号），"
        "**`LIKE` 匹配不上**，要用 `list_contains()` / `UNNEST()`。")


def sec_schema(r: Report, s: Snap, info: dict) -> None:
    r.h(2, "2. 数据模型与表清单")
    r.p("Schema 镜像 SureChEMBL 内部关系库，**不做冗余展开**——"
        "这是 2.0 相对旧版 MAP files 的主要改进，文件体积因此小很多。")
    r.p("")
    r.code("""compounds.id ──────┐
                   ├──< patent_compound_map >── fields.id
patents.id ────────┘                    (field_id)

biomedical_entities.id ──┐
                         ├──< biomedical_locations >── fields.id
patents.id ──────────────┘                    (field_id)
        │
        └── biomedical_entities.type_id ──> biomedical_types.id""")
    r.p("**没有外键约束**（parquet 不支持），完整性靠上游保证——"
        "本报告 §7 有实测的完整性检查。")
    r.p("")

    rows = []
    for t in Snap.TABLES:
        if t not in info:
            continue
        m = info[t]
        rows.append([f"`{t}`", fmt(m["rows"]), fmt(m["row_groups"]),
                     human_bytes(m["bytes"]), TABLE_DESC.get(t, "")])
    if s.fps.is_file():
        rows.append(["`fpsim2_fingerprints.h5`", "—", "—",
                     human_bytes(s.fps.stat().st_size),
                     TABLE_DESC["fpsim2_fingerprints.h5"]])
    r.table(["表", "行数", "row group", "文件大小", "说明"], rows,
            ["---", "---:", "---:", "---:", "---"])

    r.h(3, "2.1 各表字段（从 parquet schema 实读）")
    for t in Snap.TABLES:
        if t not in info:
            continue
        cols = ", ".join(f"`{n}`:{ty}" for n, ty in info[t]["cols"])
        r.p(f"- **`{t}`** — {cols}")
    r.p("")
    r.p("注意 `patents` 的 `cpc` / `ipcr` / `ipc` / `ecla` / `assignee` 是 "
        "**LIST 类型**，`compounds` 与 `patents` 的字符串列是 `large_string`（BYTE_ARRAY）。")


def sec_patents(r: Report, s: Snap, info: dict, args) -> None:
    if not s.has("patents"):
        return
    n = info["patents"]["rows"]
    r.h(2, "3. 专利层面")

    r.h(3, "3.1 专利局分布")
    rows = s.q(f"SELECT country, COUNT(*) FROM {s.p('patents')} GROUP BY 1 ORDER BY 2 DESC",
               "专利局分布")
    vmax = max((c for _, c in rows), default=1)
    r.table(["country", "专利数", "占比", "", "说明"],
            [[f"`{c}`", fmt(k), pct(k, n), bar(k, vmax),
              COUNTRY_NOTE.get(c, "**异常值，数量极少，来源不明**")] for c, k in rows],
            ["---", "---:", "---:", "---", "---"])
    odd = [c for c, k in rows if c not in COUNTRY_NOTE]
    if odd:
        r.p(f"⚠ 出现了不在官方 5 家专利局之列的 `country` 值：{'、'.join(odd)}。"
            "数量极少，做统计时应显式排除或单独核查。")
        r.p("")

    r.h(3, "3.2 公开年份分布")
    rows = s.q(f"SELECT year(publication_date) y, COUNT(*) FROM {s.p('patents')} "
               "WHERE publication_date IS NOT NULL GROUP BY 1 ORDER BY 1 DESC LIMIT 25",
               "年份分布")
    if rows:
        vmax = max(c for _, c in rows)
        r.table(["年份", "专利数", ""],
                [[y, fmt(c), bar(c, vmax)] for y, c in rows],
                ["---:", "---:", "---"])
        nulls = s.one(f"SELECT COUNT(*) FROM {s.p('patents')} WHERE publication_date IS NULL")
        r.p(f"`publication_date` 为空的有 **{fmt(nulls)}** 篇（{pct(nulls, n)}）。"
            "最新一年通常不完整（快照日期之后的还没收录）。")
        r.p("")

    r.h(3, "3.3 专利同族 family_id ⭐")
    row = s.q(f"""SELECT COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0),
                         COUNT(*) FILTER (WHERE family_id < 0),
                         COUNT(*) FILTER (WHERE family_id IS NULL),
                         MIN(family_id)
                  FROM {s.p('patents')}""", "同族统计")
    if row:
        fam, sentinel, fam_null, fmin = row[0]
        real_docs = n - sentinel - fam_null
        r.p(f"**{fmt(real_docs)} 篇有真实同族的专利文档只对应 {fmt(fam)} 个同族**"
            f"（平均每族 {real_docs / fam:.2f} 篇）。")
        r.p("")
        r.p("**这是去重的关键**：直接数专利篇数会把同一发明重复计入。"
            "做「有多少个 GKA 发明」这类统计，一律 `COUNT(DISTINCT family_id)`。")
        r.p("")
        if sentinel or fam_null:
            r.p("⚠ **但 `family_id` 有两类无效值，必须先排除，否则 `COUNT(DISTINCT)` 会被污染**：")
            r.p("")
            r.table(["值", "文档数", "含义"],
                    [[f"`{fmt(fmin)}`（哨兵）", fmt(sentinel),
                      "**不是同族编号**，是「未分配同族」的占位。"
                      "不排除的话这 " + fmt(sentinel) + " 篇会被当成同一个发明"],
                     ["`NULL`", fmt(fam_null), "字段缺失"]],
                    ["---", "---:", "---"])
            r.code("-- 正确写法\n"
                   "COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)", "sql")
        dn = s.one(f"SELECT COUNT(*) FROM {s.p('patents')} "
                   "WHERE family_id IS NULL AND publication_date IS NULL")
        if dn and dn == fam_null:
            r.p(f"另外实测：**`family_id` 为空的 {fmt(fam_null)} 篇，"
                "与 `publication_date` 为空的完全是同一批**。"
                "说明这是一整块元数据缺失的记录，不是随机缺失——"
                "做时间趋势或同族分析时它们会整体消失，要意识到这个盲区。")
            r.p("")
        rows = s.q(f"""SELECT family_id, COUNT(*) c FROM {s.p('patents')}
                       WHERE family_id > 0 GROUP BY 1 ORDER BY 2 DESC LIMIT 5""",
                   "最大同族")
        if rows:
            r.p("最大的几个真实同族（一项发明在全球公开了多少次）：")
            r.p("")
            r.table(["family_id", "文档数"], [[f, fmt(c)] for f, c in rows], ["---", "---:"])

    r.h(3, "3.4 申请人 assignee")
    r.p("`assignee` 是 LIST 类型，要 `UNNEST` 后再统计。")
    r.p("")
    if args.deep:
        rows = s.q(f"""SELECT a, COUNT(*) c FROM (
                         SELECT UNNEST(assignee) a FROM {s.p('patents')})
                       WHERE a IS NOT NULL AND a <> '' GROUP BY 1 ORDER BY 2 DESC LIMIT 15""",
                   "申请人 Top15")
        if rows:
            r.table(["申请人", "专利数"], [[a[:60], fmt(c)] for a, c in rows], ["---", "---:"])
        r.p("申请人名称**没有做机构消歧**——同一家公司会有多种写法（含子公司、"
            "不同语言转写、OCR 噪声）。要按公司统计必须自己归一。")
    else:
        r.p("（`--deep` 模式才统计，需 UNNEST 4,491 万行）")
    r.p("")


def sec_compounds(r: Report, s: Snap, info: dict, args) -> None:
    if not s.has("compounds"):
        return
    n = info["compounds"]["rows"]
    r.h(2, "4. 化合物层面")
    r.p(f"**{fmt(n)}** 个唯一化学结构。每行给 `smiles` / `inchi` / `inchi_key` / `mol_weight`，"
        "**没有活性、没有名称、没有理化性质**——要这些得回 ChEMBL 或自己算。")
    r.p("")

    r.h(3, "4.1 结构标识的完整性")
    rows = []
    for col in ("smiles", "inchi", "inchi_key"):
        miss = s.one(f"SELECT COUNT(*) FROM {s.p('compounds')} "
                     f"WHERE {col} IS NULL OR {col}=''", f"{col} 缺失")
        rows.append([f"`{col}`", fmt(n - miss), fmt(miss), pct(miss, n)])
    r.table(["字段", "有值", "缺失", "缺失率"], rows, ["---", "---:", "---:", "---:"])

    if args.deep:
        uniq_ik = s.one(f"SELECT COUNT(DISTINCT inchi_key) FROM {s.p('compounds')}",
                        "InChIKey 唯一性")
        if uniq_ik:
            r.p(f"`inchi_key` 唯一值 **{fmt(uniq_ik)}** / {fmt(n)}"
                + ("——**完全唯一，可直接当跨库主键**。" if uniq_ik == n
                   else f"——有 {fmt(n - uniq_ik)} 个重复，跨库对齐前要先去重。"))
            r.p("")

    r.h(3, "4.2 分子量分布")
    row = s.q(f"""SELECT COUNT(*), MIN(mol_weight), quantile_cont(mol_weight,0.25),
                         median(mol_weight), quantile_cont(mol_weight,0.75),
                         quantile_cont(mol_weight,0.99), MAX(mol_weight)
                  FROM {s.p('compounds')} WHERE mol_weight IS NOT NULL""", "分子量分位")
    if row:
        c, mn, q1, med, q3, p99, mx = row[0]
        r.table(["统计量", "值"],
                [["有值", fmt(c)], ["最小", fmt(mn)], ["25%", fmt(q1)],
                 ["中位", fmt(med)], ["75%", fmt(q3)], ["99%", fmt(p99)],
                 ["最大", fmt(mx)]], ["---", "---:"])
        zero = s.one(f"SELECT COUNT(*) FROM {s.p('compounds')} WHERE mol_weight = 0")
        big = s.one(f"SELECT COUNT(*) FROM {s.p('compounds')} WHERE mol_weight > 2000")
        r.p(f"⚠ **两端都有脏数据**：`mol_weight = 0` 的有 **{fmt(zero)}** 条，"
            f"`> 2000`（超出小分子范围，多为聚合物/多肽/OSR 误识别）有 **{fmt(big)}** 条，"
            f"最大值 {fmt(mx)}。**做小分子筛选务必加分子量下限与上限**。")
        r.p("")
    r.p("> 抽取是全自动的（含化学结构图 OSR 与 OCR），"
        "**必然存在错误结构**。SureChEMBL 的化合物没有经过人工审编，"
        "这与 ChEMBL 是本质区别——拿来做训练集或候选池前要自己过一遍 RDKit 合法性与合理性检查。")
    r.p("")


def sec_map(r: Report, s: Snap, info: dict, args) -> None:
    if not s.has("patent_compound_map"):
        return
    n = info["patent_compound_map"]["rows"]
    r.h(2, "5. 专利 ↔ 化合物关联（核心表）⭐")
    r.p(f"**{fmt(n)} 行**，全库最大。三列：`patent_id`、`compound_id`、`field_id`。")
    r.p("")

    r.h(3, "5.1 按出现部分 field_id 拆开 ⭐")
    fields = {fid: name for fid, name in s.q(f"SELECT id, field_name FROM {s.p('fields')}")}
    rows = s.q(f"SELECT field_id, COUNT(*) FROM {s.p('patent_compound_map')} "
               "GROUP BY 1 ORDER BY 2 DESC", "field 分布")
    vmax = max((c for _, c in rows), default=1)
    r.table(["field_id", "field_name", "关联数", "占比", ""],
            [[fid, f"`{fields.get(fid, '?')}`", fmt(c), pct(c, n), bar(c, vmax)]
             for fid, c in rows],
            ["---:", "---", "---:", "---:", "---"])
    desc_n = dict(rows).get(1, 0)
    clms_n = dict(rows).get(2, 0)
    if desc_n and clms_n:
        r.p(f"**说明书（`desc`）占了 {pct(desc_n, n)}，权利要求（`clms`）只有 "
            f"{pct(clms_n, n)}，相差 {desc_n / clms_n:.1f} 倍。**")
        r.p("")
        r.p("这就是本库最大的坑：**如果不按 `field_id` 过滤，你拿到的绝大部分是"
            "「说明书里被提到过」的化合物**——包括背景技术里引用的他人化合物、"
            "对比例、乃至试剂。判断「这篇专利要保护什么」必须 `field_id = 2`。")
        r.p("")

    r.h(3, "5.2 参与关联的实体数")
    if args.deep:
        row = s.q(f"""SELECT COUNT(DISTINCT patent_id), COUNT(DISTINCT compound_id)
                      FROM {s.p('patent_compound_map')}""", "distinct 统计")
        if row:
            np_, nc = row[0]
            r.table(["", "数量", "占该表总数"],
                    [["出现在关联表里的专利", fmt(np_), pct(np_, info["patents"]["rows"])],
                     ["出现在关联表里的化合物", fmt(nc), pct(nc, info["compounds"]["rows"])]],
                    ["---", "---:", "---:"])

        r.h(3, "5.3 每篇专利含多少化合物")
        row = s.q(f"""WITH c AS (SELECT patent_id, COUNT(DISTINCT compound_id) k
                                 FROM {s.p('patent_compound_map')} GROUP BY 1)
                      SELECT MIN(k), median(k), quantile_cont(k,0.9),
                             quantile_cont(k,0.99), MAX(k) FROM c""", "每专利化合物数（慢）")
        if row:
            mn, med, p90, p99, mx = row[0]
            r.table(["最小", "中位", "90%", "99%", "最大"],
                    [[fmt(mn), fmt(med), fmt(p90), fmt(p99), fmt(mx)]],
                    ["---:", "---:", "---:", "---:", "---:"])
            r.p(f"**分布极度右偏**：中位只有 {fmt(med)} 个，但最大一篇有 **{fmt(mx)}** 个化合物——"
                "那是马库什结构（Markush）被枚举展开的组合库专利。"
                "**这类专利会主导任何按化合物计数的统计**，做分析时要么按专利归一，"
                "要么把超大专利单列。")
            r.p("")
    else:
        r.p("（`--deep` 模式才统计：distinct 计数与分位数需要全表扫，约 1.5 分钟）")
        r.p("")


def sec_bio(r: Report, s: Snap, info: dict, args) -> None:
    if not s.has("biomedical_locations"):
        return
    r.h(2, "6. 生物医学标注")
    nloc = info["biomedical_locations"]["rows"]
    nent = info["biomedical_entities"]["rows"] if "biomedical_entities" in info else None
    r.p(f"`biomedical_entities` **{fmt(nent)}** 个实体，"
        f"`biomedical_locations` **{fmt(nloc)}** 条位置记录。")
    r.p("")

    r.h(3, "6.1 各类型实体的出现规模")
    rows = s.q(f"""SELECT t.type_name, COUNT(DISTINCT e.id) AS n_entities
                   FROM {s.p('biomedical_entities')} e
                   JOIN {s.p('biomedical_types')} t ON t.id = e.type_id
                   GROUP BY 1 ORDER BY 2 DESC""", "实体类型规模")
    r.table(["类型", "实体数", "说明"],
            [[f"`{t}`", fmt(c), TYPE_NOTE.get(t, "")] for t, c in rows],
            ["---", "---:", "---"])
    missing = [t for t in TYPE_NOTE if t not in {x[0] for x in rows}]
    if missing:
        r.p(f"⚠ 字典里定义了但**本快照一条实体都没有**的类型：{'、'.join(missing)}。"
            "字典有定义不等于有数据。")
        r.p("")

    r.h(3, "6.2 归一化率")
    rows = s.q(f"""SELECT t.type_name,
                          COUNT(*) AS total,
                          COUNT(*) FILTER (WHERE e.resolved_form IS NOT NULL
                                             AND e.resolved_form <> '') AS resolved
                   FROM {s.p('biomedical_entities')} e
                   JOIN {s.p('biomedical_types')} t ON t.id = e.type_id
                   GROUP BY 1 ORDER BY 2 DESC""", "归一化率")
    r.table(["类型", "实体数", "有 resolved_form", "归一化率"],
            [[f"`{t}`", fmt(tot), fmt(res), pct(res, tot)] for t, tot, res in rows],
            ["---", "---:", "---:", "---:"])
    r.p("**归一化率低意味着大量实体只能当自由文本用。** "
        "检索时按 `resolved_form` 锚定虽然会漏掉未归一的那部分，"
        "但按 `original_text` 匹配又会引入大量歧义——两头都有代价，"
        "正确做法是**先按 `resolved_form` 取，再人工看未归一实体里有没有该捞的**。")
    r.p("")

    r.h(3, "6.3 GCK / glucokinase 的实测（本项目锚点）⭐")
    ents = s.q(f"""SELECT id, original_text FROM {s.p('biomedical_entities')}
                   WHERE resolved_form = 'HGNC:4195' ORDER BY id""", "HGNC:4195 实体")
    if ents:
        r.p(f"`resolved_form = 'HGNC:4195'`（人 GCK 基因）对应 **{len(ents)} 个不同的 "
            "`original_text`**：")
        r.p("")
        r.code("、".join(f"{t}" for _, t in ents))
        r.p("里面有 OCR 破碎形（`glucokmase`、`gl uc ok i na se`）、"
            "ChEMBL 那边的名字（`Hexokinase 4`）、以及一批缩写。"
            "**按字符串搜 `glucokinase` 会漏掉其中大部分。**")
        r.p("")
        ids = ",".join(str(i) for i, _ in ents)
        rows = s.q(f"""SELECT f.field_name, COUNT(DISTINCT l.patent_id), SUM(l.count)
                       FROM {s.p('biomedical_locations')} l
                       JOIN {s.p('fields')} f ON f.id = l.field_id
                       WHERE l.entity_id IN ({ids})
                       GROUP BY 1 ORDER BY 2 DESC""", "GCK 按部分分布")
        if rows:
            r.table(["field", "专利数", "提及次数"],
                    [[f"`{f}`", fmt(a), fmt(b)] for f, a, b in rows],
                    ["---", "---:", "---:"])
        rows = s.q(f"""SELECT e.original_text, COUNT(DISTINCT l.patent_id) c
                       FROM {s.p('biomedical_locations')} l
                       JOIN {s.p('biomedical_entities')} e ON e.id = l.entity_id
                       WHERE e.resolved_form = 'HGNC:4195'
                       GROUP BY 1 ORDER BY 2 DESC LIMIT 12""", "GCK surface form")
        if rows:
            r.p("各写法贡献的专利数：")
            r.p("")
            r.table(["original_text", "专利数", "备注"],
                    [[f"`{t}`", fmt(c),
                      "⚠ 糖尿病文献里 `GK` 更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型），"
                      "与本领域高度混淆" if t == "GK" else
                      ("⚠ 单个数字被解析成基因，**标注错误**" if t.strip() == "4" else
                       ("细菌 glucokinase 基因名，不是人 GCK"
                        if t in ("glk", "GlkA", "gukA", "GlcK") else ""))]
                     for t, c in rows],
                    ["---", "---:", "---"])
        r.p("⚠ 另外这些**不**解析到 HGNC:4195，别混进来："
            "`Glucokinase regulator` → `HGNC:4196`（GCKR/GKRP，不是 GCK）、"
            "全大写 `GCK` → **空（未解析）**、`HPK/GCK-like kinase` → `HGNC:6866`（MAP4K 家族）。")
        r.p("")


def sec_integrity(r: Report, s: Snap, info: dict, args) -> None:
    r.h(2, "7. 数据完整性实测")
    if not args.deep:
        r.p("（`--deep` 模式才做，需要几次全表扫）")
        r.p("")
        return
    rows = []
    orph_p = s.one(f"""SELECT COUNT(*) FROM (SELECT DISTINCT patent_id
                       FROM {s.p('patent_compound_map')}
                       EXCEPT SELECT id FROM {s.p('patents')})""", "孤儿 patent_id")
    orph_c = s.one(f"""SELECT COUNT(*) FROM (SELECT DISTINCT compound_id
                       FROM {s.p('patent_compound_map')}
                       EXCEPT SELECT id FROM {s.p('compounds')})""", "孤儿 compound_id")
    no_cmp = s.one(f"""SELECT COUNT(*) FROM (SELECT id FROM {s.p('patents')}
                       EXCEPT SELECT DISTINCT patent_id
                       FROM {s.p('patent_compound_map')})""", "无化合物的专利")
    rows.append(["map 里引用了但 `patents` 中不存在的 `patent_id`", fmt(orph_p),
                 "⚠ 外键悬空" if orph_p else "✅"])
    rows.append(["map 里引用了但 `compounds` 中不存在的 `compound_id`", fmt(orph_c),
                 "⚠ 外键悬空" if orph_c else "✅"])
    rows.append(["`patents` 里一个化合物都没有的专利", fmt(no_cmp),
                 "见下方说明" if no_cmp == 0 else ""])
    r.table(["检查项", "数量", "结论"], rows, ["---", "---:", "---"])
    if no_cmp == 0:
        r.p("**「一个化合物都没有的专利 = 0」是个重要事实**：说明 bulk 导出**只收录了"
            "含化学结构的专利**，不是 SureChEMBL 系统里的全部文档。"
            "这解释了为什么本地实测的专利数（约 4,491 万）远小于官方宣传的 1.166 亿——"
            "两者口径不同。**做规模陈述以本地实测为准。**")
        r.p("")


def sec_recipes(r: Report, s: Snap) -> None:
    r.h(2, "8. 上手：可直接运行的查询")
    r.p("环境：micromamba `GKA_in_Brain`（`duckdb` + `pyarrow`）。"
        "**所有查询都直接打 parquet，不需要先导入数据库。**")
    r.p("")

    r.h(3, "8.1 基本设置")
    r.code("""import duckdb
con = duckdb.connect()
con.execute("SET threads=8; SET memory_limit='12GB';")
con.execute("SET enable_progress_bar=false;")   # 不关会把 stdout 刷爆

D = "/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17"
""", "python")

    r.h(3, "8.2 找出在权利要求里主张 glucokinase 的专利")
    r.p("注意三点：按 `resolved_form` 锚定、按 `field_id=2` 限定权利要求、按 `family_id` 去重。")
    r.code(f"""con.sql(f'''
WITH e AS (                                   -- 1. 按归一 ID 锚定，不用字符串匹配
  SELECT id FROM '{{D}}/biomedical_entities.parquet'
  WHERE resolved_form = 'HGNC:4195'
),
p AS (                                        -- 2. 只要权利要求里提到的
  SELECT DISTINCT l.patent_id
  FROM '{{D}}/biomedical_locations.parquet' l
  JOIN e ON e.id = l.entity_id
  WHERE l.field_id = 2
)
SELECT COUNT(*) AS n_docs,                    -- 3. 同族去重后才是「发明数」
       COUNT(DISTINCT pt.family_id) AS n_families
FROM p JOIN '{{D}}/patents.parquet' pt ON pt.id = p.patent_id
''').show()""", "python")

    r.h(3, "8.3 取某篇专利在权利要求里保护的化合物结构")
    r.code(f"""con.sql(f'''
SELECT c.id, c.smiles, c.inchi_key, c.mol_weight
FROM '{{D}}/patent_compound_map.parquet' m
JOIN '{{D}}/compounds.parquet' c ON c.id = m.compound_id
WHERE m.patent_id = 12345678
  AND m.field_id = 2                          -- 只要权利要求
''').show()""", "python")

    r.h(3, "8.4 反查：某个结构出现在哪些专利里")
    r.p("先用 InChIKey 定位 compound_id，再查关联表——**别拿 SMILES 做字符串匹配**，"
        "同一结构有多种 SMILES 写法。")
    r.code(f"""con.sql(f'''
WITH c AS (
  SELECT id FROM '{{D}}/compounds.parquet'
  WHERE inchi_key = 'KJSGTWFWVTYPFZ-AWEZNQCLSA-N'   -- MK-0941
)
SELECT pt.patent_number, pt.country, pt.publication_date, m.field_id
FROM '{{D}}/patent_compound_map.parquet' m
JOIN c ON c.id = m.compound_id
JOIN '{{D}}/patents.parquet' pt ON pt.id = m.patent_id
ORDER BY pt.publication_date
''').show()""", "python")

    r.h(3, "8.5 按分类号过滤（LIST 类型，不能用 LIKE）")
    r.code(f"""con.sql(f'''
SELECT COUNT(*) FROM '{{D}}/patents.parquet'
WHERE list_contains(cpc, 'A61P3/10')          -- 抗糖尿病用途
''').show()""", "python")

    r.h(3, "8.6 相似性检索用 FPSim2，不要自己算指纹")
    r.code("""# 官方已提供 3,099 万化合物的预计算指纹（1.26 GB）
from FPSim2 import FPSim2Engine
fpe = FPSim2Engine("fpsim2_fingerprints.h5")
results = fpe.similarity("CCS(=O)(=O)c1ccc(...)cn1", 0.7, n_workers=4)
# 返回 (compound_id, similarity)，再拿 compound_id 回 patent_compound_map""", "python")
    r.p("> FPSim2 需另装（`pip install FPSim2`），本项目环境暂未安装。")
    r.p("")


def sec_pitfalls(r: Report, s: Snap) -> None:
    r.h(2, "9. 用这份数据前必须知道的坑")
    r.p("按踩坑代价从大到小排。")
    r.p("")
    items = [
        ("**`field_id` 不过滤，结论就是错的**",
         "说明书关联数是权利要求的 6.5 倍。「专利里出现过某化合物」和「这篇专利要保护它」"
         "完全是两回事——说明书含背景技术，会大段引用他人化合物。"
         "**判断专利主张范围只能用 `field_id = 2`。**"),
        ("**不按 `family_id` 去重，数字会虚高；去重时不排除哨兵值，又会算少**",
         "同一项发明在多国多阶段公开，每次都是一篇独立文档，"
         "「有多少个 GKA 发明」必须 `COUNT(DISTINCT family_id)`。"
         "但 `family_id = -1` 是「未分配同族」的**哨兵值**（7 万余篇），"
         "直接 `COUNT(DISTINCT)` 会把它们错当成同一个发明。"
         "正确写法：`COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)`。"),
        ("**按 `original_text` 做文本匹配一定会漏**",
         "专利全文来自扫描件 OCR，同一概念有几十种字面写法。"
         "必须按 `resolved_form` 锚定。但反过来，缩写形（如 `GK`）会引入大量歧义，"
         "**锚定后仍需逐个 surface form 评估**。"),
        ("**这里没有活性数据**",
         "SureChEMBL 只说「化合物出现在专利里」，不说它有没有效、多强。"
         "活性要回 ChEMBL 或读专利原文。**不能把「出现在 GCK 专利里」当成「是 GKA」。**"),
        ("**结构是全自动抽取的，没有人工审编**",
         "含化学结构图识别（OSR）与 OCR，必然有错误结构。"
         "`mol_weight = 0` 与 `> 2000` 的都存在。"
         "拿来做候选池前要自己过 RDKit 合法性检查并设分子量窗口。"),
        ("**马库什结构会主导计数**",
         "组合库专利可以枚举出上万个化合物，"
         "任何「按化合物计数」的统计都会被这类专利带偏。要么按专利归一，要么单列。"),
        ("**`resolved_form` 不是单一命名空间**",
         "`HGNC:4195` 带前缀，`Q14397` 是裸 UniProt accession。解析前先判断命名空间。"),
        ("**`cpc` / `ipc` / `assignee` 是 LIST 类型**",
         "`LIKE` 匹配不上，要用 `list_contains()` / `UNNEST()`。"
         "申请人名称未做机构消歧，同一公司多种写法。"),
        ("**`WO` 不是国家**",
         "是 PCT 国际申请，与各国文档并存，不能和 `US`/`CN` 并列做「国别分布」解读。"),
        ("**规模数字以本地实测为准**",
         "官方宣传 1.166 亿专利 / 4,770 万化合物，"
         "本地 bulk 快照实测是 4,491 万 / 3,099 万——bulk 只含有化学标注的专利，口径不同。"),
        ("**快照两周一覆盖**",
         "`latest/` 会被覆盖，必须固定到日期目录，否则结论无法复现。"),
    ]
    for i, (t, d) in enumerate(items, 1):
        r.p(f"{i}. {t}")
        r.p(f"   {d}")
        r.p("")

    r.h(2, "10. 与 ChEMBL 侧的衔接")
    r.p("本项目的 ChEMBL 侧锚点是 **UniProt `P35557`**，SureChEMBL 侧是 **`HGNC:4195`**，"
        "两者是同一个基因（GCK）的不同命名空间。**跨库对齐走基因层，不要字符串比对蛋白名。**")
    r.p("")
    r.p("化合物层面的对齐用 **InChIKey**——两边都有，且都是标准 InChI 生成的。"
        "注意 InChIKey 前 14 位相同只代表骨架相同，"
        "**立体异构体会被合并**，而手性通常决定 GKA 的活性。")
    r.p("")
    r.table(["", "ChEMBL", "SureChEMBL"],
            [["靶点锚点", "`P35557`（UniProt）", "`HGNC:4195`"],
             ["化合物对齐键", "`standard_inchi_key`", "`inchi_key`"],
             ["能回答", "这个化合物**有多强**", "这个化合物**在谁的专利里**"],
             ["数据质量", "人工审编", "全自动抽取，需自己过滤"]])


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="调查 SureChEMBL 2.0 bulk parquet 快照，生成 Markdown 概览报告。")
    ap.add_argument("snapshot", type=Path, help="快照目录，如 SureChEMBL_2026-07-17/")
    ap.add_argument("-o", "--out", type=Path, default=None)
    ap.add_argument("--quick", action="store_true", help="只用 parquet 元数据，秒级")
    ap.add_argument("--deep", action="store_true", help="加做分位数与完整性全扫（约 2-3 分钟）")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="12GB")
    args = ap.parse_args()

    if not args.snapshot.is_dir():
        print(f"错误：{args.snapshot} 不是目录", file=sys.stderr)
        return 1

    s = Snap(args.snapshot, args.threads, args.memory)
    if not s.files:
        print(f"错误：{args.snapshot} 下没找到任何 parquet 文件", file=sys.stderr)
        return 1
    log(f"快照 {s.root}，找到 {len(s.files)} 张表"
        + ("（quick 模式）" if args.quick else ("（deep 模式）" if args.deep else "")))

    r = Report()
    info = sec_intro(r, s, args)
    sec_glossary(r, s)
    sec_schema(r, s, info)
    if not args.quick:
        sec_patents(r, s, info, args)
        sec_compounds(r, s, info, args)
        sec_map(r, s, info, args)
        sec_bio(r, s, info, args)
        sec_integrity(r, s, info, args)
    sec_recipes(r, s)
    sec_pitfalls(r, s)

    r.p("")
    r.p("---")
    r.p("")
    r.p(f"> 本报告由 `surechembl_profile.py` 自动生成，"
        f"耗时 {time.time() - T0:.1f} 秒"
        + ("（quick）" if args.quick else ("（deep）" if args.deep else "（默认）")) + "。")

    out = args.out or (s.root.parent / f"surechembl_{s.root.name.split('_')[-1]}_profile_report.md")
    r.write(out)
    log(f"报告已写出：{out}")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
