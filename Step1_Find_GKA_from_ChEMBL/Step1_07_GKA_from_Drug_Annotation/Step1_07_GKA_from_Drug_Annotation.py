#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_07_GKA_from_Drug_Annotation
=================================

补 Step1 的召回缺口：把**因为没有活性数据而被整条链路漏掉**的 GKA 找回来。

为什么会漏
----------
Step1_01–06 的链路每一步都挂在 `activities` 上：

    靶点 → assays → 激活 assay → activities → 分子

**零活性记录的分子在结构上不可见。** ChEMBL 的分子有两条独立入口——
从文献活性数据抽取（Step1 走的这条），和从药名/临床登记册收录
（USAN / INN / ATC / Clinical Candidates，完全没走）。

多格列艾汀（`CHEMBL4297508`，III 期，中国已上市）就是典型：`compound_records`
里 4 条来源全是 `doc_type = DATASET`，没有一篇论文，`activities = 0`。

从哪几张表找
------------
路径 A  `drug_mechanism` × `target_dictionary`
        挂在 CHEMBL3820 上的记录，直接给 `action_type`（实测 6 条全是 ACTIVATOR）。
        **这是人工审编的方向标注**——Step1_03/05 从 assay 描述辛苦推的结论，这里现成。
        安全网：另外不限 tid 扫 `mechanism_of_action` 文本，防止 GKA 被挂到别的靶点。

路径 B  `usan_stems` × `molecule_dictionary.usan_stem`
        `-gliatin` 的官方注释就是 "glucokinase activator"，按药名词干识别。

路径 C  `molecule_synonyms`
        兜住 `usan_stem` 字段没填、但同义词里有 -gliatin 名的
        （如 `CHEMBL4297399` 主名是 LY-2608204，同义词才是 Globalagliatin）。

⚠ 同一个药有多个 ChEMBL ID（MK-0941 游离碱 vs 药物条目；LY-2608204 vs 其盐酸盐），
  用 `molecule_hierarchy` 给出 `parent_chembl_id` 供归并，产物不擅自合并行。

口径
----
本步骤的分子**没有活性数值**，无法参与 Step1_05 的效力分档与排序，
因此**单独成表**、带 `source = drug_annotation`，不并入 782 个候选。

输出
----
    Step1_07_GKA_from_Drug_Annotation.csv
    Step1_07_GKA_from_Drug_Annotation.md

用法
----
    python3 Step1_07_GKA_from_Drug_Annotation.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_DB = REPO / "ChEMBL" / "ChEMBL_37" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"
DEFAULT_STEP1_06 = (HERE.parent / "Step1_06_GKA_Physicochemical_Property_Extraction"
                    / "Step1_06_GKA_Physicochemical_Properties.csv")

GCK_TARGET = "CHEMBL3820"          # Hexokinase-4，Step1_01 锚定的主靶点
USAN_STEM = "-gliatin"             # usan_stems.annotation = 'glucokinase activator'
SYNONYM_PATTERN = "%gliatin%"

# 视为「激活方向」的 action_type。ChEMBL 的 action_type 字典里激活类就这几个，
# 全列出来而不是只写 ACTIVATOR，换版本时不至于漏。
ACTIVATOR_ACTIONS = {"ACTIVATOR", "AGONIST", "PARTIAL AGONIST",
                     "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "OPENER"}

OUT_COLUMNS = [
    # --- 身份 ---
    "molecule_chembl_id", "molregno", "molecule_pref_name", "max_phase",
    "first_approval", "molecule_type", "structure_type",
    # --- 命中路径与证据（本步骤的核心）---
    "found_by", "n_paths",
    "dm_target_chembl_id", "dm_target_pref_name", "dm_mechanism_of_action",
    "dm_action_type", "dm_direct_interaction", "dm_molecular_mechanism",
    "usan_stem", "usan_stem_definition", "usan_year",
    "matched_synonyms",
    # --- 为什么以前漏了 ---
    "n_activities", "n_activities_on_gck", "record_sources", "has_literature_record",
    "in_step1_candidates",
    # --- 盐型归并 ---
    "parent_chembl_id", "is_parent",
    # --- 结构与理化（字段与 Step1_06 对齐）---
    "canonical_smiles", "standard_inchi_key",
    "mw_freebase", "full_mwt", "alogp", "hba", "hbd", "psa", "rtb",
    "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
    "full_molformula",
    # --- 其他名称 ---
    "all_synonyms",
]


def connect_readonly(db: Path) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.text_factory = lambda b: b.decode("utf-8", "replace")
    return con


def chembl_version(con) -> str:
    r = con.execute("SELECT chembl_release, creation_date FROM chembl_release "
                    "ORDER BY chembl_release_id DESC LIMIT 1").fetchone()
    return f"{r[0]}（{str(r[1])[:10]}）" if r else "未知"


# ---------------------------------------------------------------------------
# 三条锚定路径
# ---------------------------------------------------------------------------

def path_a_drug_mechanism(con) -> tuple:
    """drug_mechanism × target_dictionary。返回 (命中 dict, 安全网命中 dict)。"""
    hits = {}
    q = """
    SELECT md.molregno, md.chembl_id,
           t.chembl_id AS tgt_chembl_id, t.pref_name AS tgt_pref_name,
           dm.mechanism_of_action, dm.action_type,
           dm.direct_interaction, dm.molecular_mechanism
    FROM drug_mechanism dm
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN target_dictionary  t  ON t.tid = dm.tid
    WHERE t.chembl_id = ?
    """
    for r in con.execute(q, (GCK_TARGET,)):
        if (r["action_type"] or "").upper() not in ACTIVATOR_ACTIONS:
            continue
        hits[r["chembl_id"]] = dict(r)

    # 安全网：不限 tid，按机制文本再扫一遍
    net = {}
    q2 = """
    SELECT md.chembl_id, t.chembl_id AS tgt_chembl_id, t.pref_name AS tgt_pref_name,
           dm.mechanism_of_action, dm.action_type
    FROM drug_mechanism dm
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    LEFT JOIN target_dictionary t ON t.tid = dm.tid
    WHERE lower(dm.mechanism_of_action) LIKE '%glucokinase%'
       OR lower(dm.mechanism_of_action) LIKE '%hexokinase type iv%'
    """
    for r in con.execute(q2):
        net[r["chembl_id"]] = dict(r)
    return hits, net


def path_b_usan_stem(con) -> tuple:
    """usan_stems × molecule_dictionary.usan_stem。"""
    stem = con.execute("SELECT stem, annotation FROM usan_stems WHERE stem = ?",
                       (USAN_STEM,)).fetchone()
    hits = {}
    for r in con.execute(
            "SELECT molregno, chembl_id, usan_stem, usan_stem_definition, usan_year "
            "FROM molecule_dictionary WHERE usan_stem = ?", (USAN_STEM,)):
        hits[r["chembl_id"]] = dict(r)
    return hits, (dict(stem) if stem else None)


def path_c_synonyms(con) -> dict:
    hits = defaultdict(list)
    for r in con.execute(
            "SELECT md.chembl_id, ms.synonyms, ms.syn_type "
            "FROM molecule_synonyms ms JOIN molecule_dictionary md ON md.molregno = ms.molregno "
            "WHERE lower(ms.synonyms) LIKE ? ORDER BY md.chembl_id, ms.synonyms",
            (SYNONYM_PATTERN,)):
        s = {"name": r["synonyms"], "type": r["syn_type"]}
        if s not in hits[r["chembl_id"]]:
            hits[r["chembl_id"]].append(s)
    return dict(hits)


# ---------------------------------------------------------------------------

def enrich(con, chembl_ids: list) -> dict:
    """把身份、结构、理化、来源、活性计数、盐型一次取齐。"""
    ph = ",".join("?" * len(chembl_ids))
    out = {}
    q = f"""
    SELECT md.chembl_id, md.molregno, md.pref_name, md.max_phase, md.first_approval,
           md.molecule_type, md.structure_type, md.usan_stem, md.usan_stem_definition,
           md.usan_year,
           mh.parent_molregno, pmd.chembl_id AS parent_chembl_id,
           cs.canonical_smiles, cs.standard_inchi_key,
           cp.mw_freebase, cp.full_mwt, cp.alogp, cp.hba, cp.hbd, cp.psa, cp.rtb,
           cp.num_ro5_violations, cp.aromatic_rings, cp.heavy_atoms, cp.qed_weighted,
           cp.full_molformula
    FROM molecule_dictionary md
    LEFT JOIN molecule_hierarchy  mh  ON mh.molregno = md.molregno
    LEFT JOIN molecule_dictionary pmd ON pmd.molregno = mh.parent_molregno
    LEFT JOIN compound_structures cs  ON cs.molregno = md.molregno
    LEFT JOIN compound_properties cp  ON cp.molregno = md.molregno
    WHERE md.chembl_id IN ({ph})
    """
    for r in con.execute(q, chembl_ids):
        d = dict(r)
        d["is_parent"] = ("TRUE" if (d["parent_molregno"] is None
                                     or d["parent_molregno"] == d["molregno"]) else "FALSE")
        out[d["chembl_id"]] = d

    # 活性计数：总数 + 落在 GCK 靶点上的
    for r in con.execute(f"""
        SELECT md.chembl_id, COUNT(*) AS n FROM activities a
        JOIN molecule_dictionary md ON md.molregno = a.molregno
        WHERE md.chembl_id IN ({ph}) GROUP BY 1""", chembl_ids):
        out[r["chembl_id"]]["n_activities"] = r["n"]
    for r in con.execute(f"""
        SELECT md.chembl_id, COUNT(*) AS n FROM activities a
        JOIN molecule_dictionary md ON md.molregno = a.molregno
        JOIN assays s ON s.assay_id = a.assay_id
        JOIN target_dictionary t ON t.tid = s.tid
        WHERE md.chembl_id IN ({ph}) AND t.chembl_id = ? GROUP BY 1""",
                         chembl_ids + [GCK_TARGET]):
        out[r["chembl_id"]]["n_activities_on_gck"] = r["n"]

    # compound_records 的来源 —— 这是「为什么以前漏了」的直接证据
    srcs = defaultdict(list)
    for r in con.execute(f"""
        SELECT md.chembl_id, s.src_short_name, d.doc_type
        FROM compound_records cr
        JOIN molecule_dictionary md ON md.molregno = cr.molregno
        JOIN source s ON s.src_id = cr.src_id
        LEFT JOIN docs d ON d.doc_id = cr.doc_id
        WHERE md.chembl_id IN ({ph})""", chembl_ids):
        item = {"src": r["src_short_name"], "doc_type": r["doc_type"]}
        if item not in srcs[r["chembl_id"]]:
            srcs[r["chembl_id"]].append(item)
    for cid, v in srcs.items():
        out[cid]["record_sources"] = v
        out[cid]["has_literature_record"] = (
            "TRUE" if any(x["doc_type"] not in ("DATASET", None) for x in v) else "FALSE")

    syns = defaultdict(list)
    for r in con.execute(f"""
        SELECT md.chembl_id, ms.synonyms FROM molecule_synonyms ms
        JOIN molecule_dictionary md ON md.molregno = ms.molregno
        WHERE md.chembl_id IN ({ph}) ORDER BY ms.synonyms""", chembl_ids):
        if r["synonyms"] not in syns[r["chembl_id"]]:
            syns[r["chembl_id"]].append(r["synonyms"])
    for cid, v in syns.items():
        out[cid]["all_synonyms"] = v
    return out


def load_step1_candidates(path: Path) -> set:
    if not path.is_file():
        return set()
    with path.open(encoding="utf-8") as f:
        return {r["molecule_chembl_id"] for r in csv.DictReader(f)}


def write_csv(rows: list, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            o = {}
            for k in OUT_COLUMNS:
                v = r.get(k)
                if isinstance(v, (list, dict)):
                    v = json.dumps(v, ensure_ascii=False) if v else ""
                o[k] = "" if v is None else v
            w.writerow(o)


def write_report(rows, path, db, version, step1_ids, stem_info, net_extra) -> None:
    L = []
    n = len(rows)
    new = [r for r in rows if r["in_step1_candidates"] == "FALSE"]
    noact = [r for r in rows if not r.get("n_activities")]

    L.append("# Step1_07 从药物注释补回被漏掉的 GKA")
    L.append("")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db}`")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 命中分子：**{n}** 个 ChEMBL ID，其中 **{len(new)}** 个不在 Step1_06 的 782 个候选里")
    L.append("")
    L.append("> Step1_01–06 的链路每一步都挂在 `activities` 上，"
             "**零活性记录的分子在结构上不可见**。本步骤从药物注释侧独立锚定，把它们找回来。")
    L.append("")

    L.append("## 一、从哪几张表找的")
    L.append("")
    L.append("| 路径 | 表 | 判据 | 命中 |")
    L.append("| --- | --- | --- | ---: |")
    a = sum(1 for r in rows if "drug_mechanism" in r["found_by"])
    b = sum(1 for r in rows if "usan_stem" in r["found_by"])
    c = sum(1 for r in rows if "synonym" in r["found_by"])
    L.append(f"| **A** | `drug_mechanism` × `target_dictionary` | "
             f"`tid` → `{GCK_TARGET}` 且 `action_type` 属激活类 | {a} |")
    L.append(f"| **B** | `usan_stems` × `molecule_dictionary` | "
             f"`usan_stem = '{USAN_STEM}'` | {b} |")
    L.append(f"| **C** | `molecule_synonyms` | `synonyms LIKE '{SYNONYM_PATTERN}'` | {c} |")
    L.append("")
    if stem_info:
        L.append(f"路径 B 的依据是 `usan_stems` 表原文：**`{stem_info['stem']}` → "
                 f"「{stem_info['annotation']}」**。WHO/USAN 的命名规则里这个后缀"
                 "就是给 GKA 用的，**按药名就能认出来**。")
        L.append("")
    L.append("路径 A 的 `action_type` 实测分布：" + "、".join(
        f"`{k}` × {v}" for k, v in Counter(
            r["dm_action_type"] for r in rows if r.get("dm_action_type")).most_common()) + "。")
    L.append("")
    L.append("**路径 A 给的是人工审编的方向标注**——`mechanism_of_action` 写着 "
             "「Hexokinase type IV activator」。这正是 Step1_03/05 用规则 + LLM "
             "从 assay 描述里推的结论，ChEMBL 已经标好了，比推的可靠。")
    L.append("")
    L.append(f"**安全网**：另外不限 `tid` 扫了一遍 `mechanism_of_action` 文本"
             f"（含 glucokinase / hexokinase type IV），命中 {len(net_extra)} 条，"
             + ("与按 `tid` 筛的结果完全一致，没有 GKA 被挂到别的靶点上。"
                if not net_extra else
                f"**其中 {len(net_extra)} 条不在按 tid 的结果里，需人工核**。"))
    L.append("")

    L.append("## 二、命中清单")
    L.append("")
    L.append("| 分子 | 名称 | phase | 命中路径 | GCK 活性数 | 总活性数 | 在 782 里 |")
    L.append("| --- | --- | ---: | --- | ---: | ---: | :---: |")
    for r in sorted(rows, key=lambda x: (-(x.get("max_phase") or 0),
                                         x["molecule_chembl_id"])):
        L.append(f"| `{r['molecule_chembl_id']}` | {r.get('molecule_pref_name') or '—'} | "
                 f"{r.get('max_phase') or '—'} | {r['found_by']} | "
                 f"{r.get('n_activities_on_gck') or 0} | {r.get('n_activities') or 0} | "
                 f"{'✅' if r['in_step1_candidates'] == 'TRUE' else '**❌ 新增**'} |")
    L.append("")

    L.append("## 三、为什么以前会漏")
    L.append("")
    L.append(f"**{len(noact)} / {n}** 个分子的 `activities` 为 0。看它们的 "
             "`compound_records` 来源就明白了：")
    L.append("")
    L.append("| 分子 | 来源 | 有文献记录 |")
    L.append("| --- | --- | :---: |")
    for r in sorted(rows, key=lambda x: (x.get("n_activities") or 0)):
        srcs = r.get("record_sources") or []
        s = "、".join(f"`{x['src']}`({x['doc_type']})" for x in srcs) or "—"
        L.append(f"| `{r['molecule_chembl_id']}` {r.get('molecule_pref_name') or ''} | "
                 f"{s[:90]} | {'✅' if r.get('has_literature_record') == 'TRUE' else '❌'} |")
    L.append("")
    L.append("`doc_type = DATASET` 的来源（USAN / INN / ATC / CANDIDATES）是**药名与临床登记册**，"
             "不带活性数值。分子经这条路进 ChEMBL，`activities` 就是 0，"
             "Step1 那条 `靶点 → assay → activity → 分子` 的链路自然看不见它。")
    L.append("")

    L.append("## 四、⚠ 同一个药有多个 ChEMBL ID")
    L.append("")
    groups = defaultdict(list)
    for r in rows:
        key = r.get("parent_chembl_id") or r["molecule_chembl_id"]
        groups[key].append(r)
    dup = {k: v for k, v in groups.items() if len(v) > 1}
    known = [("MK-0941", ["CHEMBL3580737", "CHEMBL4297302"]),
             ("Globalagliatin / LY-2608204", ["CHEMBL4297399", "CHEMBL5095182"])]
    L.append("`molecule_hierarchy` 的 `parent_molregno` 能归并一部分，但**不是全部**——"
             "同一个药以「游离碱」「盐」「药物条目」等多种身份注册时，"
             "彼此之间未必有 parent 关系。实测到的：")
    L.append("")
    L.append("| 药 | 条目 | 说明 |")
    L.append("| --- | --- | --- |")
    for name, ids in known:
        present = [i for i in ids if any(r["molecule_chembl_id"] == i for r in rows)
                   or i in step1_ids]
        L.append(f"| {name} | " + "、".join(f"`{i}`" for i in ids) +
                 f" | 本表命中 {len(present)} 个 |")
    L.append("")
    L.append("**做去重统计时必须按结构（InChIKey）或人工归并，不能只按 ChEMBL ID 计数。**")
    L.append("")

    L.append("## 五、口径：不并入 782 个候选")
    L.append("")
    L.append("本表的分子**没有活性数值**，无法参与 Step1_05 的效力分档与排序，因此单独成表。"
             "下游合并时必须知道两点：")
    L.append("")
    L.append("1. 782 个候选有 `pactivity`、可排序；**这批没有**，只有方向标注")
    L.append("2. 这批的方向是 ChEMBL **人工审编**的（`action_type = ACTIVATOR`），"
             "比 Step1_05 从读数推出来的更可靠")
    L.append("")

    L.append("## 六、更新后的阳性对照集")
    L.append("")
    L.append("Step2 自检用。原 8 个 → 现 **%d** 个。" % (len(rows) + 2))
    L.append("")
    L.append("| 分子 | 名称 | phase | InChIKey | 来源 |")
    L.append("| --- | --- | ---: | --- | --- |")
    for r in sorted(rows, key=lambda x: (-(x.get("max_phase") or 0),
                                         x["molecule_chembl_id"])):
        L.append(f"| `{r['molecule_chembl_id']}` | {r.get('molecule_pref_name') or '—'} | "
                 f"{r.get('max_phase') or '—'} | `{r.get('standard_inchi_key') or '—'}` | "
                 f"Step1_07 |")
    L.append("")
    L.append("另加 Step1_06 里已有的 `CHEMBL1096435` Ro-28-1675（参比化合物，非临床药）"
             "与 `CHEMBL5072532` BMS-820132。")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="从药物注释侧补回被漏掉的 GKA。")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--step1-06", type=Path, default=DEFAULT_STEP1_06)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    if not args.db.is_file():
        print(f"错误：找不到数据库 {args.db}", file=sys.stderr)
        return 1
    con = connect_readonly(args.db)
    version = chembl_version(con)
    print(f"数据库：{args.db}")
    print(f"ChEMBL 版本：{version}\n")

    a_hits, a_net = path_a_drug_mechanism(con)
    b_hits, stem_info = path_b_usan_stem(con)
    c_hits = path_c_synonyms(con)
    print(f"路径 A  drug_mechanism（tid={GCK_TARGET}，激活类）：{len(a_hits)}")
    print(f"        安全网（不限 tid 的机制文本）：      {len(a_net)}")
    print(f"路径 B  usan_stem = '{USAN_STEM}'：            {len(b_hits)}")
    print(f"路径 C  synonyms LIKE '{SYNONYM_PATTERN}'：      {len(c_hits)}")
    net_extra = {k: v for k, v in a_net.items() if k not in a_hits}
    if net_extra:
        print(f"  ⚠ 安全网多出 {len(net_extra)} 条不在 tid 结果里：{list(net_extra)}")

    ids = sorted(set(a_hits) | set(b_hits) | set(c_hits))
    if not ids:
        print("没有命中任何分子。", file=sys.stderr)
        return 2
    print(f"\n三路合并：{len(ids)} 个 ChEMBL ID")

    info = enrich(con, ids)
    step1_ids = load_step1_candidates(args.step1_06)

    rows = []
    for cid in ids:
        d = info.get(cid, {"chembl_id": cid})
        found = []
        if cid in a_hits:
            found.append("drug_mechanism")
            a = a_hits[cid]
            d.update({"dm_target_chembl_id": a["tgt_chembl_id"],
                      "dm_target_pref_name": a["tgt_pref_name"],
                      "dm_mechanism_of_action": a["mechanism_of_action"],
                      "dm_action_type": a["action_type"],
                      "dm_direct_interaction": a["direct_interaction"],
                      "dm_molecular_mechanism": a["molecular_mechanism"]})
        if cid in b_hits:
            found.append("usan_stem")
        if cid in c_hits:
            found.append("synonym")
            d["matched_synonyms"] = c_hits[cid]
        d["molecule_chembl_id"] = cid
        d["molecule_pref_name"] = d.pop("pref_name", None)
        d["found_by"] = "+".join(found)
        d["n_paths"] = len(found)
        d["in_step1_candidates"] = "TRUE" if cid in step1_ids else "FALSE"
        d.setdefault("n_activities", 0)
        d.setdefault("n_activities_on_gck", 0)
        rows.append(d)

    new = [r for r in rows if r["in_step1_candidates"] == "FALSE"]
    noact = [r for r in rows if not r.get("n_activities")]
    print(f"  其中不在 Step1_06 的 782 个候选里：{len(new)}")
    print(f"  其中 activities = 0：             {len(noact)}")
    for r in sorted(new, key=lambda x: -(x.get("max_phase") or 0)):
        print(f"    {r['molecule_chembl_id']:<16s} {str(r.get('molecule_pref_name')):<30s} "
              f"phase={r.get('max_phase')}  活性={r.get('n_activities')}  [{r['found_by']}]")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = args.outdir / "Step1_07_GKA_from_Drug_Annotation.csv"
    out_md = args.outdir / "Step1_07_GKA_from_Drug_Annotation.md"
    write_csv(rows, out_csv)
    write_report(rows, out_md, args.db, version, step1_ids, stem_info, net_extra)
    print(f"\n主表：{out_csv}")
    print(f"报告：{out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
