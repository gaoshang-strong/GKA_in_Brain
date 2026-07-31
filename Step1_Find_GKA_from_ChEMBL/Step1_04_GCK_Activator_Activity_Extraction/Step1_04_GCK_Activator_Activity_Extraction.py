#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_04_GCK_Activator_Activity_Extraction
==========================================

从 Step1_03 已确认的「GCK 激活」assay 中提取全部 activity，关联到小分子，
并对每个分子分别给出三项结果：

    效力（potency）  EC50/AC50 类浓度指标，越小越强；同时给出 pActivity
    效能（efficacy）  % activation / fold activation，越大越强
    证据（evidence）  重复 assay 数、独立文献数、结果一致性

**三项分开记录，不合成总分；本步骤也不对 activity 排序**（README 要求）。
输出按 `molecule_chembl_id` 字典序排列，这是稳定的身份序，不含强弱含义。

输入
----
    ../Step1_03_GCK_Assay_Classification/Step1_03_Assay_Classification_final.csv

筛选条件（两条都必须满足）：

    final_category           == "GCK 激活"
    target_identity_suspect  == "FALSE"      ← 排除误挂在 CHEMBL3820 下的 MAP4K2

第二条是硬前置条件，来源见 `Step1_03_Target_Mismapping_MAP4K2.md`。

指标分流：为什么不是把 standard_type 直接当方向用
--------------------------------------------------
这批数据里的指标语义很不齐，直接按 `standard_type` 归并会算错。实测到的坑：

1. `S0.5` / `S50`（mM）**不是化合物效力**，而是酶对葡萄糖的半饱和常数。
   GKA 的作用就是把它从野生型的 ~7 mM 拉低到 0.5-2 mM，所以它是**机制读数**。
   若当成效力，会得到「EC50 = 0.6 mM」这种量级完全错误的结论。
   同理 `Km` / `Vmax` / `Kcat` 也是酶动力学参数，单列不进效力/效能。
2. `standard_units` 不可信。`Emax` 标着 `%`，值却是 1.04、1.02
   （相对参比化合物的倍数）；同一 `type + units` 组合内值域跨越 0.63–1454，
   混着百分数和倍数两种尺度。因此尺度必须**按 (assay, type) 分组按值域判定**，
   不能只看单位。
3. `Activity` / `Ratio` 是自由文本型指标（`activity_stds_lookup` 里没有标准化
   规则），含义随实验而异，**默认不跨实验汇总**。但其中确有真信息——
   `Activity`(uM) 的描述写着 "concentration required to 50% increase in enzyme
   activity"（即 AC50），`Ratio` 的描述写着 "ratio of enzyme activation in
   treated to untreated"（即 fold activation）。
   处理方式沿用 Step1_03 的思路：**带证据的受限提升**——只有 assay 描述命中显式
   模式、且单位与目标尺度自洽时才归类，命中的规则原文写进输出供复核；
   没命中的留在 `unclassified_metrics` 里，如实记录而不丢弃。

删失值
------
`EC50` 有 86 条 `standard_relation = '>'`（最高浓度仍未达到 EC50，实为未测出活性）。
这类**不参与效力汇总**，单独计数到 `n_ec50_censored_gt`——把它们当等号值平均，
会把弱化合物算强。

输出
----
    Step1_04_GCK_Activator_Activity_Extraction.csv   一行一个分子（主产物）
    Step1_04_GCK_Activator_Activities.csv            一行一条 activity（明细，可追溯）
    Step1_04_GCK_Activator_Activity_Extraction.md    报告

用法
----
    python3 Step1_04_GCK_Activator_Activity_Extraction.py
    python3 Step1_04_GCK_Activator_Activity_Extraction.py --db /path/to/chembl_37.db
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_DB = REPO / "ChEMBL" / "ChEMBL_37" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"
DEFAULT_IN = (HERE.parent / "Step1_03_GCK_Assay_Classification"
              / "Step1_03_Assay_Classification_final.csv")

ACTIVATION_CATEGORY = "GCK 激活"

# ---------------------------------------------------------------------------
# 指标分流表
# ---------------------------------------------------------------------------
# 效力：标准化的浓度型指标，值越小越强
POTENCY_STD = {"EC50", "AC50"}
# 效能：激活幅度，值越大越强（尺度另判，见 detect_scale）
EFFICACY_STD = {"Emax", "%max", "max activation"}
EFFICACY_FOLD_STD = {"FC"}
# 酶动力学：机制读数，既不是效力也不是效能，单列保留
KINETIC = {"S0.5", "S50", "Km", "Vmax", "Kcat", "kcat", "Ki"}
# 自由文本型：默认不汇总，只在描述命中显式模式时受限提升
FREE_TEXT = {"Activity", "Ratio", "Ratio EC50", "Percent Effect", "Z score"}

# 各 standard_type 的物理含义，只用于报告表的说明列，不参与任何计算。
# 写在这里而不是报告里，是为了重跑报告时说明不丢失。
STANDARD_TYPE_MEANING = {
    "EC50": "半最大效应浓度：激活幅度达到该化合物自身最大值一半时所需的浓度。**越小效力越强**",
    "AC50": "半最大激活浓度，与 `EC50` 同义的另一种写法。**越小效力越强**",
    "Activity": "ChEMBL 的自由文本指标，**本身没有固定含义**——本数据里既有浓度（nM/uM/mM）"
                "也有百分比和无量纲值。只在描述明写 half-maximal activation / 50% increase 时才提升为效力",
    "FC": "fold change：处理组酶活 ÷ 未处理对照，**基线 1**，>1 为激活",
    "Ratio": "比值，**分母不统一**：多数是未处理对照（基线 1），"
             "少数是参比激活剂如 Ro-28-1675（基线 0）；本数据里还混有 Km ratio、Vmax ratio。"
             "必须逐 assay 从描述读基线",
    "Ratio EC50": "两个条件下 EC50 的比值（如 ±4% 人血清白蛋白），衡量条件造成的效力位移，"
                  "**不是绝对效力**",
    "Emax": "最大效应：化合物浓度饱和后能达到的激活上限（回答「能激活到多强」，"
            "与 EC50 的「多低浓度就起效」是两回事）。本数据里有的以倍数记、有的以 % 记",
    "%max": "达到最大激活的百分数；**基线随 assay 而异**（0 或 100），不能设全局阈值",
    "max activation": "最大激活幅度，`Emax` 的自由写法",
    "S0.5": "酶活达到 Vmax 一半时的**葡萄糖**浓度。GCK 是正协同、非米氏动力学，故不写作 Km。"
            "GKA 把它从野生型 ~7 mM 拉到 0.5–2 mM——这是**酶的性质，不是化合物效力**",
    "S50": "同 `S0.5`，另一种写法。**酶的性质，不是化合物效力**",
    "Km": "米氏常数：酶对底物（此处葡萄糖）的半饱和浓度，GKA 使其下降。**酶的性质**",
    "Vmax": "最大反应速度：底物饱和时的催化上限。多数 GKA 只降 S0.5 不改 Vmax（K 型），"
            "改 Vmax 的是 V 型。**酶的性质**",
    "Kcat": "转换数：单个酶分子单位时间催化的底物分子数，催化效率上限。**酶的性质**",
    "kcat": "同 `Kcat`（大小写变体）。**酶的性质**",
    "Ki": "抑制常数，抑制剂与酶的解离常数。**方向与激活相反**",
    "Percent Effect": "自由文本的百分比效应，基线与满标随实验而异，不可跨实验汇总",
    "Z score": "筛选的统计量（偏离对照几个标准差），是**筛选质量指标**不是活性强度",
}

CONC_UNITS = {"nM", "uM", "mM", "M", "pM", "ug ml-1", "ng ml-1"}
# 换算到 nM
TO_NM = {"pM": 1e-3, "nM": 1.0, "uM": 1e3, "mM": 1e6, "M": 1e9}

# 受限提升的门槛：assay 描述必须命中，命中的原文写进输出
RE_PROMOTE_POTENCY = re.compile(
    r"concentration required to[^.;]{0,60}(?:increase|activat|stimulat)"
    r"|half[- ]maximal(?:ly)? (?:effective|activat)"
    r"|\bAC50\b|\bEC50\b", re.I)
RE_PROMOTE_FOLD = re.compile(
    r"ratio of[^.;]{0,60}(?:activation|activity)[^.;]{0,40}(?:treated|control)"
    r"|relative to (?:the )?(?:untreated |vehicle )?control"
    r"|fold[- ]activation|fold increase|\btimes\b[^.;]{0,20}control", re.I)
RE_PROMOTE_PCT = re.compile(
    r"percent(?:age)? (?:of )?activation|% ?activation"
    r"|increase in enzyme activity at|activation at [\d.]+ ?[num]M", re.I)

# 一致性判定阈值（pActivity 对数单位）
SPREAD_TIGHT = 0.5      # 3 倍以内
SPREAD_LOOSE = 1.0      # 10 倍以内

ACT_SQL = """
SELECT
    act.activity_id            AS activity_id,
    act.molregno               AS molregno,
    md.chembl_id               AS molecule_chembl_id,
    md.pref_name               AS molecule_pref_name,
    md.molecule_type           AS molecule_type,
    md.max_phase               AS max_phase,
    cs.canonical_smiles        AS canonical_smiles,
    cs.standard_inchi_key      AS standard_inchi_key,
    a.chembl_id                AS assay_chembl_id,
    a.description              AS assay_description,
    act.standard_type          AS standard_type,
    act.standard_relation      AS standard_relation,
    act.standard_value         AS standard_value,
    act.standard_units         AS standard_units,
    act.pchembl_value          AS pchembl_value,
    act.activity_comment       AS activity_comment,
    act.data_validity_comment  AS data_validity_comment,
    act.potential_duplicate    AS potential_duplicate,
    d.chembl_id                AS doc_chembl_id,
    d.year                     AS doc_year,
    s.src_short_name           AS src_short_name
FROM activities act
JOIN assays a                       ON a.assay_id = act.assay_id
JOIN molecule_dictionary md         ON md.molregno = act.molregno
LEFT JOIN compound_structures cs    ON cs.molregno = act.molregno
LEFT JOIN docs d                    ON d.doc_id = act.doc_id
LEFT JOIN source s                  ON s.src_id = act.src_id
WHERE a.chembl_id IN ({placeholders})
ORDER BY md.chembl_id, a.chembl_id, act.activity_id
"""

MOL_COLUMNS = [
    # --- 分子身份 ---
    "molecule_chembl_id", "molecule_pref_name", "molecule_type", "max_phase",
    "canonical_smiles", "standard_inchi_key",
    # --- 效力 ---
    "n_potency_records", "potency_types",
    "potency_nm_min", "potency_nm_median", "potency_nm_max",
    "pactivity_max", "pactivity_median", "pactivity_min", "pactivity_spread",
    "n_potency_censored_gt", "potency_censored_gt_min_nm", "n_potency_censored_lt",
    "potency_n_assays", "potency_n_docs", "potency_records",
    # --- 效能 ---
    "n_efficacy_records",
    "efficacy_fold_max", "efficacy_fold_max_source", "n_efficacy_fold",
    "efficacy_pct_max", "efficacy_pct_max_source", "n_efficacy_pct",
    "efficacy_n_assays", "efficacy_n_docs", "efficacy_records",
    # --- 酶动力学（机制读数，不计入效力/效能）---
    "n_kinetic_records", "kinetic_records",
    # --- 未归类（自由文本且描述未命中提升门槛）---
    "n_unclassified_records", "unclassified_records",
    # --- 证据 ---
    "n_activities", "n_assays", "n_docs", "doc_chembl_ids", "doc_years",
    "src_short_names", "n_potential_duplicate",
    "evidence_consistency", "evidence_level", "evidence_note",
]

ACT_COLUMNS = [
    "molecule_chembl_id", "molecule_pref_name", "assay_chembl_id",
    "activity_id", "standard_type", "standard_relation", "standard_value",
    "standard_units", "pchembl_value",
    "metric_role", "metric_scale", "value_nm", "pactivity", "role_rule",
    "glucose_mM_parsed", "activity_comment", "potential_duplicate",
    "doc_chembl_id", "doc_year", "src_short_name", "assay_description",
]

# 葡萄糖浓度只写在描述自由文本里（见 CLAUDE.md）。这里顺手解析出来放进明细表，
# 供后续按低糖/高糖分层用；**本步骤不据此做任何判定**。
RE_GLUCOSE = re.compile(r"(\d+(?:\.\d+)?)\s*mM\s+glucose", re.I)


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


def load_activation_assays(path: Path) -> list:
    """读 Step1_03 的终表，取激活且靶点身份无疑的 assay。"""
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    keep, dropped_suspect = [], 0
    for r in rows:
        if r.get("final_category") != ACTIVATION_CATEGORY:
            continue
        if r.get("target_identity_suspect") == "TRUE":
            dropped_suspect += 1
            continue
        keep.append(r)
    return keep, dropped_suspect


def to_nm(value, units: str):
    """浓度换算到 nM；非浓度单位返回 None。"""
    if value is None or units not in TO_NM:
        return None
    return float(value) * TO_NM[units]


def p_activity(value_nm):
    """pActivity = -log10(mol/L)。值越大越强。"""
    if not value_nm or value_nm <= 0:
        return None
    return round(9.0 - math.log10(value_nm), 2)


def detect_scale(group_max, units: str) -> str:
    """按 (assay, standard_type) 分组的值域判定效能尺度。

    单位不可信：`Emax` 标着 `%` 而值是 1.04（倍数）。以值域为准，
    单位只作为次要参考。
    """
    if group_max is None:
        return "unknown"
    if units == "%":
        # 百分数不会集中在 10 以下；这种是被错标单位的倍数
        return "fold" if group_max <= 10 else "percent"
    if units in (None, "", "-"):
        return "fold" if group_max <= 50 else "percent"
    return "other"


def assign_role(rec: dict, group_max) -> tuple:
    """给单条 activity 定角色。返回 (role, scale, rule)。

    role ∈ potency / efficacy / kinetic / unclassified
    rule 是判定依据；受限提升的必须带上命中的描述原文，供人工复核。
    """
    st = (rec["standard_type"] or "").strip()
    units = rec["standard_units"]
    desc = rec["assay_description"] or ""

    if st in KINETIC:
        return "kinetic", "", f"{st} 属酶动力学参数（机制读数），不计入效力/效能"

    if st in POTENCY_STD:
        if units in CONC_UNITS:
            return "potency", "concentration", f"{st} 为标准化浓度型效力指标"
        return "unclassified", "", f"{st} 缺少浓度单位（units={units!r}），无法作为效力"

    if st in EFFICACY_FOLD_STD:
        return "efficacy", "fold", f"{st}（fold change）为激活倍数"

    if st in EFFICACY_STD:
        scale = detect_scale(group_max, units)
        if scale in ("fold", "percent"):
            note = f"{st} 为激活幅度指标，按该 assay 内值域（max={group_max}）判为 {scale} 尺度"
            if units == "%" and scale == "fold":
                note += "；⚠ 原始单位标注为 % 但数值是倍数，以数值为准"
            return "efficacy", scale, note
        return "unclassified", "", f"{st} 尺度无法判定（units={units!r}, max={group_max}）"

    if st in FREE_TEXT:
        # 自由文本指标：默认不汇总，只有描述命中显式模式才提升
        m = RE_PROMOTE_POTENCY.search(desc)
        if m and units in CONC_UNITS:
            return "potency", "concentration", (
                "自由文本指标 " + st + " 受限提升为效力，描述命中：「" + m.group(0).strip() + "」")
        m = RE_PROMOTE_FOLD.search(desc)
        if m and detect_scale(group_max, units) == "fold":
            return "efficacy", "fold", (
                "自由文本指标 " + st + " 受限提升为效能(fold)，描述命中：「" + m.group(0).strip() + "」")
        m = RE_PROMOTE_PCT.search(desc)
        if m and detect_scale(group_max, units) == "percent":
            return "efficacy", "percent", (
                "自由文本指标 " + st + " 受限提升为效能(%)，描述命中：「" + m.group(0).strip() + "」")
        return "unclassified", "", (
            st + " 为自由文本指标（无官方标准化规则），assay 描述未命中提升门槛，不跨实验汇总")

    return "unclassified", "", f"未在指标分流表中登记的 standard_type：{st}"


def fetch_activities(con: sqlite3.Connection, assay_ids: list) -> list:
    ph = ",".join("?" * len(assay_ids))
    rows = [dict(r) for r in con.execute(ACT_SQL.format(placeholders=ph), assay_ids)]

    # 尺度判定要看同一 (assay, standard_type) 内的整体值域，先做一遍分组统计
    gmax: dict = {}
    for r in rows:
        if r["standard_value"] is None:
            continue
        k = (r["assay_chembl_id"], r["standard_type"])
        v = float(r["standard_value"])
        gmax[k] = v if k not in gmax else max(gmax[k], v)

    for r in rows:
        g = gmax.get((r["assay_chembl_id"], r["standard_type"]))
        role, scale, rule = assign_role(r, g)
        r["metric_role"] = role
        r["metric_scale"] = scale
        r["role_rule"] = rule
        r["value_nm"] = to_nm(r["standard_value"], r["standard_units"]) if role == "potency" else None
        if role == "potency" and r["standard_relation"] == "=":
            r["pactivity"] = (round(float(r["pchembl_value"]), 2)
                              if r["pchembl_value"] is not None
                              else p_activity(r["value_nm"]))
        else:
            r["pactivity"] = None
        m = RE_GLUCOSE.search(r["assay_description"] or "")
        r["glucose_mM_parsed"] = m.group(1) if m else ""
    return rows


def _stat(vals: list, fn):
    return round(fn(vals), 3) if vals else ""


def summarize_molecule(recs: list) -> dict:
    """把一个分子的全部 activity 汇总成效力 / 效能 / 证据三块。"""
    head = recs[0]
    out = {
        "molecule_chembl_id": head["molecule_chembl_id"],
        "molecule_pref_name": head["molecule_pref_name"] or "",
        "molecule_type": head["molecule_type"] or "",
        "max_phase": head["max_phase"] if head["max_phase"] is not None else "",
        "canonical_smiles": head["canonical_smiles"] or "",
        "standard_inchi_key": head["standard_inchi_key"] or "",
    }

    pot = [r for r in recs if r["metric_role"] == "potency"]
    eff = [r for r in recs if r["metric_role"] == "efficacy"]
    kin = [r for r in recs if r["metric_role"] == "kinetic"]
    unc = [r for r in recs if r["metric_role"] == "unclassified"]

    # ---------------- 效力 ----------------
    # 删失值（'>' 表示最高浓度未达 EC50）不参与汇总，否则会把弱化合物算强
    quant = [r for r in pot if r["standard_relation"] == "=" and r["value_nm"]]
    cens_gt = [r for r in pot if r["standard_relation"] == ">"]
    cens_lt = [r for r in pot if r["standard_relation"] == "<"]
    nms = sorted(r["value_nm"] for r in quant)
    pas = sorted(r["pactivity"] for r in quant if r["pactivity"] is not None)

    out["n_potency_records"] = len(pot)
    out["potency_types"] = "; ".join(f"{k}:{v}" for k, v in
                                     Counter(r["standard_type"] for r in pot).most_common())
    out["potency_nm_min"] = _stat(nms, min)
    out["potency_nm_median"] = _stat(nms, statistics.median)
    out["potency_nm_max"] = _stat(nms, max)
    out["pactivity_max"] = _stat(pas, max)
    out["pactivity_median"] = _stat(pas, statistics.median)
    out["pactivity_min"] = _stat(pas, min)
    out["pactivity_spread"] = round(max(pas) - min(pas), 2) if len(pas) >= 2 else ""
    out["n_potency_censored_gt"] = len(cens_gt)
    gt_nm = [r["value_nm"] for r in cens_gt if r["value_nm"]]
    out["potency_censored_gt_min_nm"] = _stat(sorted(gt_nm), min)
    out["n_potency_censored_lt"] = len(cens_lt)
    out["potency_n_assays"] = len({r["assay_chembl_id"] for r in pot})
    out["potency_n_docs"] = len({r["doc_chembl_id"] for r in pot if r["doc_chembl_id"]})
    out["potency_records"] = json.dumps(
        [{"assay": r["assay_chembl_id"], "type": r["standard_type"],
          "relation": r["standard_relation"], "value_nm": r["value_nm"],
          "pactivity": r["pactivity"], "doc": r["doc_chembl_id"], "rule": r["role_rule"]}
         for r in pot], ensure_ascii=False) if pot else ""

    # ---------------- 效能 ----------------
    fold = [r for r in eff if r["metric_scale"] == "fold" and r["standard_value"] is not None]
    pct = [r for r in eff if r["metric_scale"] == "percent" and r["standard_value"] is not None]
    out["n_efficacy_records"] = len(eff)

    def _best(items):
        if not items:
            return "", ""
        b = max(items, key=lambda r: float(r["standard_value"]))
        return round(float(b["standard_value"]), 3), f"{b['standard_type']}@{b['assay_chembl_id']}"

    out["efficacy_fold_max"], out["efficacy_fold_max_source"] = _best(fold)
    out["n_efficacy_fold"] = len(fold)
    out["efficacy_pct_max"], out["efficacy_pct_max_source"] = _best(pct)
    out["n_efficacy_pct"] = len(pct)
    out["efficacy_n_assays"] = len({r["assay_chembl_id"] for r in eff})
    out["efficacy_n_docs"] = len({r["doc_chembl_id"] for r in eff if r["doc_chembl_id"]})
    out["efficacy_records"] = json.dumps(
        [{"assay": r["assay_chembl_id"], "type": r["standard_type"],
          "scale": r["metric_scale"], "value": r["standard_value"],
          "units": r["standard_units"], "doc": r["doc_chembl_id"], "rule": r["role_rule"]}
         for r in eff], ensure_ascii=False) if eff else ""

    # ---------------- 酶动力学 / 未归类 ----------------
    out["n_kinetic_records"] = len(kin)
    out["kinetic_records"] = json.dumps(
        [{"assay": r["assay_chembl_id"], "type": r["standard_type"],
          "value": r["standard_value"], "units": r["standard_units"]} for r in kin],
        ensure_ascii=False) if kin else ""
    out["n_unclassified_records"] = len(unc)
    out["unclassified_records"] = json.dumps(
        [{"assay": r["assay_chembl_id"], "type": r["standard_type"],
          "value": r["standard_value"], "units": r["standard_units"],
          "reason": r["role_rule"]} for r in unc], ensure_ascii=False) if unc else ""

    # ---------------- 证据 ----------------
    assays = {r["assay_chembl_id"] for r in recs}
    docs = {r["doc_chembl_id"] for r in recs if r["doc_chembl_id"]}
    years = sorted({str(r["doc_year"]) for r in recs if r["doc_year"]})
    out["n_activities"] = len(recs)
    out["n_assays"] = len(assays)
    out["n_docs"] = len(docs)
    out["doc_chembl_ids"] = json.dumps(sorted(docs), ensure_ascii=False) if docs else ""
    out["doc_years"] = "; ".join(years)
    out["src_short_names"] = "; ".join(sorted({r["src_short_name"] for r in recs
                                               if r["src_short_name"]}))
    out["n_potential_duplicate"] = sum(1 for r in recs if r["potential_duplicate"])

    # 一致性只在「同一分子有 ≥2 个独立 assay 的可定量效力值」时才有意义
    spread = out["pactivity_spread"]
    if len(pas) < 2 or out["potency_n_assays"] < 2:
        consistency = "不适用（可定量效力值不足 2 个 assay）"
    elif spread <= SPREAD_TIGHT:
        consistency = f"一致（pActivity 极差 {spread}，3 倍以内）"
    elif spread <= SPREAD_LOOSE:
        consistency = f"较一致（pActivity 极差 {spread}，10 倍以内）"
    else:
        consistency = f"分歧（pActivity 极差 {spread}，超过 10 倍）"
    out["evidence_consistency"] = consistency

    # 证据强度是对「支撑量」的描述，不是活性强弱，也不参与效力/效能
    if out["n_docs"] >= 2 and out["n_assays"] >= 2 and not consistency.startswith("分歧"):
        level = "强"
    elif out["n_assays"] >= 2 or out["n_docs"] >= 2:
        level = "中"
    else:
        level = "弱"
    out["evidence_level"] = level

    notes = []
    if out["n_docs"] <= 1:
        notes.append("仅 1 篇文献支撑")
    if cens_gt:
        notes.append(f"{len(cens_gt)} 条效力为删失值（>，最高浓度未达 EC50），未计入汇总")
    if not quant and not eff:
        notes.append("无可定量的效力或效能值")
    if out["n_potential_duplicate"]:
        notes.append(f"{out['n_potential_duplicate']} 条被 ChEMBL 标为 potential_duplicate")
    if unc:
        notes.append(f"{len(unc)} 条指标未归类（见 unclassified_records）")
    out["evidence_note"] = "；".join(notes)
    return out


def write_csv(rows: list, path: Path, columns: list) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in columns})


def write_report(mols: list, acts: list, assays: list, dropped: int,
                 path: Path, src: Path, db_path: Path, version: str) -> None:
    L: list = []
    L.append("# Step1_04 GCK 激活剂活性提取")
    L.append("")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db_path}`")
    L.append(f"- 运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- 输入：`{src.name}`（Step1_03 终表）")
    L.append(f"- 筛选：`final_category == \"{ACTIVATION_CATEGORY}\"` "
             "且 `target_identity_suspect == \"FALSE\"`")
    L.append(f"- 命中激活 assay：**{len(assays)}** 个"
             f"（另有 {dropped} 个因靶点身份存疑被排除）")
    L.append(f"- 提取 activity：**{len(acts):,}** 条，涉及分子 **{len(mols):,}** 个")
    L.append("")
    L.append("> 效力、效能、证据**三项分开记录，不合成总分**；"
             "本步骤**不对 activity 排序**，输出按 `molecule_chembl_id` 字典序排列。")
    L.append("")

    L.append("## 指标分流")
    L.append("")
    L.append("直接按 `standard_type` 归并会算错，实测到三个坑，处理方式如下。")
    L.append("")
    role_cnt = Counter(r["metric_role"] for r in acts)
    L.append("| 角色 | activity 数 | 含义 |")
    L.append("| --- | ---: | --- |")
    role_desc = {
        "potency": "效力：浓度型指标，越小越强",
        "efficacy": "效能：激活幅度，越大越强",
        "kinetic": "酶动力学：机制读数（S0.5/Km/Vmax/Kcat），**不计入效力或效能**",
        "unclassified": "未归类：自由文本指标且描述未命中提升门槛，如实保留不汇总",
    }
    for role in ("potency", "efficacy", "kinetic", "unclassified"):
        if role_cnt.get(role):
            L.append(f"| `{role}` | {role_cnt[role]:,} | {role_desc[role]} |")
    L.append("")
    L.append("### 各 standard_type 的去向")
    L.append("")
    L.append("| standard_type | units | 角色 | 尺度 | activity 数 | 这个指标测的是什么 |")
    L.append("| --- | --- | --- | --- | ---: | --- |")
    combo = Counter((r["standard_type"], r["standard_units"] or "—",
                     r["metric_role"], r["metric_scale"] or "—") for r in acts)
    for (st, su, role, scale), n in combo.most_common():
        meaning = STANDARD_TYPE_MEANING.get(st, "—")
        L.append(f"| `{st}` | {su} | {role} | {scale} | {n:,} | {meaning} |")
    L.append("")
    L.append("要点：")
    L.append("")
    L.append("1. **`S0.5` / `S50`（mM）不是化合物效力**，是酶对葡萄糖的半饱和常数。"
             "GKA 的作用正是把它从野生型的 ~7 mM 拉低到 0.5–2 mM，属机制读数。"
             "误当效力会得出「EC50 = 0.6 mM」这种量级完全错误的结论。")
    L.append("2. **`standard_units` 不可信**：`Emax` 标着 `%` 而数值是 1.04（倍数），"
             "同一 `type + units` 内值域跨 0.63–1454。尺度按 (assay, type) 分组的值域判定，"
             "单位只作参考。")
    L.append("3. **`Activity` / `Ratio` 是自由文本指标**，默认不跨实验汇总；"
             "只有 assay 描述命中显式模式时才做**带证据的受限提升**，"
             "命中的原文写在 `role_rule` / `potency_records` / `efficacy_records` 里可复核。")
    L.append("")

    promoted = [r for r in acts if "受限提升" in (r["role_rule"] or "")]
    L.append(f"### 受限提升的记录（{len(promoted):,} 条）")
    L.append("")
    if promoted:
        L.append("| standard_type | units | 提升为 | activity 数 | 命中的描述原文 |")
        L.append("| --- | --- | --- | ---: | --- |")
        pc = Counter()
        for r in promoted:
            hit = (r["role_rule"].split("：「")[-1].rstrip("」")
                   if "：「" in r["role_rule"] else "")
            pc[(r["standard_type"], r["standard_units"] or "—",
                r["metric_role"] + "/" + r["metric_scale"], hit)] += 1
        for (st, su, tgt, hit), n in pc.most_common():
            L.append(f"| `{st}` | {su} | {tgt} | {n:,} | {hit.replace('|', chr(92) + '|')} |")
    else:
        L.append("无。")
    L.append("")

    L.append("## 覆盖情况")
    L.append("")
    n_pot = sum(1 for m in mols if m["pactivity_max"] != "")
    n_eff = sum(1 for m in mols if m["n_efficacy_records"])
    n_both = sum(1 for m in mols if m["pactivity_max"] != "" and m["n_efficacy_records"])
    n_none = sum(1 for m in mols if m["pactivity_max"] == "" and not m["n_efficacy_records"])
    L.append("| 项目 | 分子数 | 占比 |")
    L.append("| --- | ---: | ---: |")
    tot = len(mols)
    for name, n in [("有可定量效力（pActivity）", n_pot), ("有效能值", n_eff),
                    ("效力与效能兼有", n_both), ("两者皆无", n_none)]:
        L.append(f"| {name} | {n:,} | {100.0 * n / tot:.1f}% |")
    L.append("")
    L.append("「两者皆无」不是错误，是如实记录：这些分子在激活 assay 里只留下了"
             "酶动力学参数或自由文本指标。")
    L.append("")

    L.append("## 效力分布")
    L.append("")
    pas = [float(m["pactivity_max"]) for m in mols if m["pactivity_max"] != ""]
    if pas:
        bins = [(9, 99, "pActivity ≥ 9（EC50 ≤ 1 nM）"),
                (8, 9, "8 ≤ pActivity < 9（1–10 nM）"),
                (7, 8, "7 ≤ pActivity < 8（10–100 nM）"),
                (6, 7, "6 ≤ pActivity < 7（0.1–1 uM）"),
                (5, 6, "5 ≤ pActivity < 6（1–10 uM）"),
                (0, 5, "pActivity < 5（> 10 uM）")]
        L.append("按每个分子的**最优** pActivity 统计。分箱是描述性的，不构成排序。")
        L.append("")
        L.append("| 区间 | 分子数 |")
        L.append("| --- | ---: |")
        for lo, hi, label in bins:
            L.append(f"| {label} | {sum(1 for v in pas if lo <= v < hi):,} |")
        L.append("")
        L.append(f"pActivity 中位数 {statistics.median(pas):.2f}，"
                 f"范围 {min(pas):.2f} – {max(pas):.2f}。")
        L.append("")
    n_cens = sum(1 for m in mols if m["n_potency_censored_gt"])
    L.append(f"另有 **{n_cens:,}** 个分子存在删失效力值（`standard_relation = '>'`，"
             "最高浓度仍未达 EC50）。这类值**未计入任何汇总统计**——"
             "当作等号值会把弱化合物算强，但也不能直接丢弃，"
             "计数保留在 `n_potency_censored_gt`。")
    L.append("")

    L.append("## 效能分布")
    L.append("")
    L.append("fold 与 percent 是两种尺度，**分开统计，不做换算**。")
    L.append("")
    folds = [float(m["efficacy_fold_max"]) for m in mols if m["efficacy_fold_max"] != ""]
    pcts = [float(m["efficacy_pct_max"]) for m in mols if m["efficacy_pct_max"] != ""]
    L.append("| 尺度 | 分子数 | 最小 | 中位 | 最大 |")
    L.append("| --- | ---: | ---: | ---: | ---: |")
    for name, vals in [("fold activation", folds), ("% activation", pcts)]:
        if vals:
            L.append(f"| {name} | {len(vals):,} | {min(vals):.2f} | "
                     f"{statistics.median(vals):.2f} | {max(vals):.2f} |")
        else:
            L.append(f"| {name} | 0 | — | — | — |")
    L.append("")

    L.append("## 证据分布")
    L.append("")
    L.append("| 证据强度 | 分子数 | 判定条件 |")
    L.append("| --- | ---: | --- |")
    cond = {"强": "≥2 篇独立文献且 ≥2 个 assay，且结果不分歧",
            "中": "≥2 个 assay 或 ≥2 篇文献",
            "弱": "单一 assay 且单一文献"}
    lv = Counter(m["evidence_level"] for m in mols)
    for k in ("强", "中", "弱"):
        L.append(f"| {k} | {lv.get(k, 0):,} | {cond[k]} |")
    L.append("")
    L.append("> 证据强度描述的是**支撑量**，与活性强弱无关，不参与效力/效能。")
    L.append("")
    L.append("| 一致性 | 分子数 |")
    L.append("| --- | ---: |")
    for k, n in Counter(m["evidence_consistency"].split("（")[0]
                        for m in mols).most_common():
        L.append(f"| {k} | {n:,} |")
    L.append("")

    L.append("## 数据来源分布")
    L.append("")
    L.append("| 来源 | activity 数 |")
    L.append("| --- | ---: |")
    for k, n in Counter(r["src_short_name"] or "(空)" for r in acts).most_common():
        L.append(f"| {k} | {n:,} |")
    L.append("")
    L.append(f"涉及文献 **{len({r['doc_chembl_id'] for r in acts if r['doc_chembl_id']})}** 篇，"
             f"年份 {min(str(r['doc_year']) for r in acts if r['doc_year'])}–"
             f"{max(str(r['doc_year']) for r in acts if r['doc_year'])}。")
    L.append("")

    n_gluc = sum(1 for r in acts if r["glucose_mM_parsed"])
    L.append("## 附：从描述解析出的葡萄糖浓度")
    L.append("")
    L.append(f"葡萄糖浓度没有结构化到 `assay_parameters`，只写在描述自由文本里。"
             f"明细表的 `glucose_mM_parsed` 列做了正则解析，**{n_gluc:,} / {len(acts):,}** "
             "条 activity 解析到了值。低糖/高糖是区分 GKA 类型的关键条件，"
             "留给后续步骤使用；**本步骤不据此做任何判定**。")
    L.append("")
    gc = Counter(r["glucose_mM_parsed"] for r in acts if r["glucose_mM_parsed"])
    if gc:
        L.append("| 葡萄糖浓度 (mM) | activity 数 |")
        L.append("| ---: | ---: |")
        for k, n in sorted(gc.items(), key=lambda kv: float(kv[0])):
            L.append(f"| {k} | {n:,} |")
        L.append("")

    L.append("## 全部分子清单")
    L.append("")
    L.append("按 `molecule_chembl_id` 字典序，**不是按活性排序**。完整字段见同目录 CSV。")
    L.append("")
    L.append("| molecule | 名称 | 效力 pActivity(最优/中位) | EC50 最小(nM) | "
             "效能 fold / % | assay | 文献 | 证据 | 一致性 |")
    L.append("| --- | --- | --- | ---: | --- | ---: | ---: | --- | --- |")
    for m in mols:
        pa = (f"{m['pactivity_max']} / {m['pactivity_median']}"
              if m["pactivity_max"] != "" else "—")
        ef = "/".join([str(m["efficacy_fold_max"]) if m["efficacy_fold_max"] != "" else "—",
                       str(m["efficacy_pct_max"]) if m["efficacy_pct_max"] != "" else "—"])
        name = (m["molecule_pref_name"] or "—")[:28].replace("|", "\\|")
        L.append(f"| `{m['molecule_chembl_id']}` | {name} | {pa} | "
                 f"{m['potency_nm_min'] if m['potency_nm_min'] != '' else '—'} | {ef} | "
                 f"{m['n_assays']} | {m['n_docs']} | {m['evidence_level']} | "
                 f"{m['evidence_consistency'].split('（')[0]} |")
    L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(
        description="从 GCK 激活 assay 提取活性并按分子汇总效力/效能/证据。")
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--in-csv", type=Path, default=DEFAULT_IN,
                   help=f"Step1_03 的终表，默认 {DEFAULT_IN}")
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    if not args.db.is_file():
        print(f"错误：找不到数据库 {args.db}", file=sys.stderr)
        return 1
    if not args.in_csv.is_file():
        print(f"错误：找不到 Step1_03 的输出 {args.in_csv}，请先运行 Step1_03。",
              file=sys.stderr)
        return 1

    assays, dropped = load_activation_assays(args.in_csv)
    if not assays:
        print("错误：输入里没有符合条件的激活 assay。", file=sys.stderr)
        return 2
    assay_ids = [a["assay_chembl_id"] for a in assays]

    con = connect_readonly(args.db)
    version = chembl_version(con)
    print(f"数据库：{args.db}")
    print(f"ChEMBL 版本：{version}")
    print(f"输入：{args.in_csv}")
    print(f"激活 assay：{len(assay_ids)} 个（排除靶点身份存疑 {dropped} 个）\n")

    acts = fetch_activities(con, assay_ids)
    print(f"提取 activity：{len(acts):,} 条")
    rc = Counter(r["metric_role"] for r in acts)
    for role in ("potency", "efficacy", "kinetic", "unclassified"):
        if rc.get(role):
            print(f"  {role:<14s} {rc[role]:>6,}")

    by_mol: dict = defaultdict(list)
    for r in acts:
        by_mol[r["molecule_chembl_id"]].append(r)
    mols = [summarize_molecule(by_mol[k]) for k in sorted(by_mol)]

    print(f"\n分子：{len(mols):,} 个")
    print(f"  有可定量效力：{sum(1 for m in mols if m['pactivity_max'] != ''):,}")
    print(f"  有效能值：    {sum(1 for m in mols if m['n_efficacy_records']):,}")
    print(f"  证据强度分布：{dict(Counter(m['evidence_level'] for m in mols))}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    mol_csv = args.outdir / "Step1_04_GCK_Activator_Activity_Extraction.csv"
    act_csv = args.outdir / "Step1_04_GCK_Activator_Activities.csv"
    md_path = args.outdir / "Step1_04_GCK_Activator_Activity_Extraction.md"
    write_csv(mols, mol_csv, MOL_COLUMNS)
    write_csv(acts, act_csv, ACT_COLUMNS)
    write_report(mols, acts, assays, dropped, md_path, args.in_csv, args.db, version)

    print(f"\n分子汇总（主产物）：{mol_csv}")
    print(f"activity 明细：      {act_csv}")
    print(f"报告：               {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
