#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_05_GCK_Activator_Ranking_and_Candidate_Selection
======================================================

对 Step1_04 的 1,333 个分子做**方向判定 → 排除打标 → 排序 → 分层 → 骨架去冗余**，
输出一张可复核的 GKA 候选表。

本步骤只回答「谁是好的 GCK 激活剂」，**不判断能不能进脑**——物化性质与 BBB 留到 Step2。

输入
----
    ../Step1_04_GCK_Activator_Activity_Extraction/Step1_04_GCK_Activator_Activity_Extraction.csv
    ../Step1_04_GCK_Activator_Activity_Extraction/Step1_04_GCK_Activator_Activities.csv

方向判定：为什么不能对 fold 设全局阈值
--------------------------------------
Step1_04 的分子是「出现在激活 assay 里的分子」，不是激活剂。补方向这道门时，
最直觉的做法「fold ≤ 1.05 即无激活」是错的，实测到三种翻车方式：

1. **分母不是对照，而是参比激活剂。** `CHEMBL1825591` 的描述是
   "ratio of enzyme activation in compound **treated to Ro-28-1675**"，
   值域 -0.025 – 1.0：这里 1.0 表示「与参比一样强」，是**最好**的结果，
   基线是 0 不是 1。同类还有 `CHEMBL4272645`（relative to RO0281675）、
   `CHEMBL1056162/64`、`CHEMBL1038416/18`、`CHEMBL1058314`
   （relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-...）。
   按 1.05 判，这些 assay 里的分子会被整批误杀。

2. **percent 类的基线有 0 和 100 两种。** `CHEMBL2353000` 值域 0–127（基线 0），
   `CHEMBL4014755` 值域 452–784（基线 100），同为 `%max`，含义相反。

3. **Km / Vmax 类比值是机制读数，极性还相反。** `CHEMBL3095514` 的描述是
   "**decrease** in enzyme Km for glucose"，值 0.05–0.11——Km 降到 5–11%，
   是**强激活**；按 fold ≤ 1.05 判会被当成无激活排除。
   K 型 GKA 本来就不改 Vmax，`CHEMBL4427351` 的 Vmax ratio 0.54–1.1 同理不能判死。

因此方向判定的做法是：**逐 (assay, standard_type, scale) 组解析基线**，
规则写死在 `BASELINE_RULES` 的判定顺序里，解析结果连同证据句全量写进报告供复核；
解析不出来的组标 `不确定`，**不参与方向判定**（保险的一侧），不静默当成无激活。

只有 EC50、没有效能记录的分子**不排除**——EC50 出自激活读数的 assay，本身就是方向证据。

排序
----
**单一排序键 `pactivity_median`，不合成总分。** 中位数抗单点异常；`pactivity_max`
并列展示不参与排序。**效能与证据只作并列列**——证据等级衡量的是「在 ChEMBL 里被测了
多少次」，MK-0941、AZD-1656、PF-04991532 三个 phase 2 药的证据都是「弱」，
把证据折进排序等于系统性地把最像药的分子往后排。

输出
----
    Step1_05_GCK_Activator_Candidates.csv                        一行一个分子（主产物）
    Step1_05_GCK_Activator_Ranking_and_Candidate_Selection.md    报告

用法
----
    python3 Step1_05_GCK_Activator_Ranking_and_Candidate_Selection.py
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "Step1_04_GCK_Activator_Activity_Extraction"
DEFAULT_MOL = PREV / "Step1_04_GCK_Activator_Activity_Extraction.csv"
DEFAULT_ACT = PREV / "Step1_04_GCK_Activator_Activities.csv"
DEFAULT_PREV_MD = PREV / "Step1_04_GCK_Activator_Activity_Extraction.md"

# ---------------------------------------------------------------------------
# 一、方向判定：基线解析
# ---------------------------------------------------------------------------
# 判定顺序即优先级，命中即停。每条都要能从描述原文找到依据。

# 机制读数：Km / Vmax / S0.5 类比值。极性与 fold 相反或对 K 型激活剂不敏感，
# 一律**不参与方向判定**（CLAUDE.md：这些是酶的性质，不是化合物效力）。
RE_KINETIC_LIKE = re.compile(
    r"\b(?:Km|Vmax|S0\.5|S50)\b", re.I)

# 分母是参比激活剂而非对照 → 基线 0，值是「达到参比效果的几成」
RE_REFERENCE = re.compile(
    r"(?:relative to|compared to|treated to)\s+"
    r"(?!(?:the\s+)?(?:untreated\s+|vehicle\s+)?control\b)"
    r"([A-Za-z0-9(\[][^;]{2,90})", re.I)

# 分母是未处理对照 → fold 基线 1 / percent 基线 100（percent 还要看值域）
RE_CONTROL = re.compile(
    r"relative to\s+(?:the\s+)?(?:untreated\s+|vehicle\s+)?control"
    r"|to untreated control|versus control", re.I)

# 判定为「有激活 / 无激活」的阈值，按基线分别给
THRESH = {
    # (scale, baseline): (激活阈, 无激活阈)
    ("fold", "control"): (1.05, 0.95),      # 基线 1
    ("percent", "control"): (105.0, 95.0),  # 基线 100
    ("fold", "reference"): (0.10, 0.0),     # 基线 0，1.0 = 与参比等效
    ("percent", "reference"): (10.0, 0.0),  # 基线 0，100 = 与参比等效
    ("percent", "max"): (10.0, 0.0),        # 基线 0，100 = 该实验内的满标
}

# ---------------------------------------------------------------------------
# 二、分层与排序
# ---------------------------------------------------------------------------
# 阈值由阳性对照标定，不是按整体分布拍的。最初取 7.0 / 6.5（EC50 ≤ 100 / 316 nM），
# 自检只过 4/6：Ro-28-1675 跨 4 个实验的 EC50 是 127–690 nM（pAct 中位 6.39），
# Piraglitin 是 364–6320 nM（中位 6.145）——**已知临床 GKA 本来就在几百 nM 这一档**，
# 原阈值把临床药从中间劈开了。改按对照的实际下限（6.145）设 B 层门槛。
TIER_A_PACT = 6.5    # A 层：效力门槛 + 需要正向效能佐证（双证）
TIER_B_PACT = 6.0    # B 层：效力门槛（缺效能读数不代表更差，只是只测了一个轴）
SPREAD_DIVERGENT = 1.0   # pActivity 极差 > 1（10 倍）算分歧

# 阳性对照：已知临床/参比 GKA，规则跑完必须落进 A/B 层
POSITIVE_CONTROLS = {
    "CHEMBL1096435": "Ro-281675",
    "CHEMBL1783734": "PIRAGLIATIN",
    "CHEMBL2165615": "NERIGLIATIN",
    "CHEMBL2165620": "PF-04991532",
    "CHEMBL3219124": "AZD-1656",
    "CHEMBL3580737": "MK-0941",
}

OUT_COLUMNS = [
    # --- 身份 ---
    "molecule_chembl_id", "molecule_pref_name", "max_phase", "canonical_smiles",
    # --- 方向判定 ---
    "direction", "direction_evidence",
    "n_direction_activation", "n_direction_no_activation", "n_direction_undecidable",
    # --- 入选与排除 ---
    "included", "exclude_reason", "flags",
    # --- 排序与分层 ---
    "tier", "rank_overall", "rank_in_tier",
    "pactivity_median", "pactivity_max", "pactivity_spread", "potency_nm_min",
    # --- 效能（并列展示，不参与排序）---
    "efficacy_fold_max", "efficacy_pct_max", "n_efficacy_records",
    # --- 证据（并列展示，不参与排序）---
    "evidence_level", "evidence_consistency", "n_assays", "n_docs",
    "n_potency_censored_gt", "n_potential_duplicate",
    # --- 骨架 ---
    "murcko_scaffold", "scaffold_cluster_size", "is_scaffold_representative",
]


def fnum(s):
    """空串按缺失处理——空值是事实，不能当 0。"""
    if s is None or s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_csv(path: Path) -> list:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_provenance(md_path: Path) -> tuple:
    """从 Step1_04 报告头部取 ChEMBL 版本与数据库路径，保持链路可追溯。"""
    version, db = "未知", "未知"
    if md_path.is_file():
        for line in md_path.read_text(encoding="utf-8").splitlines()[:12]:
            if line.startswith("- ChEMBL 版本："):
                version = line.split("：", 1)[1].strip().strip("*")
            elif line.startswith("- 数据库文件："):
                db = line.split("：", 1)[1].strip().strip("`")
    return version, db


# ---------------------------------------------------------------------------
# 基线解析
# ---------------------------------------------------------------------------
def resolve_baselines(acts: list) -> dict:
    """逐 (assay, standard_type, scale) 组解析基线。

    返回 {(assay, type, scale): info}，info 含 baseline / 证据 / 值域。
    baseline 取值：control | reference | max | kinetic | uncertain
    只有前三种参与方向判定。
    """
    groups = defaultdict(list)
    for r in acts:
        if r["metric_role"] != "efficacy":
            continue
        key = (r["assay_chembl_id"], r["standard_type"], r["metric_scale"])
        groups[key].append(r)

    out = {}
    for key, rows in groups.items():
        _, _, scale = key
        desc = rows[0]["assay_description"] or ""
        vals = [v for v in (fnum(r["standard_value"]) for r in rows) if v is not None]
        vmin, vmax = (min(vals), max(vals)) if vals else (None, None)
        info = {"n": len(rows), "vmin": vmin, "vmax": vmax,
                "desc": desc, "evidence": "", "note": ""}

        m_kin = RE_KINETIC_LIKE.search(desc)
        m_ref = RE_REFERENCE.search(desc)
        m_ctl = RE_CONTROL.search(desc)

        if m_kin:
            # Km/Vmax/S0.5 类：极性相反或对 K 型激活剂不敏感，不判方向
            info["baseline"] = "kinetic"
            info["evidence"] = m_kin.group(0)
            info["note"] = "酶动力学读数（Km/Vmax/S0.5），极性与 fold 不一致，不参与方向判定"
        elif m_ref:
            info["baseline"] = "reference"
            info["evidence"] = m_ref.group(0)[:70]
            info["note"] = "分母是参比激活剂 → 基线 0，值是达到参比效果的比例"
        elif m_ctl and scale == "fold":
            info["baseline"] = "control"
            info["evidence"] = m_ctl.group(0)
            info["note"] = "分母是未处理对照 → 基线 1"
        elif m_ctl and scale == "percent":
            # percent + 对照：描述本身定不下基线是 0 还是 100，用值域裁定。
            # 依据是「% of max」装不下远超 100 的值。
            if vmax is not None and vmax > 150:
                info["baseline"] = "control"
                info["evidence"] = m_ctl.group(0)
                info["note"] = f"值域上限 {vmax:g} > 150，装不下「占满标的百分数」→ 基线 100"
            elif vmin is not None and vmin == 0 and vmax is not None and vmax <= 150:
                info["baseline"] = "max"
                info["evidence"] = m_ctl.group(0)
                info["note"] = f"值域 {vmin:g}–{vmax:g} 且触 0 → 占满标的百分数，基线 0"
            else:
                info["baseline"] = "uncertain"
                info["note"] = (f"percent + 对照，值域 {vmin:g}–{vmax:g} 两种基线都讲得通，"
                                "不参与方向判定")
        else:
            info["baseline"] = "uncertain"
            info["note"] = "描述未写明分母，不参与方向判定"
        out[key] = info
    return out


def judge_direction(rows: list, baselines: dict) -> dict:
    """按分子判方向。保守：判不了就不排除。"""
    n_act = n_no = n_und = 0
    ev_act, ev_no = [], []
    for r in rows:
        if r["metric_role"] != "efficacy":
            continue
        key = (r["assay_chembl_id"], r["standard_type"], r["metric_scale"])
        info = baselines.get(key)
        v = fnum(r["standard_value"])
        if info is None or v is None or info["baseline"] in ("kinetic", "uncertain"):
            n_und += 1
            continue
        th = THRESH.get((r["metric_scale"], info["baseline"]))
        if th is None:
            n_und += 1
            continue
        hi, lo = th
        tag = f"{r['assay_chembl_id']} {r['standard_type']}={v:g}（基线{info['baseline']}）"
        if v >= hi:
            n_act += 1
            ev_act.append(tag)
        elif v <= lo:
            n_no += 1
            ev_no.append(tag)
        else:
            n_und += 1

    # 非删失的效力值本身就是方向证据：EC50 出自激活读数的 assay
    has_quant_potency = any(
        r["metric_role"] == "potency" and r["standard_relation"] != ">"
        and fnum(r["standard_value"]) is not None for r in rows)

    if n_act and n_no:
        direction = "conflict"
    elif n_act:
        direction = "activation"
    elif n_no:
        direction = "no_activation"
    elif has_quant_potency:
        direction = "activation_by_potency"
    else:
        direction = "unknown"

    ev = []
    if ev_act:
        ev.append("激活：" + "；".join(ev_act[:3]))
    if ev_no:
        ev.append("未激活：" + "；".join(ev_no[:3]))
    if not ev and direction == "activation_by_potency":
        ev.append("无可判方向的效能读数；效力值出自激活读数的 assay")
    if not ev and direction == "unknown":
        ev.append("既无可判方向的效能读数，也无非删失效力值")
    return {"direction": direction, "direction_evidence": "｜".join(ev),
            "n_direction_activation": n_act, "n_direction_no_activation": n_no,
            "n_direction_undecidable": n_und}


# ---------------------------------------------------------------------------
# 骨架
# ---------------------------------------------------------------------------
def murcko(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# 组装
# ---------------------------------------------------------------------------
def build_rows(mols: list, acts: list, baselines: dict) -> list:
    by_mol = defaultdict(list)
    for r in acts:
        by_mol[r["molecule_chembl_id"]].append(r)

    rows = []
    for m in mols:
        mid = m["molecule_chembl_id"]
        d = judge_direction(by_mol.get(mid, []), baselines)
        pa_med = fnum(m["pactivity_median"])
        spread = fnum(m["pactivity_spread"])
        n_cens = int(m["n_potency_censored_gt"] or 0)
        n_dup = int(m["n_potential_duplicate"] or 0)

        # --- 排除（不删行，写明理由）---
        reasons = []
        if d["direction"] == "no_activation":
            reasons.append("方向判定为无激活")
        if pa_med is None and n_cens:
            reasons.append(f"效力仅有删失值（{n_cens} 条 >），不进榜")

        # --- 打标（进榜但提醒）---
        flags = []
        if d["direction"] == "conflict":
            flags.append("方向冲突")
        if spread is not None and spread > SPREAD_DIVERGENT:
            flags.append(f"跨 assay 效力分歧（极差 {spread:g}）")
        if n_dup:
            flags.append(f"{n_dup} 条 potential_duplicate")
        if n_cens and pa_med is not None:
            flags.append(f"另有 {n_cens} 条删失效力值未计入")

        rows.append({
            "molecule_chembl_id": mid,
            "molecule_pref_name": m["molecule_pref_name"],
            "max_phase": m["max_phase"],
            "canonical_smiles": m["canonical_smiles"],
            "included": "FALSE" if reasons else "TRUE",
            "exclude_reason": "；".join(reasons),
            "flags": "；".join(flags),
            "pactivity_median": m["pactivity_median"],
            "pactivity_max": m["pactivity_max"],
            "pactivity_spread": m["pactivity_spread"],
            "potency_nm_min": m["potency_nm_min"],
            "efficacy_fold_max": m["efficacy_fold_max"],
            "efficacy_pct_max": m["efficacy_pct_max"],
            "n_efficacy_records": m["n_efficacy_records"],
            "evidence_level": m["evidence_level"],
            "evidence_consistency": m["evidence_consistency"],
            "n_assays": m["n_assays"], "n_docs": m["n_docs"],
            "n_potency_censored_gt": m["n_potency_censored_gt"],
            "n_potential_duplicate": m["n_potential_duplicate"],
            "_pa_med": pa_med,
            "_pa_max": fnum(m["pactivity_max"]),
            **d,
        })
    return rows


def assign_tiers(rows: list) -> None:
    """分层 + 层内排序。排序键单一：pactivity_median。"""
    for r in rows:
        if r["included"] != "TRUE":
            r["tier"] = ""
            continue
        pa = r["_pa_med"]
        if pa is None:
            r["tier"] = "无效力（只有效能/机制读数）"
            continue
        pos_eff = r["n_direction_activation"] > 0
        if pa >= TIER_A_PACT and pos_eff:
            r["tier"] = "A"
        elif pa >= TIER_B_PACT:
            r["tier"] = "B"
        else:
            r["tier"] = "C"

    ranked = [r for r in rows if r["included"] == "TRUE" and r["_pa_med"] is not None]
    ranked.sort(key=lambda r: (-r["_pa_med"], -(r["_pa_max"] or 0),
                               r["molecule_chembl_id"]))
    for i, r in enumerate(ranked, 1):
        r["rank_overall"] = i
    per_tier = defaultdict(int)
    for r in ranked:
        per_tier[r["tier"]] += 1
        r["rank_in_tier"] = per_tier[r["tier"]]
    for r in rows:
        r.setdefault("rank_overall", "")
        r.setdefault("rank_in_tier", "")


def assign_scaffolds(rows: list) -> None:
    for r in rows:
        r["murcko_scaffold"] = murcko(r["canonical_smiles"])
    clusters = defaultdict(list)
    for r in rows:
        if r["murcko_scaffold"]:
            clusters[r["murcko_scaffold"]].append(r)
    for scaf, members in clusters.items():
        for r in members:
            r["scaffold_cluster_size"] = len(members)
        # 代表：入选且效力最优者；没有可定量效力的簇不设代表
        cands = [r for r in members
                 if r["included"] == "TRUE" and r["_pa_med"] is not None]
        for r in members:
            r["is_scaffold_representative"] = "FALSE"
        if cands:
            best = min(cands, key=lambda r: (-r["_pa_med"], r["molecule_chembl_id"]))
            best["is_scaffold_representative"] = "TRUE"
    for r in rows:
        r.setdefault("scaffold_cluster_size", "")
        r.setdefault("is_scaffold_representative", "")


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------
def write_csv(rows: list, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_report(rows, acts, baselines, path, mol_csv, act_csv, version, db) -> None:
    L = []
    inc = [r for r in rows if r["included"] == "TRUE"]
    ranked = sorted((r for r in inc if r["_pa_med"] is not None),
                    key=lambda r: r["rank_overall"])

    L.append("# Step1_05 GCK 激活剂排序与候选选择")
    L.append("")
    L.append(f"- ChEMBL 版本：**{version}**")
    L.append(f"- 数据库文件：`{db}`")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 输入：`{mol_csv.name}`、`{act_csv.name}`（Step1_04 产物）")
    L.append(f"- 输入分子：**{len(rows):,}** 个 → 入选 **{len(inc):,}** 个，"
             f"排除 **{len(rows) - len(inc):,}** 个")
    L.append("")
    L.append("> 本步骤只回答「谁是好的 GCK 激活剂」，**不判断能不能进脑**——"
             "物化性质与 BBB 留到 Step2。")
    L.append("> 排序键单一（`pactivity_median`），**效能与证据只并列展示，不合成总分**。")
    L.append("")

    # --- 方向判定 ---
    L.append("## 一、方向判定")
    L.append("")
    L.append("Step1_04 的分子是「出现在激活 assay 里的分子」，不是激活剂。这里补上这道门。")
    L.append("")
    dc = Counter(r["direction"] for r in rows)
    dd = {
        "activation": "效能读数显示酶活升高",
        "activation_by_potency": "无可判方向的效能读数，但有非删失效力值——"
                                 "EC50 出自激活读数的 assay，本身即方向证据",
        "conflict": "同一分子既有激活也有降低的读数，**进榜但打标**",
        "no_activation": "效能读数显示酶活未升高，**排除**",
        "unknown": "既无可判方向的效能读数，也无非删失效力值",
    }
    L.append("| direction | 分子数 | 含义 | 处理 |")
    L.append("| --- | ---: | --- | --- |")
    for k in ("activation", "activation_by_potency", "conflict", "no_activation", "unknown"):
        if dc.get(k):
            act = "排除" if k == "no_activation" else "保留"
            L.append(f"| `{k}` | {dc[k]:,} | {dd[k]} | {act} |")
    L.append("")
    L.append("### 基线解析表（方向判定的依据，逐组可复核）")
    L.append("")
    L.append("直觉做法「fold ≤ 1.05 即无激活」是错的：**比值型读数的分母不统一**。"
             "下表是逐 (assay, standard_type, scale) 组的解析结果，"
             "`baseline` 为 `kinetic` / `uncertain` 的组**不参与方向判定**。")
    L.append("")
    L.append("| assay | type | scale | n | 值域 | baseline | 判定依据 |")
    L.append("| --- | --- | --- | ---: | --- | --- | --- |")
    order = {"control": 0, "reference": 1, "max": 2, "kinetic": 3, "uncertain": 4}
    for key in sorted(baselines, key=lambda k: (order[baselines[k]["baseline"]],
                                                -baselines[k]["n"])):
        i = baselines[key]
        rng = f"{i['vmin']:g} – {i['vmax']:g}" if i["vmin"] is not None else "—"
        note = i["note"].replace("|", chr(92) + "|")
        L.append(f"| `{key[0]}` | {key[1]} | {key[2]} | {i['n']} | {rng} | "
                 f"**{i['baseline']}** | {note} |")
    L.append("")
    nb = Counter(i["baseline"] for i in baselines.values())
    L.append(f"共 {len(baselines)} 组：" + "、".join(
        f"`{k}` {n} 组" for k, n in nb.most_common()) + "。")
    L.append("")
    L.append("三个必须逐组读描述的实例：")
    L.append("")
    L.append("1. `CHEMBL1825591` 的分母是**参比激活剂 Ro-28-1675**，值域 -0.025 – 1.0，"
             "其中 1.0 表示「与参比等效」是最好的结果。按 1.05 判会把整组误杀。")
    L.append("2. `CHEMBL2353000`（%max，0–127，基线 0）与 `CHEMBL4014755`"
             "（%max，452–784，基线 100）同为 percent，含义相反。")
    L.append("3. `CHEMBL3095514` 的描述是 \"**decrease** in enzyme Km for glucose\"，"
             "值 0.05–0.11 表示 Km 降到 5–11%，是**强激活**；"
             "按 fold ≤ 1.05 判会被当成无激活排除。Km/Vmax 类一律不参与方向判定。")
    L.append("")

    # --- 排除与打标 ---
    L.append("## 二、排除与打标")
    L.append("")
    L.append("不静默丢行：排除的分子留在 CSV 里，`included = FALSE` 且写明 `exclude_reason`。")
    L.append("")
    ec = Counter()
    for r in rows:
        for x in filter(None, r["exclude_reason"].split("；")):
            ec[re.sub(r"\d+", "N", x)] += 1
    L.append("| 排除理由 | 分子数 |")
    L.append("| --- | ---: |")
    for k, n in ec.most_common():
        L.append(f"| {k} | {n:,} |")
    if not ec:
        L.append("| （无） | 0 |")
    L.append("")
    fc = Counter()
    for r in rows:
        for x in filter(None, r["flags"].split("；")):
            fc[re.sub(r"[\d.]+", "N", x)] += 1
    L.append("| 打标（进榜但提醒） | 分子数 |")
    L.append("| --- | ---: |")
    for k, n in fc.most_common():
        L.append(f"| {k} | {n:,} |")
    L.append("")

    # --- 分层 ---
    L.append("## 三、分层与排序")
    L.append("")
    L.append(f"排序键 **`pactivity_median`**（单一，不合成总分）。"
             f"用中位数不用最大值：多数分子只有一个 assay 时两者相等，"
             f"多 assay 时中位数抗单点异常。")
    L.append("")
    tc = Counter(r["tier"] for r in inc)
    td = {"A": f"pActivity 中位 ≥ {TIER_A_PACT}（EC50 ≤ ~320 nM）**且**有正向效能佐证",
          "B": f"pActivity 中位 ≥ {TIER_B_PACT}（EC50 ≤ 1 µM）",
          "C": "其余有可定量效力的",
          "无效力（只有效能/机制读数）": "有激活证据但没有可定量效力值，无法排序"}
    L.append("| 层 | 分子数 | 条件 |")
    L.append("| --- | ---: | --- |")
    for k in ("A", "B", "C", "无效力（只有效能/机制读数）"):
        if tc.get(k):
            L.append(f"| {k} | {tc[k]:,} | {td[k]} |")
    L.append("")
    L.append("**阈值是被阳性对照标定出来的，不是按整体分布拍的。** 最初取 7.0 / 6.5"
             "（EC50 ≤ 100 / 316 nM），自检只过 4/6——Ro-28-1675 跨 4 个实验的 EC50 是"
             " 127–690 nM（pAct 中位 6.39），Piraglitin 是 364–6320 nM（中位 6.145）。"
             "**已知临床 GKA 本来就在几百 nM 这一档**，原阈值把临床药从中间劈开了，"
             "说明阈值错了而不是药不好。现按对照的实际下限重设。")
    L.append("")
    L.append("**B 层不等于比 A 层差**：A 要求效力与效能双证，落在 B 的分子里有一部分"
             "（如 MK-0941、AZD-1656）只是 ChEMBL 里没有效能读数，只测了效力这一个轴。")
    L.append("")

    # --- 骨架 ---
    L.append("## 四、骨架去冗余")
    L.append("")
    clusters = Counter(r["murcko_scaffold"] for r in rows if r["murcko_scaffold"])
    reps = sum(1 for r in rows if r["is_scaffold_representative"] == "TRUE")
    L.append(f"全部 {len(rows):,} 个分子归入 **{len(clusters):,}** 个 Murcko 骨架，"
             f"最大一簇 **{max(clusters.values())}** 个分子，"
             f"单例骨架 **{sum(1 for v in clusters.values() if v == 1):,}** 个，"
             f"前 20 个骨架覆盖 **{sum(n for _, n in clusters.most_common(20)):,}** 个分子。")
    L.append("")
    L.append(f"**加列不删行**：`murcko_scaffold` / `scaffold_cluster_size` / "
             f"`is_scaffold_representative`（簇内入选且效力最优者，共 {reps:,} 个）。"
             "不做骨架聚类，top-50 会是两篇 SAR 论文的同系物列表——"
             "看着 50 个，其实 2 个化学起点。")
    L.append("")
    L.append("A/B 层里骨架最集中的 10 簇：")
    L.append("")
    ab = [r for r in inc if r["tier"] in ("A", "B")]
    abc = Counter(r["murcko_scaffold"] for r in ab if r["murcko_scaffold"])
    L.append("| 骨架 | A/B 层分子数 | 代表分子 |")
    L.append("| --- | ---: | --- |")
    for scaf, n in abc.most_common(10):
        rep = next((r["molecule_chembl_id"] for r in ab
                    if r["murcko_scaffold"] == scaf
                    and r["is_scaffold_representative"] == "TRUE"), "—")
        L.append(f"| `{scaf[:60]}` | {n} | `{rep}` |")
    L.append("")

    # --- 自检 ---
    L.append("## 五、阳性对照自检")
    L.append("")
    L.append("6 个已知临床/参比 GKA。规则跑完它们**必须落进 A/B 层**——"
             "落不进说明规则有问题，不是数据有问题。")
    L.append("")
    L.append("| 分子 | 名称 | phase | direction | 层 | 层内名次 | pAct 中位 | 证据 |")
    L.append("| --- | --- | --- | --- | --- | ---: | ---: | --- |")
    ok = True
    for mid, name in POSITIVE_CONTROLS.items():
        r = next((x for x in rows if x["molecule_chembl_id"] == mid), None)
        if r is None:
            L.append(f"| `{mid}` | {name} | — | **不在输入里** | — | — | — | — |")
            ok = False
            continue
        tier = r["tier"] or "（排除）"
        if tier not in ("A", "B"):
            ok = False
        L.append(f"| `{mid}` | {name} | {r['max_phase'] or '—'} | {r['direction']} | "
                 f"**{tier}** | {r['rank_in_tier'] or '—'} | "
                 f"{r['pactivity_median'] or '—'} | {r['evidence_level']} |")
    L.append("")
    L.append(f"**自检结论：{'通过' if ok else '未通过'}**"
             f"（{'6/6' if ok else str(sum(1 for m in POSITIVE_CONTROLS if next((x for x in rows if x['molecule_chembl_id'] == m), {}).get('tier') in ('A', 'B')))+'/6'}"
             " 落在 A/B 层）。")
    L.append("")
    L.append("> 注意三个 phase 2 药的证据等级是「弱」，参比化合物 Ro-28-1675 是「强」。"
             "证据衡量的是「在 ChEMBL 里被测了多少次」，不是分子有多好——"
             "**这就是证据不进排序键的原因**。")
    L.append("")

    # --- 榜单 ---
    L.append("## 六、A 层完整名单")
    L.append("")
    L.append("按 `pactivity_median` 降序。效能与证据并列展示，不参与排序。")
    L.append("")
    L.append("| # | molecule | 名称 | pAct 中位/最优 | EC50 最小(nM) | 效能 fold/% | "
             "证据 | 骨架簇 | 代表 | 打标 |")
    L.append("| ---: | --- | --- | --- | ---: | --- | --- | ---: | :---: | --- |")
    for r in [x for x in ranked if x["tier"] == "A"]:
        ef = "/".join([r["efficacy_fold_max"] or "—", r["efficacy_pct_max"] or "—"])
        L.append(f"| {r['rank_in_tier']} | `{r['molecule_chembl_id']}` | "
                 f"{(r['molecule_pref_name'] or '—')[:20]} | "
                 f"{r['pactivity_median']} / {r['pactivity_max']} | "
                 f"{r['potency_nm_min'] or '—'} | {ef} | {r['evidence_level']} | "
                 f"{r['scaffold_cluster_size']} | "
                 f"{'✓' if r['is_scaffold_representative'] == 'TRUE' else ''} | "
                 f"{r['flags'] or ''} |")
    L.append("")
    L.append("## 七、B 层前 30")
    L.append("")
    L.append("| # | molecule | pAct 中位 | EC50 最小(nM) | 效能 fold/% | 证据 | 骨架簇 | 代表 |")
    L.append("| ---: | --- | ---: | ---: | --- | --- | ---: | :---: |")
    for r in [x for x in ranked if x["tier"] == "B"][:30]:
        ef = "/".join([r["efficacy_fold_max"] or "—", r["efficacy_pct_max"] or "—"])
        L.append(f"| {r['rank_in_tier']} | `{r['molecule_chembl_id']}` | "
                 f"{r['pactivity_median']} | {r['potency_nm_min'] or '—'} | {ef} | "
                 f"{r['evidence_level']} | {r['scaffold_cluster_size']} | "
                 f"{'✓' if r['is_scaffold_representative'] == 'TRUE' else ''} |")
    L.append("")
    L.append("完整名单（含 C 层、被排除的分子及其理由）见同目录 "
             "`Step1_05_GCK_Activator_Candidates.csv`。")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(
        description="对 Step1_04 的分子做方向判定、排序与候选选择。")
    p.add_argument("--mol-csv", type=Path, default=DEFAULT_MOL)
    p.add_argument("--act-csv", type=Path, default=DEFAULT_ACT)
    p.add_argument("--prev-md", type=Path, default=DEFAULT_PREV_MD,
                   help="Step1_04 报告，用来取 ChEMBL 版本与数据库路径")
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    for f in (args.mol_csv, args.act_csv):
        if not f.is_file():
            print(f"错误：找不到 {f}，请先运行 Step1_04。", file=sys.stderr)
            return 1

    version, db = read_provenance(args.prev_md)
    mols = read_csv(args.mol_csv)
    acts = read_csv(args.act_csv)
    print(f"ChEMBL 版本：{version}")
    print(f"输入：{args.mol_csv.name}（{len(mols):,} 分子）、"
          f"{args.act_csv.name}（{len(acts):,} activity）\n")

    baselines = resolve_baselines(acts)
    nb = Counter(i["baseline"] for i in baselines.values())
    print(f"基线解析：{len(baselines)} 组 "
          + "，".join(f"{k} {n}" for k, n in nb.most_common()))

    rows = build_rows(mols, acts, baselines)
    dc = Counter(r["direction"] for r in rows)
    print("\n方向判定：")
    for k, n in dc.most_common():
        print(f"  {k:<24s} {n:>6,}")

    assign_tiers(rows)
    assign_scaffolds(rows)

    inc = [r for r in rows if r["included"] == "TRUE"]
    print(f"\n入选 {len(inc):,} / {len(rows):,}")
    tc = Counter(r["tier"] for r in inc)
    for k in ("A", "B", "C", "无效力（只有效能/机制读数）"):
        if tc.get(k):
            print(f"  {k:<24s} {tc[k]:>6,}")
    print(f"骨架 {len({r['murcko_scaffold'] for r in rows if r['murcko_scaffold']}):,} 个")

    bad = [m for m in POSITIVE_CONTROLS
           if next((x for x in rows if x["molecule_chembl_id"] == m), {}).get("tier")
           not in ("A", "B")]
    print(f"\n阳性对照自检：{6 - len(bad)}/6 落在 A/B 层"
          + ("" if not bad else f"  ← 未通过：{bad}"))

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = args.outdir / "Step1_05_GCK_Activator_Candidates.csv"
    out_md = args.outdir / "Step1_05_GCK_Activator_Ranking_and_Candidate_Selection.md"
    rows.sort(key=lambda r: (r["rank_overall"] if r["rank_overall"] != "" else 10 ** 9,
                             r["molecule_chembl_id"]))
    write_csv(rows, out_csv)
    write_report(rows, acts, baselines, out_md, args.mol_csv, args.act_csv, version, db)

    print(f"\n候选主表：{out_csv}")
    print(f"报告：    {out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
