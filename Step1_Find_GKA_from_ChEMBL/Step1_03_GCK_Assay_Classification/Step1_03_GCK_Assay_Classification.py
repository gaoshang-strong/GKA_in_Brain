#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_03_GCK_Assay_Classification（第 1 阶段：规则分类）
======================================================

判断每个 GCK assay 实际测量的内容，分为：

    GCK 激活 / GCK 抑制 / GCK 结合 / GCK–GKRP 相互作用 / 细胞或表型效应 / 无法判断

本脚本只做**规则分类**，对应 README 的第 1 步。规则命中不明确的记录一律标记
`review_required = True`，交给第 2 步（LLM）与第 4 步（人工）处理。
脚本不会替模糊记录强行选一个标签。

判定依据
--------
    target_chembl_id + assay_description + assay_type + bao_format
    + 该 assay 下 activities 的 standard_type / units / relation / activity_comment

设计要点：从真实数据里总结出的几个坑
------------------------------------
1. `bao_format` 不可靠。本数据中多个纯酶实验被标为 "tissue-based format"
   （如 CHEMBL5048855，描述是 E. coli 重组表达的酶），因此**不用它判定细胞效应**，
   只作为参考信息记录。
2. "expressed in Sf9/CHO cells" 是表达宿主，不是细胞实验。需要把
   "expressed in ... cells" 从细胞实验的判据里排除。
3. "affinity" 不能作为结合的关键词："Activation of ... assessed as enzyme
   affinity for glucose"（CHEMBL4427350）测的是激活。只认 "binding affinity"。
4. "binding" 同理："Inhibition of ... by filter-binding assay"（CHEMBL3829755）
   测的是抑制，binding 只是方法名。
5. GKRP 必须先于 inhibition 判断："Inhibition of GK interaction with GKRP"
   实为 GKRP 相互作用实验，虽然它挂在 CHEMBL3820 上。
   但 "in absence of human GKRP"（CHEMBL3751340）不是 —— 这类进人工复核。

输入 / 输出
-----------
输入：`../Step1_02_GCK_Assay_Mapping/Step1_02_GCK_Assay_Mapping.csv`
输出：在 assay 主表基础上增加 4 列（README 要求）+ 1 列方法标记
      `assay_category` / `classification_confidence` / `classification_reason`
      / `review_required` / `classification_method`

    Step1_03_GCK_Assay_Classification.csv
    Step1_03_GCK_Assay_Classification.md
    Step1_03_pending_llm.csv        需要第 2 步 LLM 处理的子集

注意：ChEMBL 原生也有一个 `assay_category` 字段（本数据中 228 行全为空）。
为避免与 README 要求的新列同名，原生列在输出中改名为 `assay_category_chembl`。

用法
----
    python3 Step1_03_GCK_Assay_Classification.py
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE.parent / "Step1_02_GCK_Assay_Mapping" / "Step1_02_GCK_Assay_Mapping.csv"

# ---------------------------------------------------------------------------
# 类别常量
# ---------------------------------------------------------------------------
ACTIVATION = "GCK 激活"
INHIBITION = "GCK 抑制"
BINDING = "GCK 结合"
GKRP = "GCK–GKRP 相互作用"
CELLULAR = "细胞或表型效应"
UNKNOWN = "无法判断"

CATEGORY_ORDER = [ACTIVATION, INHIBITION, BINDING, GKRP, CELLULAR, UNKNOWN]

PPI_TARGET = "CHEMBL3885579"      # Glucokinase/Glucokinase regulatory protein

# ---------------------------------------------------------------------------
# 正则
# ---------------------------------------------------------------------------
RE_GKRP = re.compile(r"\bGKRP\b|glucokinase regulatory protein", re.I)
RE_GKRP_INTERACTION = re.compile(
    r"(interaction|complex|dissociat|disrupt|binding)\W{0,40}GKRP"
    r"|GKRP\W{0,40}(interaction|complex|dissociat|disrupt|binding)", re.I)
# "in absence/presence of GKRP" —— 是实验条件，不是在测 GK-GKRP 相互作用
RE_GKRP_CONDITION = re.compile(r"\b(in|with|without)\s+(the\s+)?(absence|presence)\s+of\b"
                               r"[^.;]{0,40}GKRP", re.I)

RE_ACTIVATION = re.compile(r"\bactivat(?:ion|or|ors|ing|es|ed|e)\b", re.I)
RE_INDUCTION = re.compile(r"\binduction of\b[^.;]{0,80}\bactivity\b", re.I)
RE_INHIBITION = re.compile(r"\binhibit(?:ion|or|ors|ory|ing|s|ed)\b", re.I)

# 结合：只认明确的措辞，避免 "filter-binding assay"、"enzyme affinity for glucose" 误伤
RE_BINDING = re.compile(r"\bbinding affinit(?:y|ies)\b|\bdisplacement of\b"
                        r"|\bsurface plasmon resonance\b|\bSPR\b", re.I)
RE_BINDING_FALSE = re.compile(r"filter[- ]binding|binding assay\b", re.I)

# 细胞/表型：先剔除 "expressed in ... cells" 这类表达宿主的说法
RE_EXPRESSION_HOST = re.compile(r"expressed in[^.;]{0,60}?\bcells?\b", re.I)
RE_CELLULAR = re.compile(
    r"\bin\s+[\w\-/]+\s+cells\b|\bhepatocytes?\b|\bislets?\b|\bin vivo\b"
    r"|\bmice\b|\brats\b|\bmouse\b|insulin secretion|glucose uptake"
    r"|glucose lowering|translocation", re.I)
# 明确的表型读数（不是在测 GCK 酶活本身）
RE_PHENOTYPE = re.compile(r"insulin secretion|glucose uptake|glucose lowering"
                          r"|translocation|blood glucose", re.I)

# standard_type → 方向倾向
ST_ACTIVATION = {"EC50", "Emax", "%max", "S0.5", "S50", "AC50", "FC", "Activation",
                 "Ratio EC50", "Activity ratio", "Fold change"}
ST_INHIBITION = {"IC50", "Inhibition", "% inhibition", "INHIBITION", "Ki"}
ST_BINDING = {"Kd", "KD", "k_off", "kon", "k_on"}
# 方向中性的酶动力学参数：单独出现时不足以定方向
ST_NEUTRAL = {"Km", "Vmax", "K", "Activity", "Ratio", "TIME", "T1/2", "n", "Hill coefficient"}


def parse_standard_types(s: str) -> Counter:
    """把 Step1_02 写出的 'EC50:51; Emax:3' 还原成 Counter。"""
    c: Counter = Counter()
    for part in filter(None, (s or "").split("; ")):
        if ":" in part:
            k, v = part.rsplit(":", 1)
            try:
                c[k] = int(v)
            except ValueError:
                pass
    return c


def snippet(text: str, m: re.Match, width: int = 60) -> str:
    """截取命中处前后的原文，作为判定证据。"""
    a = max(0, m.start() - width // 2)
    b = min(len(text), m.end() + width)
    return ("…" if a > 0 else "") + text[a:b].strip() + ("…" if b < len(text) else "")


def st_direction(st: Counter) -> tuple[set, str]:
    """由 standard_type 分布推断方向倾向。返回 (方向集合, 描述串)。"""
    dirs, notes = set(), []
    for k, v in st.items():
        if k in ST_ACTIVATION:
            dirs.add(ACTIVATION); notes.append(f"{k}:{v}→激活")
        elif k in ST_INHIBITION:
            dirs.add(INHIBITION); notes.append(f"{k}:{v}→抑制")
        elif k in ST_BINDING:
            dirs.add(BINDING); notes.append(f"{k}:{v}→结合")
    return dirs, "、".join(notes)


def classify(row: dict) -> dict:
    """对单个 assay 做规则分类。返回 4 个新字段 + 方法标记。

    规则按优先级顺序执行，先命中先定类；每条规则都要留下可复核的证据。
    """
    desc = row.get("assay_description") or ""
    st = parse_standard_types(row.get("standard_types", ""))
    st_dirs, st_note = st_direction(st)
    reasons: list[str] = []
    conflict = False

    # 表达宿主的 "expressed in ... cells" 不算细胞实验，判定细胞语境前先剔除
    desc_wo_host = RE_EXPRESSION_HOST.sub(" ", desc)
    m_cell = RE_CELLULAR.search(desc_wo_host)
    m_pheno = RE_PHENOTYPE.search(desc_wo_host)
    cell_ctx = bool(m_cell)
    if cell_ctx:
        reasons.append(f"检测到细胞/体内语境：「{snippet(desc_wo_host, m_cell, 40)}」")

    category, conf = None, None

    # --- R1 GKRP 相互作用（必须先于 inhibition 判断）---
    if row.get("target_chembl_id") == PPI_TARGET:
        category, conf = GKRP, "high"
        reasons.insert(0, f"靶点为 PPI 靶点 {PPI_TARGET}（GCK/GKRP），直接判定")
    elif RE_GKRP.search(desc):
        if RE_GKRP_CONDITION.search(desc) and not RE_GKRP_INTERACTION.search(desc):
            # "in absence of GKRP" —— GKRP 是实验条件而非被测对象
            conflict = True
            reasons.insert(0, "描述出现 GKRP，但属于「absence/presence of GKRP」实验条件表述，"
                              "并非在测 GK–GKRP 相互作用，需人工确认")
        elif RE_GKRP_INTERACTION.search(desc):
            m = RE_GKRP_INTERACTION.search(desc)
            category, conf = GKRP, "high"
            reasons.insert(0, f"描述表明测的是 GK–GKRP 相互作用：「{snippet(desc, m)}」")
        else:
            conflict = True
            reasons.insert(0, "描述提及 GKRP 但未构成明确的相互作用表述，需人工确认")

    # --- R2 结合 ---
    if category is None:
        m = RE_BINDING.search(desc)
        if m and not RE_BINDING_FALSE.search(desc):
            category, conf = BINDING, "high"
            reasons.insert(0, f"描述命中结合类措辞：「{snippet(desc, m)}」")
        elif m and RE_BINDING_FALSE.search(desc):
            reasons.append("出现 binding 字样但属于方法名（filter-binding / binding assay），未据此判定为结合")

    # --- R3 激活 ---
    if category is None:
        m = RE_ACTIVATION.search(desc) or RE_INDUCTION.search(desc)
        if m:
            category, conf = ACTIVATION, "high"
            reasons.insert(0, f"描述命中激活类措辞：「{snippet(desc, m)}」")

    # --- R4 抑制 ---
    if category is None:
        m = RE_INHIBITION.search(desc)
        if m:
            category, conf = INHIBITION, "high"
            reasons.insert(0, f"描述命中抑制类措辞：「{snippet(desc, m)}」")

    # --- R5 表型读数：读的是细胞表型而非 GCK 本身 ---
    if m_pheno:
        if category is None:
            category, conf = CELLULAR, "medium"
            reasons.insert(0, f"描述测量的是细胞/整体表型：「{snippet(desc_wo_host, m_pheno)}」")
        else:
            # 既有方向措辞又有表型读数 —— 复合实验，标记复核
            conflict = True
            reasons.append(f"同时含方向措辞与表型读数「{snippet(desc_wo_host, m_pheno, 40)}」，"
                           "属复合实验，需人工确认归类")

    # --- R6 描述无信号时，退回用 standard_type 推断 ---
    if category is None:
        if len(st_dirs) == 1:
            category = next(iter(st_dirs))
            conf = "medium"
            reasons.insert(0, f"描述无方向措辞，依据实测 standard_type 推断（{st_note}）")
        elif len(st_dirs) > 1:
            category, conf = UNKNOWN, "low"
            conflict = True
            reasons.insert(0, f"描述无方向措辞，且 standard_type 指向多个方向（{st_note}）")
        else:
            category, conf = UNKNOWN, "low"
            reasons.insert(0, "描述无方向措辞，standard_type 也全为方向中性"
                              f"（{'、'.join(f'{k}:{v}' for k, v in st.most_common()) or '无活性数据'}）")

    # --- R7 描述与 standard_type 冲突检测 ---
    if conf == "high" and st_dirs and category in (ACTIVATION, INHIBITION, BINDING):
        if category not in st_dirs:
            conflict = True
            conf = "low"
            reasons.append(f"⚠ 描述判定为「{category}」，但实测 standard_type 指向其他方向（{st_note}）")
        elif len(st_dirs) > 1:
            conf = "medium"
            reasons.append(f"standard_type 同时指向多个方向（{st_note}），已按描述判定")
        else:
            reasons.append(f"standard_type 与描述一致（{st_note}）")

    if conflict and conf == "high":
        conf = "medium"

    review = conf != "high" or category == UNKNOWN or conflict

    return {
        "assay_category": category,
        "classification_confidence": conf,
        "classification_reason": "；".join(reasons),
        "review_required": "TRUE" if review else "FALSE",
        "classification_method": "rule",
    }


def main() -> int:
    p = argparse.ArgumentParser(description="规则法对 GCK assay 分类（第 1 阶段）。")
    p.add_argument("--in-csv", type=Path, default=DEFAULT_IN,
                   help=f"Step1_02 的输出，默认 {DEFAULT_IN}")
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    if not args.in_csv.is_file():
        print(f"错误：找不到输入 {args.in_csv}，请先运行 Step1_02。", file=sys.stderr)
        return 1

    with args.in_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        base_cols = list(reader.fieldnames or [])
        rows = list(reader)

    # ChEMBL 原生 assay_category（本数据全空）让位给 README 要求的新列名
    for r in rows:
        r["assay_category_chembl"] = r.pop("assay_category", "")
    base_cols = ["assay_category_chembl" if c == "assay_category" else c for c in base_cols]

    for r in rows:
        r.update(classify(r))

    new_cols = ["assay_category", "classification_confidence", "classification_reason",
                "review_required", "classification_method"]
    out_cols = base_cols + new_cols

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = args.outdir / "Step1_03_GCK_Assay_Classification.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in out_cols})

    # 需要 LLM / 人工处理的子集，单独出一份，字段裁剪到判断所需
    pending = [r for r in rows if r["review_required"] == "TRUE"]
    pend_cols = ["assay_chembl_id", "target_chembl_id", "assay_description",
                 "assay_type", "bao_label", "standard_types", "activity_comments",
                 "assay_category", "classification_confidence", "classification_reason"]
    pend_csv = args.outdir / "Step1_03_pending_llm.csv"
    with pend_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=pend_cols)
        w.writeheader()
        for r in pending:
            w.writerow({k: r.get(k, "") for k in pend_cols})

    write_report(rows, pending, args.outdir / "Step1_03_GCK_Assay_Classification.md",
                 args.in_csv)

    # 控制台摘要
    print(f"输入：{args.in_csv}")
    print(f"assay 总数：{len(rows)}\n")
    cat = Counter(r["assay_category"] for r in rows)
    for c in CATEGORY_ORDER:
        if cat.get(c):
            act = sum(int(r["n_activities"]) for r in rows if r["assay_category"] == c)
            print(f"  {c:<18s} {cat[c]:>4d} 个 assay，{act:>6,} 条活性")
    print(f"\n置信度：{dict(Counter(r['classification_confidence'] for r in rows))}")
    print(f"需复核（review_required=TRUE）：{len(pending)} / {len(rows)}")
    print(f"\n分类结果：{out_csv}")
    print(f"待复核子集：{pend_csv}")
    return 0


def write_report(rows: list[dict], pending: list[dict], path: Path, src: Path) -> None:
    L: list[str] = []
    cat = Counter(r["assay_category"] for r in rows)
    tot_act = sum(int(r["n_activities"]) for r in rows)

    L.append("# Step1_03 GCK assay 分类结果（第 1 阶段：规则）")
    L.append("")
    L.append(f"- 输入：`{src.name}`")
    L.append(f"- 运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- assay 总数：**{len(rows)}**，活性合计 **{tot_act:,}**")
    L.append(f"- 需要进入第 2 步（LLM）/ 第 4 步（人工）的：**{len(pending)}**")
    L.append("")
    L.append("> 本阶段只用规则分类。规则命中不明确的一律标 `review_required = TRUE`，"
             "**不强行给标签**。每条判定都附证据原文，可逐条复核。")
    L.append("")

    L.append("## 分类分布")
    L.append("")
    L.append("| 类别 | assay 数 | 活性数 | 其中需复核 |")
    L.append("| --- | ---: | ---: | ---: |")
    for c in CATEGORY_ORDER:
        if not cat.get(c):
            continue
        sub = [r for r in rows if r["assay_category"] == c]
        L.append(f"| {c} | {len(sub)} | {sum(int(r['n_activities']) for r in sub):,} | "
                 f"{sum(1 for r in sub if r['review_required'] == 'TRUE')} |")
    L.append("")

    L.append("## 置信度分布")
    L.append("")
    L.append("| 置信度 | assay 数 | 含义 |")
    L.append("| --- | ---: | --- |")
    notes = {"high": "描述中有明确措辞，且与实测 standard_type 不冲突",
             "medium": "描述无措辞靠 standard_type 推断，或存在复合/多方向情形",
             "low": "信号缺失或相互矛盾"}
    for lv in ("high", "medium", "low"):
        n = sum(1 for r in rows if r["classification_confidence"] == lv)
        if n:
            L.append(f"| {lv} | {n} | {notes[lv]} |")
    L.append("")

    if pending:
        L.append("## 待第 2 步处理的 assay")
        L.append("")
        L.append("下列记录规则无法明确判定，已导出到 `Step1_03_pending_llm.csv`。"
                 "按 README，第 2 步交由 LLM 判断，且**必须输出证据句与置信度**；"
                 "仍不能确定的进入人工审核。")
        L.append("")
        L.append("| assay_chembl_id | 规则暂定 | 置信度 | 规则给出的理由 | 描述 |")
        L.append("| --- | --- | --- | --- | --- |")
        for r in pending:
            reason = r["classification_reason"].replace("|", "\\|")
            desc = (r["assay_description"] or "")[:110].replace("|", "\\|")
            L.append(f"| `{r['assay_chembl_id']}` | {r['assay_category']} | "
                     f"{r['classification_confidence']} | {reason} | {desc} |")
        L.append("")

    L.append("## 全部分类明细")
    L.append("")
    L.append("| assay_chembl_id | 类别 | 置信度 | 复核 | 活性数 | 描述 |")
    L.append("| --- | --- | --- | --- | ---: | --- |")
    for r in sorted(rows, key=lambda x: (CATEGORY_ORDER.index(x["assay_category"]),
                                         -int(x["n_activities"]))):
        desc = (r["assay_description"] or "")[:110].replace("|", "\\|")
        L.append(f"| `{r['assay_chembl_id']}` | {r['assay_category']} | "
                 f"{r['classification_confidence']} | {r['review_required']} | "
                 f"{int(r['n_activities']):,} | {desc} |")
    L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
