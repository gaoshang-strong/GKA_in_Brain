#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_Integrate_Candidates.py — 汇总 Step1 全链路的候选分子与理化性质

把两条**来源不同**的候选合并成一张总表，供 Step2 直接使用：

    Step1_06  活性路径   782 个   有效力值（pActivity）、可排序、方向靠读数推
    Step1_07  药物注释路径  9 个   无活性值、不可排序、方向是 ChEMBL 人工审编的

**两条路径的口径不同，合并后必须能分辨**——`source` 列标明每一行从哪来，
`pactivity_median` 等排序相关字段对药物注释路径的分子为空，**这是事实不是缺失**。

去重键：为什么用 parent_chembl_id 而不是 InChIKey
------------------------------------------------
同一个药会以「游离碱 + 盐」两个 chembl_id 存在，**两者 InChIKey 完全不同**：

    MK-0941       CHEMBL3580737 游离碱 KJSGTWFWVTYPFZ-…
                  CHEMBL4297302 甲磺酸盐 PIDNRTWDGDJKSQ-…   (SMILES 带 .CS(=O)(=O)O)
    Globalagliatin CHEMBL4297399 LY-2608204 QIIVJLHCZUTGSD-…
                  CHEMBL5095182 盐酸盐  FRUQQNDJVRDIRH-…   (SMILES 带 Cl.)

按 InChIKey 去重会把同一个药算成两个。`molecule_hierarchy.parent_molregno`
能正确归并这两对，所以**药物层去重用 `parent_chembl_id`**。
（跨库对齐如与 SureChEMBL 对结构时仍用 InChIKey——两个场景用不同的键。）

产物保留**一行一个 chembl_id**（保持一行一个实体），另给 `dedup_group` 列标注归属，
`is_dedup_representative` 标出每组的代表，不擅自删行。

输出
----
    Step1_GKA_Candidates_with_Properties.csv   整合总表（本目录，Step2 的输入）
    Step1_GKA_Candidates_Summary.md            汇总说明

用法
----
    python3 Step1_Integrate_Candidates.py
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
SRC_06 = (HERE / "Step1_06_GKA_Physicochemical_Property_Extraction"
          / "Step1_06_GKA_Physicochemical_Properties.csv")
SRC_07 = (HERE / "Step1_07_GKA_from_Drug_Annotation"
          / "Step1_07_GKA_from_Drug_Annotation.csv")

# 阳性对照（CLAUDE.md 的 11 个），用于在总表里打标
POSITIVE_CONTROLS = {
    "CHEMBL4297508": "Dorzagliatin",
    "CHEMBL1783734": "Piraglitin",
    "CHEMBL2165615": "Neriglitin",
    "CHEMBL2165620": "PF-04991532",
    "CHEMBL3219124": "AZD-1656",
    "CHEMBL3580737": "MK-0941 (free base)",
    "CHEMBL4297302": "MK-0941 (mesylate)",
    "CHEMBL4297399": "LY-2608204 / Globalagliatin",
    "CHEMBL5095182": "Globalagliatin HCl",
    "CHEMBL5095262": "Cadisegliatin",
    "CHEMBL5072532": "BMS-820132",
    "CHEMBL1096435": "Ro-28-1675 (参比化合物)",
}

# 从 Step1_06 直接继承的列（活性路径独有的排序相关字段 + 理化性质）
FROM_06 = [
    "molecule_pref_name", "max_phase", "first_approval", "molecule_type",
    "structure_type", "canonical_smiles", "standard_inchi_key", "standard_inchi",
    "mw_freebase", "full_mwt", "alogp", "hba", "hbd", "psa", "rtb", "ro3_pass",
    "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
    "full_molformula", "np_likeness_score",
    "n_structural_alerts", "structural_alert_sets",
    "n_ligand_eff", "le_max", "bei_max", "sei_max", "lle_max",
    "priority", "potency_band", "pactivity_median", "potency_nm_min", "direction",
    "murcko_scaffold", "scaffold_cluster_size", "is_scaffold_representative",
    "parent_chembl_id", "is_parent", "n_synonyms", "synonyms",
]

OUT_COLUMNS = [
    # --- 身份 ---
    "molecule_chembl_id", "molecule_pref_name", "max_phase", "first_approval",
    "molecule_type", "structure_type",
    # --- 来源与口径（本表核心）---
    "source", "activity_evidence", "curated_direction", "curated_mechanism",
    "found_by_drug_annotation",
    # --- 去重 ---
    "parent_chembl_id", "dedup_group", "is_dedup_representative", "is_salt",
    # --- 活性路径的排序字段（药物注释路径为空，是事实不是缺失）---
    "priority", "potency_band", "pactivity_median", "potency_nm_min", "direction",
    # --- 结构 ---
    "canonical_smiles", "standard_inchi_key",
    # --- 理化性质 ---
    "mw_freebase", "full_mwt", "alogp", "hba", "hbd", "psa", "rtb",
    "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
    "full_molformula",
    # --- 结构警示与配体效率 ---
    "n_structural_alerts", "structural_alert_sets", "lle_max",
    # --- 骨架 ---
    "murcko_scaffold", "scaffold_cluster_size", "is_scaffold_representative",
    # --- 对照 ---
    "is_positive_control", "control_name",
]


def read_csv(p: Path) -> list:
    if not p.is_file():
        sys.exit(f"错误：找不到 {p}")
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def murcko(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:                                   # noqa: BLE001
        return ""


def fnum(s):
    if s in (None, ""):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def build(rows06: list, rows07: list) -> list:
    by_id = {}

    for r in rows06:
        cid = r["molecule_chembl_id"]
        d = {k: r.get(k, "") for k in FROM_06}
        d["molecule_chembl_id"] = cid
        d["source"] = "activity"
        d["activity_evidence"] = "TRUE"
        d["curated_direction"] = ""
        d["curated_mechanism"] = ""
        d["found_by_drug_annotation"] = ""
        by_id[cid] = d

    for r in rows07:
        cid = r["molecule_chembl_id"]
        if cid in by_id:                       # 两条路径都命中
            d = by_id[cid]
            d["source"] = "activity+drug_annotation"
        else:                                  # 只有药物注释路径 —— 新补回来的
            d = {k: "" for k in FROM_06}
            d["molecule_chembl_id"] = cid
            d["source"] = "drug_annotation"
            d["activity_evidence"] = "FALSE"
            for k in ("molecule_pref_name", "max_phase", "first_approval",
                      "molecule_type", "structure_type", "canonical_smiles",
                      "standard_inchi_key", "mw_freebase", "full_mwt", "alogp",
                      "hba", "hbd", "psa", "rtb", "num_ro5_violations",
                      "aromatic_rings", "heavy_atoms", "qed_weighted",
                      "full_molformula", "parent_chembl_id", "is_parent"):
                d[k] = r.get(k, "")
            d["murcko_scaffold"] = murcko(r.get("canonical_smiles", ""))
            by_id[cid] = d
        d["curated_direction"] = r.get("dm_action_type", "")
        d["curated_mechanism"] = r.get("dm_mechanism_of_action", "")
        d["found_by_drug_annotation"] = r.get("found_by", "")

    # --- 去重分组：药物层用 parent_chembl_id（盐与游离碱的 InChIKey 不同，归并不了）---
    for d in by_id.values():
        parent = d.get("parent_chembl_id") or ""
        d["dedup_group"] = parent if parent else d["molecule_chembl_id"]
        smi = d.get("canonical_smiles") or ""
        d["is_salt"] = "TRUE" if "." in smi else "FALSE"

    groups = defaultdict(list)
    for d in by_id.values():
        groups[d["dedup_group"]].append(d)
    for g, members in groups.items():
        for d in members:
            d["is_dedup_representative"] = "FALSE"
        # 代表优先取：非盐 > 有活性证据 > chembl_id 字典序
        best = sorted(members, key=lambda x: (x["is_salt"] == "TRUE",
                                              x["activity_evidence"] != "TRUE",
                                              x["molecule_chembl_id"]))[0]
        best["is_dedup_representative"] = "TRUE"

    for d in by_id.values():
        cid = d["molecule_chembl_id"]
        d["is_positive_control"] = "TRUE" if cid in POSITIVE_CONTROLS else "FALSE"
        d["control_name"] = POSITIVE_CONTROLS.get(cid, "")

    rows = list(by_id.values())
    rows.sort(key=lambda x: (x["source"] != "activity+drug_annotation",
                             {"P1": 0, "P2": 1, "P3": 2, "P4": 3}.get(x.get("priority"), 9),
                             -(fnum(x.get("pactivity_median")) or 0),
                             x["molecule_chembl_id"]))
    return rows


def write_csv(rows: list, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in OUT_COLUMNS})


def write_report(rows: list, path: Path, s6: Path, s7: Path) -> None:
    L = []
    n = len(rows)
    src = Counter(r["source"] for r in rows)
    L.append("# Step1 候选分子与理化性质整合表")
    L.append("")
    L.append(f"- 运行时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    L.append(f"- 输入：`{s6.parent.name}/{s6.name}`、`{s7.parent.name}/{s7.name}`")
    L.append(f"- 合计 **{n:,}** 个 ChEMBL ID")
    L.append("")
    L.append("> 本表是 **Step1 全链路的最终产物**，Step2 直接读它。"
             "两条来源不同的路径合并在一起，`source` 列标明每一行从哪来。")
    L.append("")

    L.append("## 一、两条来源，口径不同")
    L.append("")
    L.append("| source | 分子数 | 有效力值 | 方向来自 | 能否排序 |")
    L.append("| --- | ---: | :---: | --- | :---: |")
    L.append(f"| `activity` | {src.get('activity', 0):,} | ✅ | 从 assay 读数推（Step1_05） | ✅ |")
    L.append(f"| `activity+drug_annotation` | {src.get('activity+drug_annotation', 0):,} | ✅ | "
             "**两者都有**，人工审编可复核推断结果 | ✅ |")
    L.append(f"| `drug_annotation` | {src.get('drug_annotation', 0):,} | ❌ | "
             "**ChEMBL 人工审编**（`action_type`） | ❌ |")
    L.append("")
    L.append("**`drug_annotation` 那批没有 `pactivity_median` / `priority` / `potency_band`——"
             "这是事实不是缺失。** 它们在 ChEMBL 里一条活性记录都没有（或没有打在 GCK 上），"
             "无法参与效力分档与排序。但它们的方向是**人工审编**的，比推断的可靠。")
    L.append("")
    only_da = [r for r in rows if r["source"] == "drug_annotation"]
    if only_da:
        L.append("只靠药物注释路径进来的（Step1_07 补回的）：")
        L.append("")
        L.append("| 分子 | 名称 | phase | 命中路径 | 审编方向 |")
        L.append("| --- | --- | ---: | --- | --- |")
        for r in sorted(only_da, key=lambda x: -(fnum(x.get("max_phase")) or 0)):
            L.append(f"| `{r['molecule_chembl_id']}` | {r.get('molecule_pref_name') or '—'} | "
                     f"{r.get('max_phase') or '—'} | {r.get('found_by_drug_annotation') or '—'} | "
                     f"{r.get('curated_direction') or '—'} |")
        L.append("")

    L.append("## 二、⚠ 去重键：`parent_chembl_id`，不是 InChIKey")
    L.append("")
    groups = defaultdict(list)
    for r in rows:
        groups[r["dedup_group"]].append(r)
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    L.append(f"**{n:,} 个 chembl_id → {len(groups):,} 个去重组**，"
             f"其中 {len(multi)} 组含多个 ID。")
    L.append("")
    if multi:
        L.append("| 去重组 | 成员 | 说明 |")
        L.append("| --- | --- | --- |")
        for g, members in sorted(multi.items(), key=lambda kv: -len(kv[1]))[:10]:
            ids = "、".join(
                f"`{m['molecule_chembl_id']}`" + ("（盐）" if m["is_salt"] == "TRUE" else "")
                for m in members)
            names = {m.get("control_name") for m in members if m.get("control_name")}
            L.append(f"| `{g}` | {ids} | {'、'.join(names) if names else ''} |")
        L.append("")
    L.append("**同一个药以「游离碱 + 盐」两个 chembl_id 存在时，两者 InChIKey 完全不同**——")
    L.append("")
    L.append("```")
    L.append("MK-0941        CHEMBL3580737 游离碱   KJSGTWFWVTYPFZ-AWEZNQCLSA-N")
    L.append("               CHEMBL4297302 甲磺酸盐 PIDNRTWDGDJKSQ-UQKRIMTDSA-N  ← SMILES 带 .CS(=O)(=O)O")
    L.append("Globalagliatin CHEMBL4297399 游离碱   QIIVJLHCZUTGSD-CUBQBAPOSA-N")
    L.append("               CHEMBL5095182 盐酸盐   FRUQQNDJVRDIRH-JOFLZTHPSA-N  ← SMILES 带 Cl.")
    L.append("```")
    L.append("")
    L.append("按 InChIKey 去重会把同一个药算成两个。"
             "`molecule_hierarchy.parent_molregno` 能正确归并，所以**药物层去重用 "
             "`parent_chembl_id`**；跨库对齐结构（如与 SureChEMBL）才用 InChIKey，"
             "且要先归到母体。**两个场景用不同的键。**")
    L.append("")
    L.append("产物保留一行一个 chembl_id，`is_dedup_representative` 标出每组代表"
             "（优先非盐、有活性证据），**不擅自删行**。")
    L.append("")

    L.append("## 三、理化性质覆盖")
    L.append("")
    L.append("| 字段 | 有值 | 缺失 |")
    L.append("| --- | ---: | ---: |")
    for f in ("canonical_smiles", "standard_inchi_key", "mw_freebase", "alogp",
              "psa", "hbd", "hba", "rtb", "qed_weighted", "murcko_scaffold"):
        have = sum(1 for r in rows if r.get(f))
        L.append(f"| `{f}` | {have:,} | {n - have:,} |")
    L.append("")
    for label, f, win, ok in (("MW", "mw_freebase", "≤ 450", lambda v: v <= 450),
                              ("ALogP", "alogp", "1–5", lambda v: 1 <= v <= 5),
                              ("PSA", "psa", "≤ 90", lambda v: v <= 90)):
        vals = [fnum(r.get(f)) for r in rows]
        vals = [v for v in vals if v is not None]
        if vals:
            L.append(f"- **{label}**：中位 {statistics.median(vals):.2f}，"
                     f"落在常用 CNS 窗口（{win}）的 {sum(1 for v in vals if ok(v)):,} / {len(vals):,}")
    L.append("")
    L.append("> 性质窗口只作描述，**本表不做任何脑暴露筛选**——留给 Step2。"
             "另注意 ChEMBL 37 没有 `cx_logd` / `molecular_species`，`alogp` 不是 logD。")
    L.append("")

    L.append("## 四、阳性对照在表内的位置")
    L.append("")
    ctl = [r for r in rows if r["is_positive_control"] == "TRUE"]
    L.append(f"**{len(ctl)} / {len(POSITIVE_CONTROLS)}** 个对照在表内。")
    L.append("")
    L.append("| 分子 | 名称 | phase | source | 效力档 | pAct 中位 |")
    L.append("| --- | --- | ---: | --- | --- | ---: |")
    for r in sorted(ctl, key=lambda x: -(fnum(x.get("max_phase")) or 0)):
        L.append(f"| `{r['molecule_chembl_id']}` | {r.get('control_name')} | "
                 f"{r.get('max_phase') or '—'} | {r['source']} | "
                 f"{r.get('potency_band') or '—'} | {r.get('pactivity_median') or '—'} |")
    miss = [k for k in POSITIVE_CONTROLS if not any(
        r["molecule_chembl_id"] == k for r in rows)]
    if miss:
        L.append("")
        L.append("⚠ 不在表内的对照：" + "、".join(f"`{m}`" for m in miss))
    L.append("")

    L.append("## 五、下游怎么用")
    L.append("")
    L.append("| 想做的事 | 该怎么筛 |")
    L.append("| --- | --- |")
    L.append("| 按效力挑候选 | `source` 含 `activity` 且 `priority` / `potency_band` 非空 |")
    L.append("| 要方向最可靠的 | `curated_direction = 'ACTIVATOR'`（ChEMBL 人工审编） |")
    L.append("| 结构去冗余 | `is_scaffold_representative = TRUE` |")
    L.append("| 药物层计数 | 按 `dedup_group` 或 `is_dedup_representative = TRUE` |")
    L.append("| 与 SureChEMBL 对齐 | `standard_inchi_key`，**先归母体再取 key** |")
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="整合 Step1 全链路的候选分子与理化性质。")
    ap.add_argument("--step1-06", type=Path, default=SRC_06)
    ap.add_argument("--step1-07", type=Path, default=SRC_07)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    r6, r7 = read_csv(args.step1_06), read_csv(args.step1_07)
    print(f"Step1_06 活性路径：    {len(r6):,} 个")
    print(f"Step1_07 药物注释路径：{len(r7):,} 个")

    rows = build(r6, r7)
    src = Counter(r["source"] for r in rows)
    print(f"\n整合后：{len(rows):,} 个 ChEMBL ID")
    for k in ("activity", "activity+drug_annotation", "drug_annotation"):
        if src.get(k):
            print(f"  {k:<26s} {src[k]:>5,}")
    groups = {r["dedup_group"] for r in rows}
    print(f"去重组（parent_chembl_id）：{len(groups):,}")
    ctl = sum(1 for r in rows if r["is_positive_control"] == "TRUE")
    print(f"阳性对照在表内：{ctl}/{len(POSITIVE_CONTROLS)}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = args.outdir / "Step1_GKA_Candidates_with_Properties.csv"
    out_md = args.outdir / "Step1_GKA_Candidates_Summary.md"
    write_csv(rows, out_csv)
    write_report(rows, out_md, args.step1_06, args.step1_07)
    print(f"\n整合总表：{out_csv}")
    print(f"汇总说明：{out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
