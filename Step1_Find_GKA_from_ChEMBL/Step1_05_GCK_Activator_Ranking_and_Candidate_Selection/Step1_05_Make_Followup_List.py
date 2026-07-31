#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_05_Make_Followup_List
===========================

从 Step1_05 的候选主表切出**进入后续实验的清单**，给每个分子一个 `priority`。

这里只做「切」，不做任何新的判断——所有列都来自
`Step1_05_GCK_Activator_Candidates.csv`，切法写死在下面的 `PRIORITY_RULES` 里。

切的口径
--------
入选池 = 通过 Step1_05 两道门（方向不为无激活、效力非仅删失）**且**有可定量效力值
        **且** pActivity 中位 ≥ 6.0

6.0 这条线由阳性对照标定：6 个已知临床/参比 GKA 全部在其之上，最低的 Piraglitin
是 6.145。不是按分布拍的。

| priority | 条件 | 用途 |
|---|---|---|
| P1 | pAct ≥ 7.0 + 骨架代表 + `direction == activation` | 首批实验，方向是实测的 |
| P2 | pAct ≥ 7.0 + 骨架代表（方向靠 EC50 推定） | 次批，需先确认方向 |
| P3 | 6.0 ≤ pAct < 7.0 + 骨架代表 | 结构起点储备 |
| P4 | pAct ≥ 6.0，非骨架代表 | 同系物，备查不单独下单 |

**骨架代表**是每个 Murcko 簇里效力最优的那一个。按效力取 top-N 会拿到同一篇 SAR
论文的一串同系物——1,333 个分子只有 521 个骨架，最大一簇 46 个。

P1/P2 的区别是**方向证据的来源**，不是强弱：`activation` 有实测的效能读数支持，
`activation_by_potency` 只是「EC50 出自激活读数的 assay」这一推定。

**阳性对照无论落在哪一档都单独标 `is_positive_control`**，实验时应带上作基准。

不在这张表里的口径
------------------
效能、证据强度、物化性质、脑暴露**都没有参与切分**。脑暴露留给 Step2——
这张表回答的是「谁是好的 GCK 激活剂」，不是「谁能进脑」。

用法
----
    python3 Step1_05_Make_Followup_List.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / "Step1_05_GCK_Activator_Candidates.csv"

POOL_MIN_PACT = 6.0    # 由阳性对照标定（Piraglitin 6.145 是 6 个里最低的）
P1_MIN_PACT = 7.0

POSITIVE_CONTROLS = {
    "CHEMBL1096435": "Ro-28-1675",
    "CHEMBL1783734": "PIRAGLIATIN",
    "CHEMBL2165615": "NERIGLIATIN",
    "CHEMBL2165620": "PF-04991532",
    "CHEMBL3219124": "AZD-1656",
    "CHEMBL3580737": "MK-0941",
}

PRIORITY_NOTE = {
    "P1": "首批：效力 ≥7.0，骨架代表，方向有实测效能读数支持",
    "P2": "次批：效力 ≥7.0，骨架代表，方向靠 EC50 推定，需先确认方向",
    "P3": "储备：效力 6.0–7.0，骨架代表，作结构起点",
    "P4": "备查：效力 ≥6.0 但非骨架代表，是已入选骨架的同系物",
}

OUT_COLUMNS = [
    "priority", "priority_note", "is_positive_control", "control_name",
    "molecule_chembl_id", "molecule_pref_name", "max_phase", "canonical_smiles",
    "rank_overall", "potency_band",
    "pactivity_median", "pactivity_max", "potency_nm_min",
    "direction", "has_efficacy_corroboration",
    "efficacy_fold_max", "efficacy_pct_max",
    "evidence_level", "n_assays", "n_docs", "flags",
    "murcko_scaffold", "scaffold_cluster_size", "is_scaffold_representative",
]


def fnum(s):
    if s is None or s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def priority_of(row) -> str:
    pa = fnum(row["pactivity_median"])
    is_rep = row["is_scaffold_representative"] == "TRUE"
    if pa >= P1_MIN_PACT and is_rep:
        return "P1" if row["direction"] == "activation" else "P2"
    if is_rep:
        return "P3"
    return "P4"


def main() -> int:
    p = argparse.ArgumentParser(description="从 Step1_05 候选主表切出后续实验清单。")
    p.add_argument("--in-csv", type=Path, default=DEFAULT_IN)
    p.add_argument("--outdir", type=Path, default=HERE)
    args = p.parse_args()

    if not args.in_csv.is_file():
        print(f"错误：找不到 {args.in_csv}，请先运行 Step1_05。", file=sys.stderr)
        return 1

    with args.in_csv.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    pool = [r for r in rows
            if r["included"] == "TRUE"
            and fnum(r["pactivity_median"]) is not None
            and fnum(r["pactivity_median"]) >= POOL_MIN_PACT]

    out = []
    for r in pool:
        pri = priority_of(r)
        mid = r["molecule_chembl_id"]
        rec = {k: r.get(k, "") for k in OUT_COLUMNS}
        rec["priority"] = pri
        rec["priority_note"] = PRIORITY_NOTE[pri]
        rec["is_positive_control"] = "TRUE" if mid in POSITIVE_CONTROLS else "FALSE"
        rec["control_name"] = POSITIVE_CONTROLS.get(mid, "")
        out.append(rec)

    out.sort(key=lambda r: (r["priority"], int(r["rank_overall"] or 10 ** 9)))

    args.outdir.mkdir(parents=True, exist_ok=True)
    dest = args.outdir / "Step1_05_Followup_Candidates.csv"
    with dest.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in out:
            w.writerow(r)

    print(f"输入：{args.in_csv.name}（{len(rows):,} 行）")
    print(f"入选池：pActivity ≥ {POOL_MIN_PACT} 且有可定量效力 → {len(out):,} 个\n")
    pc = Counter(r["priority"] for r in out)
    for k in ("P1", "P2", "P3", "P4"):
        if pc.get(k):
            n_ctl = sum(1 for r in out
                        if r["priority"] == k and r["is_positive_control"] == "TRUE")
            print(f"  {k}  {pc[k]:>4,}  {PRIORITY_NOTE[k]}"
                  + (f"  [含 {n_ctl} 个阳性对照]" if n_ctl else ""))
    missing = [n for m, n in POSITIVE_CONTROLS.items()
               if not any(r["molecule_chembl_id"] == m for r in out)]
    print(f"\n阳性对照：{len(POSITIVE_CONTROLS) - len(missing)}/{len(POSITIVE_CONTROLS)} 在表内"
          + (f"  ← 缺：{missing}" if missing else ""))
    print(f"骨架：{len({r['murcko_scaffold'] for r in out if r['murcko_scaffold']}):,} 个")
    print(f"\n后续实验清单：{dest}")
    print(f"生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
