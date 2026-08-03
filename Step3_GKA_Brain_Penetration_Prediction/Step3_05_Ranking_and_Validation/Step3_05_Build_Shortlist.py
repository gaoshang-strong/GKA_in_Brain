#!/usr/bin/env python
"""
Step3_05：把图 1–3 里「跑到好的那一侧」的 GKA 候选单独摘出来

⚠⚠ **这不是 Step3_05 的候选名单，是一张待讨论的观察清单。**
两条规则用的都是**文献习惯值**（`>0.5`、`MW ≤ 400`、`TPSA ≤ 90`），
**没有在本项目的入脑对照集上标定过**——按 CLAUDE.md 的方法论，
阈值必须用对照标定，这件事属 Step3_05 的正题，还没做。
所以本表的用途是「把要讨论的分子摆上桌」，不是「筛出了谁」。

两条规则（都写进 `rule` 列，一个分子可同时命中）：

  A_both_tools  两个工具都判能入脑：SwissADME `BBB permeant = Yes`
                且 ADMETlab `BBB > 0.5`。对应回答问题 1。
  B_cns_space   落在图 3 那个「CNS 常见区」且模型也说能进：
                MW ≤ 400 且 TPSA ≤ 90 且 ADMETlab `BBB > 0.5`
                且 SwissADME 判**非** P-gp 底物。对应回答问题 2。
  C_gka_control 12 个 **GKA 身份对照**（已知的临床/参比 GKA），**无条件全部写入**。

⚠⚠ C 组是**参照物，不是候选**，而且回答的是另一个问题：
它们证明的是「这个分子是不是 GKA」，**不是「能不能进脑」**——
本项目没有它们任何一个的实测脑暴露数据，所以它们
**既不能当入脑阳性对照，也不能当入脑阴性对照**（CLAUDE.md 明确记过）。
放进本表只有一个用途：**看看已知的 GKA 药在这套预测里落在什么位置**，
好让上面 A/B 两组的数值有个熟悉的参照。**不要拿 C 组去标定任何阈值。**

**A/B 两组不做人工挑选**：命中规则的全部写出来，不按谁"看着更好"取舍。
列的含义见 `Step3_05_Candidate_Shortlist.md`。

用法：
  micromamba run -n GKA_in_Brain python Step3_05_Build_Shortlist.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = (HERE.parent / "Step3_04_Result_Integration"
       / "Step3_04_Integrated_Brain_Penetration_Results.csv")
OUT = HERE / "Step3_05_Candidate_Shortlist.csv"

# 只保留每个工具最要紧的几列，其余留在 Step3_04 那张 261 列的总表里
COLUMNS = [
    # —— 身份 ——
    "mol_id", "source_id", "compound_name", "gka_max_phase",
    # —— 命中了哪条规则；C 组是参照物不是候选 ——
    "rule", "gka_is_positive_control", "gka_control_name",
    # —— SwissADME：项目口径下只有这两列真正有用 ——
    "swissadme_bbb_permeant", "swissadme_pgp_substrate",
    # —— ADMETlab：入脑概率 + 外排概率 + CNS MPO 要用的两个量 ——
    "admetlab_bbb", "admetlab_pgp_sub", "admetlab_logd", "admetlab_pka_basic",
    # —— CNS MPO：总分、换口径的波动、以及唯一两项拉低分数的 T0 ——
    "cnsmpo_score", "cnsmpo_score_variant_spread",
    "cnsmpo_t0_mw", "cnsmpo_t0_tpsa",
    # —— RDKit 本地算的结构性质 ——
    "mw", "tpsa", "clogp", "hbd",
    # —— GKA 这一侧：效力与骨架，入脑之外的另外两个维度 ——
    "gka_priority", "gka_potency_band", "gka_pactivity_median",
    "scaffold_cluster_size", "murcko_scaffold",
    # —— 结构 ——
    "std_smiles", "inchikey",
]


def main() -> None:
    d = pd.read_csv(SRC, low_memory=False)
    g = d[(d.set == "gka_candidate") & d.swissadme_ok & d.admetlab_ok].copy()
    print(f"GKA 候选（两个工具都有结果）：{len(g)}")

    sa_bbb = g.swissadme_bbb_permeant == "Yes"
    sa_pgp = g.swissadme_pgp_substrate == "Yes"
    ad_bbb = g.admetlab_bbb > 0.5

    rule_a = sa_bbb & ad_bbb
    rule_b = (g.mw <= 400) & (g.tpsa <= 90) & ad_bbb & ~sa_pgp
    rule_c = g.gka_is_positive_control.fillna(False).astype(bool)
    print(f"  A_both_tools  两个工具都判能入脑：{int(rule_a.sum())}")
    print(f"  B_cns_space   CNS 常见区 + 概率>0.5 + 非 P-gp 底物：{int(rule_b.sum())}")
    print(f"  C_gka_control GKA 身份对照（参照物，非候选）：{int(rule_c.sum())}")
    if int(rule_c.sum()) != 12:
        raise SystemExit(f"[FATAL] GKA 身份对照应为 12 个，实得 {int(rule_c.sum())}")

    hit = rule_a | rule_b | rule_c
    g = g[hit].copy()
    g["rule"] = [
        "+".join([n for n, m in (("A_both_tools", a), ("B_cns_space", b),
                                 ("C_gka_control", c)) if m])
        for a, b, c in zip(rule_a[hit], rule_b[hit], rule_c[hit])
    ]
    # 排序只为看着方便，**不是排名**：A/B 在前、C 组殿后，组内按入脑概率
    g["_is_ctrl"] = rule_c[hit].astype(int)
    g["_n_rule"] = g.rule.str.count(r"\+") + 1 - g._is_ctrl
    g = g.sort_values(["_is_ctrl", "_n_rule", "admetlab_bbb"],
                      ascending=[True, False, False])

    out = g[COLUMNS].copy()
    out.to_csv(OUT, index=False)
    print(f"\n并集 {len(out)} 个分子 → {OUT.name}")
    ab = g[g._is_ctrl == 0]
    print(f"  A/B 观察清单 {len(ab)} 个，其中同时命中 A 和 B 的 {int((ab._n_rule == 2).sum())} 个")
    print(f"  C 组参照 {int(g._is_ctrl.sum())} 个")
    print(f"  骨架数：{ab.murcko_scaffold.nunique()}（A/B 的 {len(ab)} 个分子）")
    print(f"  效力档：{ab.gka_priority.value_counts().sort_index().to_dict()}")
    print("\n⚠ A/B 是待讨论的观察清单，不是筛选结果——规则用的都是未标定的文献习惯值。")
    print("⚠ C 组是「这个分子是不是 GKA」的对照，对入脑这条轴沉默，不得用来标定阈值。")


if __name__ == "__main__":
    main()
