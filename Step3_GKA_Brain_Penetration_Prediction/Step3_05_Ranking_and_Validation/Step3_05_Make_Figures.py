#!/usr/bin/env python
"""
Step3_05 的描述性图：入脑对照 vs GKA 候选

读 Step3_04 的整合表，画四张图。**图里没有任何判据**——
不画阈值线、不标「合格/不合格」，只把数据摆出来看。阈值属 Step3_05 的讨论。

分组口径（四张图统一）：
  实测能进脑  = B3DB control_positive ∪ Fridén Kp,uu ≥ 0.3
  实测不能进脑 = B3DB control_negative ∪ Fridén Kp,uu ≤ 0.05
  中间带      = Fridén 0.05 < Kp,uu < 0.3（只在图 2 出现，画成灰色）
  GKA 候选    = 787 个

用法：
  micromamba run -n GKA_in_Brain python Step3_05_Make_Figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = (HERE.parent / "Step3_04_Result_Integration"
       / "Step3_04_Integrated_Brain_Penetration_Results.csv")

# 取自 dataviz skill 的参考调色板：分类槽 1/2/3（文档已验证 all-pairs 双模式通过）
C_GKA = "#2a78d6"        # slot 1 blue   —— GKA 候选（本项目的主角）
C_NEG = "#eb6834"        # slot 2 orange —— 实测不能进脑
C_POS = "#1baf7a"        # slot 3 aqua   —— 实测能进脑
C_MID = "#898781"        # muted         —— 中间带（不是分类槽，故意用中性灰）

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

mpl.rcParams.update({
    "font.family": ["Noto Sans CJK SC", "Noto Sans CJK JP", "DejaVu Sans"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": MUTED, "ytick.color": INK2,
    "axes.grid": False, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "figure.dpi": 200,
})


def tidy(ax, xgrid=False, ygrid=False):
    ax.tick_params(length=0)
    if xgrid:
        ax.set_axisbelow(True)
        ax.xaxis.grid(True, color=GRID, lw=0.8)
    if ygrid:
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID, lw=0.8)


def load() -> pd.DataFrame:
    d = pd.read_csv(SRC, low_memory=False)
    grp = pd.Series("", index=d.index, dtype=object)
    grp[d.set == "gka_candidate"] = "GKA 候选"
    grp[d.control_class == "control_positive"] = "实测能进脑"
    grp[d.control_class == "control_negative"] = "实测不能进脑"
    fr = d.set == "bbb_control_friden"
    grp[fr & (d.kpuu_brain >= 0.3)] = "实测能进脑"
    grp[fr & (d.kpuu_brain <= 0.05)] = "实测不能进脑"
    grp[fr & (d.kpuu_brain > 0.05) & (d.kpuu_brain < 0.3)] = "中间带"
    grp[grp == ""] = "未分组"          # 实测值缺失的，不硬塞进任何一组
    d["grp"] = grp
    d = d[d.swissadme_ok & d.admetlab_ok].copy()
    un = d[d.grp == "未分组"]
    if len(un):
        print(f"  ⚠ {len(un)} 个对照没有实测数值、不入任何一组："
              f"{', '.join(un.mol_id + '(' + un.compound_name.fillna('?') + ')')}")
    return d


# --------------------------------------------------------------------------- #
def fig1_headline(d: pd.DataFrame) -> None:
    """两个工具各自判「能进脑」的比例。三组并排。"""
    order = ["实测能进脑", "实测不能进脑", "GKA 候选"]
    colors = [C_POS, C_NEG, C_GKA]
    sa = [(d[d.grp == g].swissadme_bbb_permeant == "Yes").mean() * 100 for g in order]
    ad = [(d[d.grp == g].admetlab_bbb > 0.5).mean() * 100 for g in order]
    ns = [int((d.grp == g).sum()) for g in order]

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, vals, title in zip(axes, (sa, ad),
                               ("SwissADME　判为「能进脑」", "ADMETlab　入脑概率 > 0.5")):
        y = np.arange(len(order))[::-1]
        ax.barh(y, vals, height=0.55, color=colors, zorder=3)
        for yy, v in zip(y, vals):
            ax.text(v + 2, yy, f"{v:.0f}%", va="center", ha="left",
                    color=INK2, fontsize=11)
        ax.set_yticks(y, [f"{g}\n({n} 个)" for g, n in zip(order, ns)], fontsize=10)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100], ["0", "25%", "50%", "75%", "100%"])
        ax.set_title(title, color=INK, pad=12, loc="left")
        tidy(ax, xgrid=True)
    fig.suptitle("两个独立工具给出的比例：GKA 候选远低于两组对照",
                 x=0.02, ha="left", fontsize=15, weight="bold")
    fig.text(0.02, -0.04, "对照分子的进脑与否是实验测过的；GKA 候选没有实测数据，"
                          "柱高是模型预测的比例。", color=MUTED, fontsize=9.5)
    fig.tight_layout(rect=[0, 0.02, 1, 0.90])
    fig.savefig(HERE / "Step3_05_Fig1_Headline.png", bbox_inches="tight")
    plt.close(fig)


def fig2_resolution(d: pd.DataFrame) -> None:
    """每个分子一个点：模型给的入脑概率。看对照分不分得开、GKA 落在哪。"""
    rows = [("实测能进脑", C_POS), ("中间带", C_MID),
            ("实测不能进脑", C_NEG), ("GKA 候选", C_GKA)]
    fig, ax = plt.subplots(figsize=(10, 3.8))
    rng = np.random.default_rng(0)
    for i, (g, c) in enumerate(rows):
        v = d.loc[d.grp == g, "admetlab_bbb"].dropna().values
        y = len(rows) - 1 - i + rng.uniform(-0.17, 0.17, len(v))
        ax.scatter(v, y, s=14, color=c, alpha=0.45, linewidths=0, zorder=3)
        med = np.median(v)
        ax.plot([med, med], [len(rows) - 1 - i - 0.30, len(rows) - 1 - i + 0.30],
                color=c, lw=2.5, solid_capstyle="round", zorder=4)
        ax.text(1.02, len(rows) - 1 - i, f"中位 {med:.2f}", va="center",
                color=INK2, fontsize=10)
    ax.set_yticks(range(len(rows))[::-1],
                  [f"{g}\n({int((d.grp == g).sum())} 个)" for g, _ in rows], fontsize=10)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0],
                  ["0\n进不去", "0.25", "0.5", "0.75", "1.0\n进得去"])
    ax.set_xlabel("ADMETlab 给的入脑概率", labelpad=8)
    ax.set_title("两组对照落在概率轴的两端，GKA 候选整体贴在最左边",
                 loc="left", pad=12, fontsize=15, weight="bold")
    tidy(ax, xgrid=True)
    fig.text(0.01, -0.10, "一个点是一个分子，竖线是各组中位数。"
                          "「中间带」是实测值卡在两者之间的 9 个分子。",
             color=MUTED, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "Step3_05_Fig2_Resolution.png", bbox_inches="tight")
    plt.close(fig)


def fig3_space(d: pd.DataFrame) -> None:
    """为什么会这样：分子的「大小」和「极性」把三组分开了。"""
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.add_patch(plt.Rectangle((40, 100), 50, 360 - 100, facecolor="#1baf7a",
                               alpha=0.07, zorder=1, linewidth=0))
    ax.annotate("CNS MPO 认为最理想的一角\nTPSA 40–90 且 MW ≤ 360",
                xy=(90, 360), xytext=(128, 250), ha="left", color=INK2, fontsize=9.5,
                zorder=5, arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    for g, c, s, a in (("GKA 候选", C_GKA, 13, 0.40),
                       ("实测能进脑", C_POS, 15, 0.55),
                       ("实测不能进脑", C_NEG, 15, 0.60)):
        sub = d[d.grp == g]
        ax.scatter(sub.tpsa, sub.mw, s=s, color=c, alpha=a, linewidths=0,
                   label=f"{g}（{len(sub)}）", zorder=3)
    ax.set_xlim(0, 220)
    ax.set_ylim(100, 700)
    ax.set_xlabel("TPSA　极性表面积（越大越「怕油」，越难穿过细胞膜）", labelpad=8)
    ax.set_ylabel("MW　分子量（越大越难挤进去）", labelpad=8)
    ax.set_title("GKA 候选又大又极性，落在两组对照之外",
                 loc="left", pad=12, fontsize=15, weight="bold")
    leg = ax.legend(frameon=False, loc="upper right", fontsize=10,
                    handletextpad=0.4, labelspacing=0.6)
    for h in leg.legend_handles:
        h.set_alpha(0.9)
        h.set_sizes([34])
    tidy(ax, xgrid=True, ygrid=True)
    fig.text(0.01, -0.03, "一个点是一个分子。GKA 候选不是落在两类对照之间，"
                          "而是整体越过了阴性对照——这是外推，不是内插。",
             color=MUTED, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "Step3_05_Fig3_ChemicalSpace.png", bbox_inches="tight")
    plt.close(fig)


def fig4_mpo(d: pd.DataFrame) -> None:
    """CNS MPO 拆成六项，看 GKA 到底输在哪两项。"""
    props = [("cnsmpo_t0_mw", "分子量 MW"), ("cnsmpo_t0_tpsa", "极性 TPSA"),
             ("cnsmpo_t0_clogd", "脂溶性 LogD"), ("cnsmpo_t0_clogp", "脂溶性 LogP"),
             ("cnsmpo_t0_hbd", "氢键给体 HBD"), ("cnsmpo_t0_pka", "碱性 pKa")]
    groups = [("实测能进脑", C_POS), ("实测不能进脑", C_NEG), ("GKA 候选", C_GKA)]
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    h = 0.24
    y0 = np.arange(len(props))[::-1]
    for k, (g, c) in enumerate(groups):
        vals = [d.loc[d.grp == g, col].median() for col, _ in props]
        ax.barh(y0 + (1 - k) * h, vals, height=h - 0.04, color=c, zorder=3,
                label=f"{g}（总分中位 {d.loc[d.grp == g, 'cnsmpo_score'].median():.1f}）")
    ax.set_yticks(y0, [n for _, n in props], fontsize=10.5)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.5, 1.0], ["0\n最不理想", "0.5", "1.0\n最理想"])
    ax.set_xlabel("该项的得分（0–1，六项相加就是 CNS MPO 总分 0–6）", labelpad=8)
    ax.set_title("六项拆开看：差距集中在「分子量」和「极性」这两项",
                 loc="left", pad=42, fontsize=15, weight="bold")
    ax.legend(frameon=False, fontsize=10, ncol=3, loc="lower left",
              bbox_to_anchor=(0, 1.005), columnspacing=1.6, handletextpad=0.5)
    tidy(ax, xgrid=True)
    fig.text(0.01, -0.05, "柱高是各组的中位得分。其余四项三组几乎一样高——"
                          "差距集中在两项上，不是全面落后。", color=MUTED, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "Step3_05_Fig4_CNS_MPO_Breakdown.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    d = load()
    print(d.grp.value_counts().to_dict())
    fig1_headline(d)
    fig2_resolution(d)
    fig3_space(d)
    fig4_mpo(d)
    for p in sorted(HERE.glob("Step3_05_Fig*.png")):
        print(f"  {p.name}  {p.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
