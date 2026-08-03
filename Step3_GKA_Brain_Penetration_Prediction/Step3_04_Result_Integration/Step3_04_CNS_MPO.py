#!/usr/bin/env python
"""
CNS MPO 打分函数（Pfizer，Wager 2010）

**拐点与公式全部逐字取自原文 PDF**，不是凭记忆写的：

  Wager TT, Hou X, Verhoest PR, Villalobos A.
  "Moving beyond Rules: The Development of a Central Nervous System Multiparameter
  Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties."
  ACS Chem Neurosci 2010;1(6):435-449. doi:10.1021/cn100008c  (PMID 22778837)
  → 本仓库 Step3_GKA_Brain_Penetration_Prediction/cn100008c.pdf

出处逐条对应：

* **Table 1（p.438）** "The CNS MPO Properties, Functions, Weighting, Value Range and
  Parameter Ranges" —— 六项的变换类型、权重（全部 1.0）、最优与最差区间
* **Figure 4（p.439）** —— 六条曲线上标出的拐点数值（绿箭头 = T0 1.0，红箭头 = T0 0.0）
* **Methods / eq 1（p.447）** —— 拐点之间是**线性**插值的分段函数：

      T(x) = y1                                        x ≤ x1
           = y_{i-1} + (y_i - y_{i-1})/(x_i - x_{i-1}) * (x - x_{i-1})   x_{i-1} < x ≤ x_i
           = y_n                                        x > x_n

* **eq 2（p.447）** —— 总分 D = Σ w_k · T_k，六项 **w = 1.0**，故 0–6。

⚠ 原文明确写了一句话，Step3_05 用这个分数时不能忽略（p.446 右栏）：
  "the algorithm is not intended to be used purely as a predictor of CNS penetration"
  ——CNS MPO 是**成药性对齐**的设计工具，不是入脑预测器。

原文用的性质计算包（Methods, p.448）与本项目的替代来源：

| 项 | 原文用 | 本项目用 | 备注 |
|---|---|---|---|
| ClogP | BioByte CLOGP | ADMETlab `logP`（主）/ RDKit Crippen / SwissADME | 三种口径都算，见 VARIANTS |
| ClogD7.4 | ACD/Labs | ADMETlab `logD` | 唯一来源 |
| pKa（最碱性） | ACD/Labs | ADMETlab `pka_basic` | 唯一来源 |
| TPSA | **Ertl 2000（ref 9）** | RDKit `tpsa`（N/O 口径，主）| 含 S/P 的口径另算一版 |
| MW | — | RDKit `mw`（**平均分子量**，不是单同位素） | |
| HBD | — | RDKit `hbd` | |

**TPSA 口径的选择要显式**：原文 ref 9 是 Ertl/Rohde/Selzer 2000，通行实现（含 RDKit 默认）
只计 N/O，所以主口径取 `tpsa`；含 S/P 的 `tpsa_sandp` 另算一版做敏感性。
这与 BOILED-Egg 相反（那个明确用含 S/P 的），两处别搞混。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# 拐点：Table 1 + Figure 4。列表是 (x, y) 序对，x 递增，y 为该点的 T0。
# 单调下降 = 两个拐点；TPSA 的驼峰 = 四个拐点。
# --------------------------------------------------------------------------- #
INFLECTIONS: dict[str, list[tuple[float, float]]] = {
    "clogp": [(3.0, 1.0), (5.0, 0.0)],          # Table 1: ≤3 最优，>5 最差
    "clogd": [(2.0, 1.0), (4.0, 0.0)],          # Table 1: ≤2 最优，>4 最差
    "mw":    [(360.0, 1.0), (500.0, 0.0)],      # Table 1: ≤360 最优，>500 最差
    "tpsa":  [(20.0, 0.0), (40.0, 1.0),         # Table 1: 40<TPSA≤90 最优，
              (90.0, 1.0), (120.0, 0.0)],       #          ≤20 与 >120 最差
    "hbd":   [(0.5, 1.0), (3.5, 0.0)],          # Table 1: ≤0.5 最优，>3.5 最差
    "pka":   [(8.0, 1.0), (10.0, 0.0)],         # Table 1: ≤8 最优，>10 最差
}
WEIGHTS = {k: 1.0 for k in INFLECTIONS}         # Table 1: 六项权重全是 1.0
PROPERTIES = list(INFLECTIONS)


def transform(x, prop: str):
    """eq 1 的分段线性变换。x 可以是标量或 array/Series，NaN 原样传出。"""
    pts = INFLECTIONS[prop]
    xs = np.array([p[0] for p in pts], dtype=float)
    ys = np.array([p[1] for p in pts], dtype=float)
    v = np.asarray(x, dtype=float)
    # np.interp 在区间内做线性插值，两端自动取 ys[0] / ys[-1]，正与 eq 1 一致
    out = np.interp(v, xs, ys)
    return np.where(np.isnan(v), np.nan, out)


def score_frame(mw, clogp, clogd, tpsa, hbd, pka) -> pd.DataFrame:
    """给六项输入，返回六个 T0 列 + 总分列（eq 2，w 全为 1）。

    任一项缺失则总分为 NaN——**不做补零**：补零会把「没测到」说成「最差」。
    """
    t0 = pd.DataFrame({
        "t0_clogp": transform(clogp, "clogp"),
        "t0_clogd": transform(clogd, "clogd"),
        "t0_mw": transform(mw, "mw"),
        "t0_tpsa": transform(tpsa, "tpsa"),
        "t0_hbd": transform(hbd, "hbd"),
        "t0_pka": transform(pka, "pka"),
    })
    cols = [f"t0_{p}" for p in PROPERTIES]
    t0["score"] = sum(WEIGHTS[p] * t0[f"t0_{p}"] for p in PROPERTIES)
    t0.loc[t0[cols].isna().any(axis=1), "score"] = np.nan
    return t0


# --------------------------------------------------------------------------- #
# 自检：拿原文里给了「输入值 → T0」的三张表回算，对不上就报错。
# 这是 CLAUDE.md 那条方法论的落实——自检写死在脚本里，不靠人记得跑。
# --------------------------------------------------------------------------- #
# Table 4（p.446）"Active CNS MPO Calculator"：**唯一一组输入未被四舍五入的算例**，
# 六项 T0 与总分都能对到小数点后两位。
PAPER_TABLE4 = {
    "name": "Table 4 CNS MPO Calculator",
    "inputs": dict(clogp=3.1, clogd=1.7, tpsa=21.0, mw=392.5, hbd=1.0, pka=9.2),
    "t0": dict(clogp=0.95, clogd=1.00, tpsa=0.05, mw=0.77, hbd=0.83, pka=0.40),
    "score": 4.0, "tol": 0.005,
}
# Table 3（p.440）三个辉瑞候选：表里印的输入值本身是四舍五入过的（如 ClogP 3.8 → T0 0.58
# 反推真值 3.84），所以容差放到 0.03。
PAPER_TABLE3 = [
    {"name": "PF-02545920 (PDE10)",
     "inputs": dict(clogp=3.8, clogd=3.5, tpsa=52.8, mw=392.5, hbd=0.0, pka=4.3),
     "t0": dict(clogp=0.58, clogd=0.24, tpsa=1.00, mw=0.77, hbd=1.00, pka=1.00),
     "score": 4.6, "tol": 0.03},
    {"name": "PF-03654746 (H3)",
     "inputs": dict(clogp=2.4, clogd=0.0, tpsa=32.3, mw=322.4, hbd=1.0, pka=9.2),
     "t0": dict(clogp=1.00, clogd=1.00, tpsa=0.62, mw=1.00, hbd=0.83, pka=0.42),
     "score": 4.9, "tol": 0.03},
    {"name": "PF-04447943 (PDE9)",
     "inputs": dict(clogp=-1.5, clogd=-0.7, tpsa=101.9, mw=395.4, hbd=1.0, pka=7.9),
     "t0": dict(clogp=1.00, clogd=1.00, tpsa=0.60, mw=0.75, hbd=0.83, pka=1.00),
     "score": 5.2, "tol": 0.03},
]


def selfcheck(verbose: bool = True) -> list[str]:
    """用原文给的算例回算六条曲线。返回不通过的说明列表（空 = 全过）。"""
    problems = []
    for case in [PAPER_TABLE4] + PAPER_TABLE3:
        tol = case["tol"]
        got = {p: float(transform(case["inputs"][p], p)) for p in PROPERTIES}
        total = sum(got.values())
        for p in PROPERTIES:
            if abs(got[p] - case["t0"][p]) > tol:
                problems.append(f"{case['name']} 的 T0_{p}：算得 {got[p]:.3f}，"
                                f"原文 {case['t0'][p]:.2f}（容差 {tol}）")
        if abs(total - case["score"]) > tol * 6:
            problems.append(f"{case['name']} 的总分：算得 {total:.3f}，原文 {case['score']}")
        if verbose and not problems:
            print(f"      ✓ {case['name']}：六项 T0 与总分 {total:.2f} 均与原文一致")
    return problems


if __name__ == "__main__":
    print("CNS MPO 自检（对照 cn100008c.pdf 的 Table 3 / Table 4）：")
    probs = selfcheck()
    if probs:
        for x in probs:
            print(f"  ⚠ {x}")
        raise SystemExit(1)
    print("全部通过。")
