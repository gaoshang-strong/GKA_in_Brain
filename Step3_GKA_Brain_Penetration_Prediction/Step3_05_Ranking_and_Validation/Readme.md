# Step3_05：候选排序与流程验收 —— 🚧 判据未做

## 现在这个目录里有什么

**只有描述性材料**，没有任何判据、阈值或候选名单：

```text
Step3_05_Make_Figures.py            作图脚本
Step3_05_Fig1_Headline.png          两个工具给出的「能进脑」比例，三组并排
Step3_05_Fig2_Resolution.png        每分子一个点：ADMETlab 入脑概率的分布
Step3_05_Fig3_ChemicalSpace.png     MW × TPSA，GKA 候选落在两组对照之外
Step3_05_Fig4_CNS_MPO_Breakdown.png CNS MPO 六项拆解，差距集中在两项
Step3_05_Figures_Explained.md       四张图的解读
Step3_05_Build_Shortlist.py         观察清单的生成脚本
Step3_05_Candidate_Shortlist.csv    23 × 28：A/B 两条规则命中的 11 个 + 12 个 GKA 身份对照
Step3_05_Candidate_Shortlist.md     上表的逐列说明
```

图与清单都**只读** `../Step3_04_Result_Integration/Step3_04_Integrated_Brain_Penetration_Results.csv`。

⚠ **`Step3_05_Candidate_Shortlist.csv` 不是候选名单，是待讨论的观察清单**：
两条规则（`>0.5`、`MW ≤ 400`、`TPSA ≤ 90`）用的都是**未经标定的文献习惯值**。
表里另外 12 行是 **GKA 身份对照**，它们回答的是「是不是 GKA」而非「能不能进脑」，
**只作参照，不得用来标定阈值**。

## 还没做的（本步的正题）

按 CLAUDE.md 的方法论，**必须先验收、再排序**；分不开就是流程不合格，
回去查工具和参数，不是调阈值把结果凑出来。`selfcheck_bbb_controls()` 要写死在脚本里。

待讨论、本轮不预设答案的问题：

1. 用什么判据认定「对照集被分开了」？分离度怎么量、达到多少算合格？
2. 多个工具结论不一致时怎么办？
   BBB 项两个工具高度一致，**P-gp 项近乎相反**（SwissADME 判 GKA 候选 51.8% 是底物，
   ADMETlab `pgp_sub` 中位 0.001），不能简单取交集或多数。
3. `CNS MPO ≥ 4` 这类文献阈值在本项目对照集上成不成立？
   ⚠ 实测它**几乎分不开** B3DB 阳性（中位 4.78）与阴性（4.48），
   且原文自己写了 "not intended to be used purely as a predictor of CNS penetration"。
4. 排序按分子还是按骨架？787 个分子只有 521 个骨架，直接取 top-N 会拿到一串同系物。
5. **外推问题**（图 3）：GKA 候选的 MW 与 TPSA 超过了阴性对照，
   流程在对照集上分得开 ≠ 在 GKA 上同样可靠。
6. CNS MPO 的**口径敏感性**：GKA 候选换口径时分数中位波动 0.49，对照只有 0.16。
   用分数排序前要先看 `cnsmpo_score_variant_spread`。

## 复现图

```bash
/home/sgao30/micromamba/bin/micromamba run -n GKA_in_Brain python Step3_05_Make_Figures.py
```
