对 Step1_04 的 1,333 个分子做方向判定、排序和候选选择，输出一张可复核的 GKA 候选表。

这一步只回答「谁是好的 GCK 激活剂」，**不判断能不能进脑**——物化性质和 BBB 留到 Step2。

输入：`Step1_04_GCK_Activator_Activity_Extraction.csv`（分子汇总）
     `Step1_04_GCK_Activator_Activities.csv`（activity 明细，判方向要用）

## 一、方向判定（并进本步骤当门槛）

Step1_04 的分子是「出现在激活 assay 里的分子」，不是激活剂。先补这道门。

逐分子看 activity 明细，给出 `direction`：

- 效能读数显示酶活未升高 → `no_activation`，**排除**（实测 41 个）
- 同一分子既有激活又有降低 → `conflict`，进榜但打标，不静默取最优（7 个）
- 其余 → `activation`（506 个）

只有 EC50、没有效能记录的分子**不排除**，记 `activation_by_potency`（654 个）。
EC50 出自激活读数的 assay，本身就是方向证据。这条要写进报告，否则口径会被误读。

**判方向不能对 fold 设全局阈值**，比值型读数的分母不统一（见 CLAUDE.md）。
做法是逐 (assay, standard_type, scale) 组解析基线，35 个组的解析结果连同证据句
全量写进报告供复核；解析不出来的标 `uncertain`，不参与方向判定（保险的一侧）。
三类必须区分：

- `control` 21 组：分母是未处理对照 → fold 基线 1 / percent 基线 100
- `reference` 8 组：分母是参比激活剂（如 Ro-28-1675）→ 基线 0，1.0 表示与参比等效
- `kinetic` 3 组：Km / Vmax 类，极性相反或对 K 型激活剂不敏感 → 不判方向。
  `CHEMBL3095514` 是 "decrease in enzyme Km"，值 0.05–0.11 是**强激活**，
  按 fold ≤ 1.05 会被误杀

## 二、排除与打标

不静默丢行，每条都写明理由，保留在输出里可复核：

| 情况 | 处理 | 规模 |
|---|---|---|
| `no_activation` | 排除 | 41 |
| 只有删失 `>` 效力值 | 不进榜，单列保留 | 74 |
| 跨 assay 效力分歧 >10 倍 | 进榜 + 打标 | 29 |
| `potential_duplicate` | 进榜 + 打标 | 30 |
| 方向冲突 | 进榜 + 打标 | 7 |

## 三、排序

**单一排序键 `pactivity_median`，不合成总分。**

用中位数不用最大值：1,088 个分子只有一个 assay，两者相等；对其余多 assay 的分子，
中位数抗单点异常。`pactivity_max` 并列展示，不参与排序。

**效能和证据只作并列列，不折进排序键。** 证据等级衡量的是「在 ChEMBL 里被测了多少次」，
不是分子有多好——MK-0941、AZD-1656、PF-04991532 三个 phase 2 药的证据都是「弱」，
而参比化合物 Ro-28-1675 是「强」。把证据算进排序，等于系统性地把最像药的分子往后排。

## 四、分层

不排一个大榜，分三层，层内按排序键排：

- A 层：pAct 中位 ≥ 6.5 **且**有正向效能佐证（169 个）
- B 层：pAct 中位 ≥ 6.0（613 个）
- C 层：其余有可定量效力的（225 个）
- 另有 215 个有激活证据但无可定量效力值，单列不排序

**阈值由阳性对照标定，不按整体分布拍。** 最初取 7.0 / 6.5，自检只过 4/6——
Ro-28-1675 的 EC50 跨实验是 127–690 nM（pAct 中位 6.39），Piraglitin 是
364–6320 nM（中位 6.145），已知临床 GKA 本来就在几百 nM 这一档，原阈值把临床药
从中间劈开了。

**A 层不是「比 B 层好」，是「双证」**：MK-0941、AZD-1656 效力更强（8.02 / 7.55）
但 ChEMBL 里没有效能读数，只测了一个轴，落在 B。

## 五、骨架去冗余

521 个 Murcko 骨架，最大一簇 46 个分子，单例骨架 347 个，前 20 个骨架覆盖 348 个分子。
不做骨架聚类，top-50 会是两篇 SAR 论文的同系物列表——看着 50 个，其实 2 个化学起点。

**加列不删行**（保持一行一个分子）：`murcko_scaffold`、`scaffold_cluster_size`、
`is_scaffold_representative`（簇内效力最优者）。

## 六、自检

6 个已知临床 GKA 作阳性对照，必须落进 A/B 层：

`Ro-281675`、`PIRAGLIATIN`、`NERIGLIATIN`、`PF-04991532`、`AZD-1656`、`MK-0941`

落不进说明规则有问题，不是数据有问题。结果写进报告。

## 输出

- `Step1_05_GCK_Activator_Candidates.csv` —— 候选主表，一行一个分子，
  含 `direction` / `direction_evidence` / `exclude_reason` / `tier` / 排序键 / 骨架列
- `Step1_05_GCK_Activator_Ranking_and_Candidate_Selection.md` —— 报告，
  含出处、各门槛的排除计数、分层规模、骨架分布、阳性对照自检结果

代码：`Step1_05_GCK_Activator_Ranking_and_Candidate_Selection.py`
