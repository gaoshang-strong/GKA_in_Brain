# Step3_04：SwissADME 与 ADMETlab 结果合并

## 出处

| | |
|---|---|
| 运行时间 | 2026-08-03 08:11:31 |
| 脚本 | `Step3_04_Merge_Results.py` |
| 骨架输入 | `Step3_01_RDKit_Processed.csv`（1,274 行 × 67 列） |
| SwissADME | <https://www.swissadme.ch/>，网页版**不标版本号**，7 批结果文件 mtime 2026-08-02 22:59–2026-08-03 00:00 |
| ADMETlab | <https://admetlab3.scbdd.com/> **3.0**，16 批结果文件 mtime 2026-08-03 00:0x–01:0x |
| 环境 | micromamba `GKA_in_Brain` / python 3.11.15 / rdkit 2026.03.4 / pandas 3.0.5 |

各批原始结果文件的行数、mtime 与 SHA256 逐批记在 `Step3_04_Batch_Provenance.csv`（23 行）。

**⚠ 本步的产物是数据，不是结论。** 这里不下「流程是否可用」的判定，
也不下「哪些 GKA 能入脑」的判定——两者都属 Step3_05。
下文所有分布数字都只是描述，没有任何阈值被应用。

---

## 产物

```text
Step3_04_Integrated_Brain_Penetration_Results.csv   1,274 行 × 261 列   ← 主产物
Step3_04_Batch_Provenance.csv                       23 行   每批出处与校验计数
Step3_04_Verification_Failures.csv                  1 行    被剔除的结果行
Step3_04_Anchor_Drift.csv                           40 行   锚点重复间的分歧量
Step3_04_summary.json                                       机器可读汇总
```

261 列的构成：

| 块 | 列数 | 内容 |
|---|---:|---|
| Step3_01 骨架 | 67 | `mol_id` / `set` / GKA 元数据 / 对照实测值 / 标准化结构 / RDKit 描述符 |
| `swissadme_*` | 54 | 48 个结果列 + 6 个覆盖与追溯列 |
| `admetlab_*` | 128 | 122 个结果列 + 6 个覆盖与追溯列 |
| CNS MPO | 12 | 六项 `cnsmpo_t0_*` + 总分 + 3 个换口径版本 + spread + 齐备标记 |

每个工具的 6 个追溯列：`*_ok`（是否拿到结果）、`*_n_replicates`（该结构被跑了几次）、
`*_batches`（出现在哪几批，JSON 列表）、`*_result_from`（结果来自哪个 mol_id 的提交）、
`*_shared_result_with`（与谁共用一次提交）、`*_missing_reason`（没结果时写明原因）。

---

## 做法

### 1. 逐行核对结构，不符的整行剔除

CLAUDE.md 记着一条事故：SwissADME 曾返回**从未提交过的分子**
（batch6 的 `B3D_0441` 返回了樟脑），名称是对的、只有结构错了，光核对名称发现不了。
所以两个工具都走同一道校验：

* **SwissADME**：按 `(批次, Molecule 名称)` 对回提交清单 → 比 InChIKey。
* **ADMETlab**：不接受名称字段，按 `(批次, 行序)` 对回提交清单 → 用它原样返回的
  `raw_smiles` 比 InChIKey。行序与结构是两条独立路径，行序错位会被结构校验抓住。

### 2. 锚点重复合并

20 个对照分子在每批重复提交，同一结构因此有 6–16 次结果。
合并规则：**数值取中位数、文本取众数**，重复间的分歧量逐列写进 `Step3_04_Anchor_Drift.csv`。

### 3. 回填到 1,274 行

用 Step3_02/03 各自的 `Result_Join_Map.csv` 里的 `submitted_as` 列贴回。
去盐后结构相同的分子共用一次提交，回填时按结构分发，**行不合并**。

### 4. 加列不删行

1,274 行进、1,274 行出。没拿到结果的行留在表里，`*_ok = False` 且写明原因，
**不填任何默认值**。

---

## 校验结果

### 结构核对

| 工具 | 返回行数 | 通过 | 剔除 | 唯一结构 |
|---|---:|---:|---:|---:|
| SwissADME | 1,391 | 1,390 | **1** | 1,271 |
| ADMETlab | 1,571 | 1,571 | 0 | 1,271 |

唯一被剔除的一行（`Step3_04_Verification_Failures.csv`）：

| 批 | mol_id | 提交的结构 | 返回的结构 |
|---|---|---|---|
| SwissADME batch6 | `B3D_0441` trifluopromazine | `XSCGXQMFQXDFCW-UHFFFAOYSA-N` | `LHXDLQBQYFFVNW-XCBNKYQSSA-N`（樟脑） |

这正是 CLAUDE.md 里记录的那一行，脚本独立重现了它。
**剔除它之后有实质影响**：`B3D_0441` 是 `control_positive`，
错误的那一行给出 `BBB permeant = Yes`，其余 6 次重复都是 `No`；
若不剔除，众数合并会把它拉成 Yes。

ADMETlab 侧 1,571 行**全部通过**：行序 16/16 批逐行对上，`raw_smiles` 与提交的
`std_smiles` **逐字一致 1,571/1,571**，InChIKey 全同。

### 锚点漂移

| 工具 | 锚点数 | 重复次数 | 比较列数 | 有分歧的锚点 | 最大绝对差 |
|---|---:|---:|---:|---:|---|
| SwissADME | 20 | 7（`B3D_0441` 因剔除只有 6） | 48 | **0** | 0 |
| ADMETlab | 20 | 16 | 122 | 2 | **5.96e-08** |

ADMETlab 那两处分歧落在 `MRP1` / `EI` / `Neurotoxicity-DI`，差在小数点后第 8 位，
是浮点噪声不是模型漂移。**两个工具在本轮的取数窗口内都没有发生批次漂移。**

### 覆盖

| set | 行数 | SwissADME | ADMETlab |
|---|---:|---:|---:|
| `gka_candidate` | 787 | 787 | 787 |
| `bbb_control_b3db` | 445 | 444 | 444 |
| `bbb_control_friden` | 42 | 42 | 42 |
| **合计** | **1,274** | **1,273** | **1,273** |

唯一的缺口是 `B3D_0012`（mdl 63,246，MW 1802.7 / 127 重原子），
**人工排除，不是工具失败**——记录在 `Step3_02_Manual_Exclusions.csv`，
`*_missing_reason` 列写明「人工排除：分子过大，超出 SwissADME 适用范围」。

两行与别人共用一次提交（去盐后同结构，上游 `dup_group` 已判定）：

| 行 | 共用对象 | 关系 |
|---|---|---|
| `GKA_0783` MK-0941（药物条目） | `GKA_0057` MK-0941 FREE BASE | 盐 / 游离碱 |
| `GKA_0786` Globalagliatin 盐酸盐 | `GKA_0784` LY-2608204 | 盐 / 游离碱 |

### 脚本自检

全部通过，`selfcheck_problems` 为空：行数 1,274 进 1,274 出、`mol_id` 无重复、
12 个 GKA 身份对照一个不少、每个缺结果的行都有原因、
两个工具的 MW 与本地 RDKit 值逐行一致。

---

## ⚠ 本步新发现的一个坑：两个工具的 MW 口径不同

初版自检报了「`admetlab_mw` 与本地 `mw` 相差 >1 Da 的有 73 行，回填可能错位」。
查下来**不是错位**：

* **SwissADME 的 `MW` 是平均分子量**（与 RDKit `Descriptors.MolWt` 一致，最大差 0.02 Da）
* **ADMETlab 的 `MW` 是单同位素质量**（与 RDKit `ExactMolWt` 一致，最大差 0.005 Da）

73 行差异全部落在含 Cl / Br 的分子上（Cl 45 个、Br 32 个），最大差 3.07 Da——
氯溴的天然同位素分布让平均质量明显高于单同位素质量。

**影响**：CNS MPO 的第一项是 MW，用哪个口径要显式选择；
差 1–3 Da 对 MPO 分数影响极小，但**拿 `admetlab_mw` 当平均分子量去比 Lipinski 500
之类的硬边界，含卤分子会系统性偏低**。整合表里两个口径都在（`mw` / `mw_exact` /
`swissadme_mw` / `admetlab_mw`），选哪个是 Step3_05 的事。

---

## CNS MPO：已算（拐点取自原文 PDF）

拐点与公式**逐条取自本地 PDF** `cn100008c.pdf`（Table 1 + Figure 4 给拐点，
Methods eq 1/2 给「拐点间线性插值、六项等权求和」），实现在 `Step3_04_CNS_MPO.py`：

| 项 | 形状 | T0 = 1.0 | T0 = 0.0 | 本项目取哪一列 | 非空 |
|---|---|---|---|---|---:|
| cLogP | 单调下降 | ≤ 3 | > 5 | `admetlab_logp` | 1,273 |
| cLogD7.4 | 单调下降 | ≤ 2 | > 4 | `admetlab_logd` | 1,273 |
| MW | 单调下降 | ≤ 360 | > 500 | `mw`（**平均**分子量） | 1,274 |
| TPSA | **驼峰** | 40 < TPSA ≤ 90 | ≤ 20 或 > 120 | `tpsa`（**N/O 口径**） | 1,274 |
| HBD | 单调下降 | ≤ 0.5 | > 3.5 | `hbd` | 1,274 |
| 最碱性 pKa | 单调下降 | ≤ 8 | > 10 | `admetlab_pka_basic` | 1,273 |

**六项同时齐备 1,273 行、出分 1,273 行**（缺的一行就是人工排除的 `B3D_0012`）。
任一项缺失总分即 NaN，**不补零**——补零等于把「没测到」说成「最差」。

**曲线自检写死在脚本里**：原文 **Table 4（p.446）那个算例的输入没有被四舍五入**，
六项 T0（0.95 / 1.00 / 0.05 / 0.77 / 0.83 / 0.40）与总分 **4.0** 全部复现到小数点后两位；
再加 Table 3 的三个辉瑞候选（输入是舍入过的，容差 0.03）。四组任一项对不上，脚本直接退出。

**三处必须显式选择的口径**（原文用 BioByte ClogP、ACD logD 与 pKa，本项目都没有）：

| 项 | 选了什么 | 为什么 |
|---|---|---|
| ClogP | ADMETlab `logP` | 与同源的 logD / pKa 保持一致，不跨工具混脂溶性标尺 |
| TPSA | **N/O 口径** `tpsa` | 原文 ref 9 是 Ertl 2000，通行实现只计 N/O；**与 BOILED-Egg 相反**，别混 |
| MW | 平均分子量 `mw` | **不是** ADMETlab 那个单同位素质量（见上一节） |

换口径的三个版本一并算出：`cnsmpo_score_tpsa_sandp` / `_logp_rdkit` / `_logp_swiss`，
并给逐分子的 `cnsmpo_score_variant_spread`。
⚠ **GKA 候选的口径敏感性明显更高**（spread 中位 0.49，B3DB 对照只有 0.16、Fridén 0.00）——
Step3_05 拿 MPO 分数排序前应先看这一列。

### 分数分布（数据，不作判定）

| 集合 | 中位 | ≥ 4 | ≥ 5 |
|---|---:|---:|---:|
| B3DB 对照 | 4.74 | 75.3% | 41.1% |
| Fridén | 4.80 | 69.0% | 33.3% |
| **GKA 候选** | **3.78** | **38.5%** | **6.9%** |

六项 T0 的中位显示**差距只集中在两项**：MW（对照 1.00 / GKA 0.27）与 TPSA（1.00 / 0.49）；
LogP 与 pKa 三组都是 1.00，HBD 上 GKA 反而更好（0.83 vs 阴性对照 0.50）。

⚠⚠ **CNS MPO 在本项目的入脑对照上分辨力有限**：
B3DB `control_positive` 中位 4.78、`control_negative` 4.48——两组几乎重叠
（而两个工具的 BBB 项把这两组分得很开）。这不是实现问题，原文自己写过（p.446）：

> "the algorithm is **not intended to be used purely as a predictor of CNS penetration**"

它衡量的是**成药性对齐**（性质像不像 CNS 药），不是入脑的物理预测。
**Step3_05 若要用它，必须先面对这件事。**

原文（题名与卷期页已在 PubMed 逐条核对，PDF 在 Step3 根目录）：

* Wager TT, Hou X, Verhoest PR, Villalobos A.
  **"Moving beyond Rules: The Development of a Central Nervous System Multiparameter
  Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties."**
  *ACS Chem Neurosci* 2010;1(6):435–449. doi:10.1021/cn100008c（PMID 22778837）
* Wager TT, Hou X, Verhoest PR, Villalobos A.
  **"Central Nervous System Multiparameter Optimization Desirability:
  Application in Drug Discovery."**
  *ACS Chem Neurosci* 2016;7(6):767–775. doi:10.1021/acschemneuro.6b00029（PMID 26991242）
* 同期配套文（性质空间的数据来源）：Wager TT, Chandrasekaran RY, Hou X, Troutman MD,
  Verhoest PR, Villalobos A, Will Y. **"Defining Desirable Central Nervous System Drug
  Space through the Alignment of Molecular Properties, in Vitro ADME, and Safety
  Attributes."** *ACS Chem Neurosci* 2010;1(6):420–434. doi:10.1021/cn100007x

本地 PDF：`../cn100008c.pdf`（MPO 原文）与 `../cn100007x.pdf`（配套文）。

---

## ⚠ SwissADME 有 6 列整列为空

`iLOGP`、`CYP1A2/2C19/2C9/2D6/3A4 inhibitor` 六列在 1,390 条结果里**全部是 `n/d`**，
脚本把 `n/d` 转成空值，列保留在表里（列名见 `summary.json` 的 `all_nd_columns`）。
这是工具侧的既有状况（其首页公告 ChemAxon 停止对学术网站支持后的改版），不是本步丢数。
**副作用**：`Consensus Log P` 现在是四者均值（XLOGP3/WLOGP/MLOGP/Silicos-IT），不再是五者。

---

## 描述性统计（**只是数据，判定属 Step3_05**）

以下都基于两个工具都有结果的 1,273 行。

### 三个集合并排

| 集合 | n | SwissADME BBB+ | Pgp 底物 | BBB+ 且非 Pgp | ADMETlab BBB 中位 | BBB > 0.5 | `pgp_sub` 中位 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B3DB 对照 | 444 | 63.5% | 47.1% | 32.7%（145） | 0.899 | 64.0% | 0.201 |
| Fridén | 42 | 47.6% | 50.0% | 31.0%（13） | 0.311 | 47.6% | 0.710 |
| **GKA 候选** | **787** | **3.7%（29）** | 51.8% | **3.2%（25）** | **0.014** | **10.8%** | **0.001** |

### 对照按实测值分层

| 分层 | n | SwissADME BBB+ | ADMETlab BBB 中位 | BBB > 0.5 |
|---|---:|---:|---:|---:|
| Fridén `Kp,uu ≥ 0.3` | 18 | 12 | 0.986 | 77.8% |
| Fridén `0.05–0.3` | 9 | 6 | 0.058 | 33.3% |
| Fridén `Kp,uu ≤ 0.05` | 14 | 2 | 0.002 | 21.4% |
| B3DB `control_positive` | 351 | 261 | 0.961 | 75.2% |
| B3DB `control_negative` | 93 | 21 | 0.033 | 21.5% |

⚠ 这些数字与之前临时统计时记在 CLAUDE.md 的版本有 1–2 个分子的出入
（如 B3DB 阴性侧 ADMETlab `>0.5` 从 19.3% 变成 21.5%）。
**以本表为准**——临时统计没有剔除那条樟脑记录，也没有按锚点重复做中位数合并。

### CNS MPO 里两项新拿到的量

| 集合 | `logD` 中位 | `pka_basic` 中位 |
|---|---:|---:|
| B3DB 对照 | 2.67 | 6.56 |
| Fridén | 1.72 | 5.87 |
| GKA 候选 | 2.89 | 4.65 |

---

## 留给 Step3_05 的东西

本步只把数据摆齐，以下都没做、也不该在这里做：

1. **阈值与验收**——用什么判据认定「对照集被分开了」、分离度怎么量、多少算合格。
2. **多工具分歧怎么处理**——BBB 项两个工具走向一致，**P-gp 项近乎相反**
   （SwissADME 判 GKA 候选 51.8% 是底物，ADMETlab 的 `pgp_sub` 中位只有 0.001），
   不能简单取交集或多数。
3. **CNS MPO 怎么用**——分数已出，但它在本项目对照上分辨力有限（阳 4.78 / 阴 4.48），
   且 GKA 侧对口径敏感（spread 中位 0.49）。当判据前要先回答这两点。
4. **MW 口径的选择**——平均 vs 单同位素，见上文。
5. **外推问题**——GKA 候选的 MW / TPSA 超过了阴性对照，
   流程在对照集上分得开不等于在 GKA 上同样可靠。
6. **按分子还是按骨架排序**——787 个分子只有 521 个 Murcko 骨架。
