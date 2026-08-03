# Step3：GKA 候选分子入脑预测

## 本步的定位

**这一步的目标是把「计算预测分子入脑能力」这条流程跑通。**
不是在这一轮交出「可入脑 GKA」的名单。

理由：Step1 产出的是**已知 GKA 的化学空间**，而这批分子当初是冲肝和胰腺做的
（GCK 的治疗相关表达在肝细胞与胰岛 β 细胞），入脑对糖尿病适应症不但没有价值，
还会因中枢葡萄糖感受被激活而增加低血糖风险。实测这 787 个分子
**TPSA 中位 105.25、MW 中位 462.55**，作为一个整体明显偏离 CNS 药物的性质区间。
所以**本来就不该期望有很高的入脑命中率**——命中数少是预期结果，
不是脚本出错的信号，不要靠松阈值去凑。

### 本轮范围：Step3_00 → Step3_04

**终点是一张表**：787 个 GKA 候选与入脑对照集，在同一版本工具下跑完，
各项预测值并排放进 `Step3_04_Integrated_Brain_Penetration_Results.csv`。

**⚠ `Step3_05`（候选排序与流程验收）本轮不做。**
怎么定阈值、怎么判定流程合格、怎么排序，是需要专门讨论的问题，**留到下一轮**。
本轮只负责**把数据取全、取对、可比**，不下任何判断。

这个切分是有意的：对照集数据必须先在手上，才谈得上讨论怎么用它标定阈值——
反过来先定阈值再找对照来验，就是 CLAUDE.md 里踩过两次的那个坑。

### 后续工作（均不在本轮）

* Step3_05 排序与验收规则的设计与讨论
* `Step2_02_GKA_Compound_Pool.csv`（21,488 个专利化合物）的入脑评估
* 实测脑暴露数据回填、针对性的入脑改造设计

---

## 输入

### 主输入：全部 GKA 候选

```text
Step1_Find_GKA_from_ChEMBL/Step1_GKA_Candidates_with_Properties.csv    # 787 × 42
```

**这 787 个分子就是 ChEMBL 侧的全部 GKA 候选，全部进入本步。**
它合并了两条来源（`source` 列）：

| source | 数量 | 来源 |
|---|---:|---|
| `activity` | 778 | Step1_01–06，从 GCK 激活 assay 的活性数据走下来 |
| `drug_annotation` | 5 | Step1_07，`drug_mechanism` / `-gliatin` 词干补回的零活性临床药 |
| `activity+drug_annotation` | 4 | 两条路径都命中 |

⚠ **不要用 `Step1_05_Followup_Candidates.csv`（782）**，那是合表前的中间产物，
缺 Step1_07 补回的 5 个零活性临床药（含唯一已上市的多格列艾汀 `CHEMBL4297508`），
也没有 `is_positive_control` 列。
⚠ **也不要只用 `Step1_07_GKA_from_Drug_Annotation.csv`（9）**，
那 9 个只是临床药子集，已经并在 787 里了。

核心字段：

| 字段 | 用途 |
|---|---|
| `molecule_chembl_id` | 分子唯一标识 |
| `canonical_smiles` | 模型输入结构，787 个全有、RDKit 全部可解析 |
| `molecule_pref_name` / `max_phase` | 有名字的分子在报告里用名字，结果解读时看临床阶段 |
| `parent_chembl_id` / `dedup_group` / `is_dedup_representative` / `is_salt` | 盐–游离碱关系，**沿用上游结论，不要自己重判** |
| `is_positive_control` / `control_name` | 12 个 GKA 身份对照（见下节，**不用于入脑阈值**） |
| `source` / `curated_direction` | 候选的来源与方向依据 |
| `psa` / `full_mwt` / `mw_freebase` / `alogp` / `hba` / `hbd` / `rtb` | 上游已算的理化性质，用于与本步重算值交叉校验 |
| `priority` / `potency_band` / `pactivity_median` | GKA 效力，**与入脑无关**，仅作最终解读时的另一个维度 |
| `murcko_scaffold` / `scaffold_cluster_size` | 骨架，结果按骨架看比按分子看更有意义 |

### 第二输入：入脑对照集

```text
Step3_00_BBB_Control_Set.csv        # 本步自建，见下节
```

---

## ⚠ 两套对照集，回答两个不同的问题，不能混用

这是本步最容易出错的地方。项目里有两套对照，**标定的是完全不同的轴**：

| | **GKA 身份对照**（787 里的 12 个） | **入脑对照**（本步新建） |
|---|---|---|
| 回答的问题 | 这个分子**是不是 GKA** | 这个分子**能不能进脑** |
| 来源 | Step1_05 / Step1_07 审编，`is_positive_control` 列 | 独立外部来源，实测 BBB 数据 |
| 在本轮的作用 | **只校验输入完整性**——787 里 12 个必须一个不少 | 与候选**一起跑完取数**，供 Step3_05 讨论时使用 |

**GKA 身份对照不参与入脑阈值标定。** 它们对入脑这条轴是沉默的——本项目没有它们任何一个的
实测 logBB / Kp,uu。既不能当入脑阳性对照，**也不能当入脑阴性对照**：
「它们是外周药所以应该判为不入脑」是未经验证的推断，不是事实。

⚠ 本轮**不用入脑对照集去定任何阈值**——那是 Step3_05 的事。
本轮对它的唯一要求是：**和候选分子在完全相同的条件下跑完，结果可比。**
CLAUDE.md 里「阈值必须用阳性对照标定、自检写死在脚本里」这条方法论到 Step3_05 才启用，
但对照数据必须现在就备好，否则到时候又会变成「先定阈值再找对照来验」。

---

## Step3_00：构建入脑对照集

**阳性对照 = 明确知道能入脑的小分子药物；阴性对照 = 明确知道不能入脑的小分子。**

⚠ **这些分子与 GKA、与葡萄糖激酶没有任何关系，也不需要有关系。**
它们唯一的资格是**脑暴露这件事本身已经被实测确证**。
不要试图从 787 个 GKA 候选里挑对照，也不要因为某个对照「和本项目靶点无关」
就把它换掉——正是因为无关，它才能独立地检验这条流程。

### ✅ 已完成 —— 见 `Step3_00_BBB_Control_Set/`

数据源：**ChEMBL 37 assay `CHEMBL1798466`**
= Fridén et al. *J Med Chem* 2009;52:6233（PMID 19764786），
大鼠 Sprague-Dawley 实测 **`K(p,uu,brain)`**（游离脑/游离血浆），42 个化合物。

| 分组 | 数量 |
|---|---:|
| `control_positive`（Kp,uu ≥ 0.30） | **18** |
| `control_negative`（Kp,uu ≤ 0.05） | **14** |
| `intermediate` | 9 |
| `no_value` | 1 |

Kp,uu 跨度 **0.006–2.0（333 倍）**。**8 个自检锚点方向全部正确**，写死在脚本里。

### 判据选了 Kp,uu，另一个候选是 logBB

另一个候选来源是 **B3DB**（Meng et al. *Sci Data* 2021，50 个已发表数据集的汇编，
7,807 个化合物，判据 `logBB > −1`）。两个判据测的不是同一个量：
`logBB` 是脑/血浆**总**浓度比，`Kp,uu` 是**游离**浓度比。
**本步选 Kp,uu，理由是只有游离药物能作用于靶点。**

36 个共同化合物里 23 个分组相同、**4 个分组不同**：

| 化合物 | Kp,uu | 本表 | B3DB | B3DB 有数值 logBB？ |
|---|---:|---|---|---|
| Loperamide | 0.007 | 阴性 | BBB+ | 有（−0.25） |
| Paclitaxel | 0.007 | 阴性 | BBB+ | 有（−0.55） |
| Nelfinavir | 0.019 | 阴性 | BBB+ | 有（−0.93） |
| Propranolol | 0.610 | 阳性 | BBB− | **无**（标签继承自上游文献） |

差异有两种成因，性质不同：前三个都是 P-gp 底物，**属于两个量在外排底物上的定义差异**；
Propranolol 那行 B3DB 中无数值，属于**数据可得性**问题。

⚠ **这是判据差异，不是谁对谁错。** 判据取决于问题——
若关心组织总蓄积、或需要大样本训练 ML 模型，B3DB 的规模正是所需。
差异全部保留在 `Step3_00_B3DB_Comparison.csv` 中可复核。

**这一项本身是给 Step3_05 的证据**：若某预测工具用 logBB 类数据训练，
它的行为可能更接近 B3DB 的分组——到时候要能分辨
「工具错了」还是「工具在用另一个判据」。

### ⭐ 两组匹配对

| 组 | 跨度 | 内容 |
|---|---|---|
| **β 阻滞剂** | 25× | Metoprolol 0.64 → Propranolol 0.61 → Pindolol 0.50 → Alprenolol 0.38 → Oxprenolol 0.20 → Nadolol 0.037 → **Atenolol 0.026** |
| **阿片类** | 147× | Oxycodone 1.03 → Codeine 0.89 → Oxymorphone 0.79 → Morphine 0.15 → M3G 0.011 → **Loperamide 0.007** |

β 阻滞剂那组是**一条连续梯度而不只是一对**——七个同类药按亲脂性单调排开。
预测流程若有真实分辨力，应该能重现这个次序，
这比「阳性组均值 > 阴性组均值」严格得多。

**Loperamide（Kp,uu 0.007）是最关键的一条**：它理化性质像 CNS 药
（高亲脂、TPSA 低、MW 适中），**纯理化模型必然误判为入脑**，
只有把 P-gp 外排纳入的模型才判得对——
它是「CNS MPO 单独够不够用」的判别探针。

**42 行全部进入 Step3_01–04**；二分类检验只用 32 行
（`use_for_separation_test = True`），但中间带 9 行的预测值在 Step3_05 有用，
不要提前丢掉。

---

## Step3_01：结构标准化

使用 RDKit：

1. 检查 SMILES 是否可解析，**失败的单独记录，不得静默丢弃**。
2. 去除盐和无机对离子，保留主要有机结构。
3. 生成 `standardized_smiles`。
4. 保留原始 SMILES 和处理记录（`standardization_note`）。

⚠ **787 里有 3 个多组分结构**（`canonical_smiles` 含 `.`）：

| ChEMBL ID | 说明 | 去盐后 |
|---|---|---|
| `CHEMBL1204008` | 母体 `CHEMBL575092` | 单独结构 |
| `CHEMBL4297302` | MK-0941 甲磺酸盐 | 得到游离碱，与 `CHEMBL3580737` 一致 |
| `CHEMBL5095182` | Globalagliatin 盐酸盐 | 与 `CHEMBL4297399` 一致 |

**787 行全部保留**（加列不删行）。去盐后出现相同结构时，用
上游的 `dedup_group` / `is_dedup_representative` / `parent_chembl_id` 标注关系，
**不要自己按 InChIKey 重判**——CLAUDE.md 明确记着盐与游离碱的 InChIKey 完全不同
（`KJSGTWFWVTYPFZ-…` vs `PIDNRTWDGDJKSQ-…`），药物层去重与跨库对齐是两个场景、两个键。

另注：787 个 InChIKey 全唯一，但**有 28 组分子前 14 位相同**（骨架同、立体不同）。
**不得按前 14 位去重**——GKA 的手性中心通常决定活性。

对照集分子走同一套标准化，保证与候选分子口径一致。

输出：

```text
Step3_01_Standardized_Candidates.csv
```

---

## Step3_02：基础理化性质计算与 CNS MPO

RDKit 本地可算：MW、cLogP、TPSA、HBD、HBA、Rotatable bonds、芳香环数。

⚠ **`LogD7.4` 和 `pKa` RDKit 算不出来**，而 **CNS MPO 的 6 个参数里正好有这两个**
（MW / cLogP / cLogD7.4 / TPSA / HBD / pKa_most_basic）。缺了就算不出标准 CNS MPO。

处理方式，按可得性择一，**并在输出中用 `cns_mpo_source` 列标明实际用的是哪种**：

1. ADMETlab 3.0 直接给 logD7.4 与 pKa（与 Step3_03 同一次提交，顺带拿到）；
2. 本地若有 ChemAxon `cxcalc` 或 ACD/Percepta，优先用它们；
3. 都没有则退化：用 cLogP 代 logD、pKa 按官能团规则粗判，
   **此时 CNS MPO 记为近似值并在报告中明确标注口径**。

⚠ 退化方案有已知的系统性偏差：**787 里 115 个分子含羧酸（14.6%）**，
生理 pH 下带负电，CLAUDE.md 明确警告过「拿 alogp 当 logD 用会高估膜通透性」。
入脑对照集里的 Cetirizine、Fexofenadine 也是酸，同样受影响。
**退化方案下这批分子要单独标注，不能与其余分子同等看待。**

与上游 `psa` / `full_mwt` / `alogp` 交叉校验，**差异显著的分子单独列出**，不要默默覆盖。

输出：

```text
Step3_02_Physicochemical_Properties.csv
```

---

## Step3_03：入脑和外排预测

| 工具 | 主要输出 | 可得性 | 本轮 |
|---|---|---|---|
| SwissADME | BBB、P-gp、BOILED-Egg | 网页，免费，无官方 API，单批约 200 个 | 必做 |
| ADMETlab 3.0 | BBB、P-gp、BCRP、Caco-2、PPB、logD、pKa | 网页，免费，支持 CSV 批量 | 必做 |
| QikProp | QPlogBB、CNS activity、MDCK | **需 Schrödinger license** | 有则做 |
| ACD/Percepta | LogBB、CNS Access Score | **需 license** | 有则做 |

所有工具使用 `standardized_smiles` 作为输入。

### ⚠ 前两个是免费网页工具，无官方批量 API

需要人工在浏览器里提交并下载结果，脚本不能假设「跑一下就有结果」，必须拆成两半：

* **导出**：Step3_01 之后生成分批的提交文件 `Step3_03_submission_<tool>_batch<N>.smi`，
  附提交说明（网址、勾哪些选项）。
* **导入**：给每个工具写 `parse_<tool>()`，读下载回来的结果文件，
  **校验分子数与标识对得上**（网页工具常改列名、丢行、改 ID 大小写），
  对不上直接报错而不是静默丢行。

这样 Step3_01/02 可以先本地跑通并出中间产物，Step3_03 的结果到位后直接并入，
**不必重跑上游**。

### ⚠ 787 个分子必须分批，分批就有批次效应

SwissADME 单批约 200 个，787 + 对照要分 5 批左右。网页模型会静默更新，
**跨批次的结果未必可比**。两条应对，都要做：

1. **入脑对照集完整地放进每一批**（每批只多约 20 个分子，成本极低）。
   同一个对照在不同批次给出不同结果，就是批次漂移的直接证据。
2. **记录每一批的提交时间与工具版本**，写进 `batch_id` / `submitted_at` / `tool_version` 列。

对照只放在其中一批、然后拿它的阈值去套别的批次，是本步最隐蔽的错误来源。

---

## Step3_04：结果合并

按 `molecule_chembl_id` 合并（对照集用各自的标识列）。不同工具的字段加前缀：

```text
rdkit_tpsa
chemaxon_cns_mpo
swissadme_bbb
swissadme_pgp
admetlab_bbb
admetlab_bcrp
qikprop_qplogbb
```

保留这几列，让对照与候选在同一张表里可比、可追溯：

* `set`：`candidate` / `control_positive` / `control_negative`
* `batch_id` / `submitted_at` / `tool_version`：批次追溯
* 各工具的覆盖标记：**哪个分子在哪个工具上没拿到结果，如实记录，不要填默认值**

输出：

```text
Step3_04_Integrated_Brain_Penetration_Results.csv
```

---

## Step3_05：候选排序与流程验收 —— 🚧 本轮不做，待讨论

**这一步非常重要，需要专门讨论后再定，不要在本轮顺手实现。**

已经明确、可作为讨论起点的几点：

* **必须先验收、再排序**。流程分不开对照集，任何排序都没有意义。
* `selfcheck_bbb_controls()` 要写死在脚本里（CLAUDE.md 方法论：
  Step1_05 写进脚本所以守住了，Step2_02 没写所以没守住）。
* 分不开就是流程不合格，**回去查工具和参数，不是调阈值把结果凑出来**。
* 输出遵循「加列不删行」：`keep` / `exclude_reason`，被排除的行留在 CSV 里。

**尚待讨论、本轮不预设答案的问题：**

1. 用什么判据认定「对照集被分开了」？分离度怎么量、达到多少算合格？
2. 多个工具结论不一致时怎么办——取交集、取多数、还是分别报告？
3. CNS MPO ≥ 4 这类文献阈值，在本项目的对照集上到底成不成立？
4. 787 个分子只有 521 个 Murcko 骨架，**排序应该按分子还是按骨架**？
   直接取 top-N 会拿到同一篇 SAR 论文的一串同系物（CLAUDE.md 已验证过的坑）。
5. 现有工具都是**被动扩散 + 外排**的视角，GKA 若有转运体介导的摄取，这套框架看不见。

这些问题的答案会改变脚本结构，所以**先出 Step3_04 的数据，再回来讨论**。

---

## 约定

* 不得修改原始输入。
* 不得静默删除失败分子——解析失败、预测失败的分子单独记录并说明原因。
* 一律「加列不删行」，787 行自始至终保留。
* 报告开头记录出处：输入文件与行数、各外部工具的版本与访问日期、运行时间。
* 报告中必须写明：**本轮的产物是数据，不是结论**——
  既不下「流程是否可用」的判定，也不下「哪些 GKA 能入脑」的判定，两者都属 Step3_05。
