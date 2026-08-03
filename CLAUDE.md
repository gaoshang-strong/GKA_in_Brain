# GKA_in_Brain

从 ChEMBL 中挖掘葡萄糖激酶激活剂（GKA, glucokinase activator）。

## 数据

ChEMBL 37 SQLite，本地路径：

```
ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db     # 29 GB，不入 git
```

不入 git 的还有官方压缩包与 EBI 随附文件，重新下载与 SHA256 校验方式见
`ChEMBL/ChEMBL_37/README_ChEMBL_37.md`。数据库结构与字段含义速查见
`ChEMBL/chembl_37_profile_report.md`（由 `ChEMBL/chembl_profile.py` 生成）。

**一律以只读模式打开**：`sqlite3.connect(f"file:{db}?mode=ro", uri=True)`。

SureChEMBL 2.0 专利化学，`2026-07-17` 全量快照，本地路径：

```
SureChEMBL/SureChEMBL_2026-07-17/     # 17 GB，9 个文件，不入 git
```

**每两周发一版并覆盖 `latest/`，所以固定在日期目录上**，任何结论都对应这一版。
下载/校验/表结构见 `SureChEMBL/SureChEMBL_2026-07-17/README_SureChEMBL_2026-07-17.md`，
库的结构与内容速查见 `SureChEMBL/surechembl_2026-07-17_profile_report.md`
（由 `SureChEMBL/surechembl_profile.py` 生成）。

**大表绝不能整表读进内存**（`patent_compound_map` 15.4 亿行），用 duckdb 直查 parquet：

```python
con = duckdb.connect()
con.execute("SET threads=8; SET memory_limit='12GB'; SET enable_progress_bar=false;")
```

Step3 用的两份外部数据，**都放在 `Step3_00_BBB_Control_Set/`，各带一个 README 记出处与 SHA256**：

```
Friden2009_CHEMBL1798466_raw.tsv    42×36   ← 入脑对照集的**唯一**数据源
B3DB_classification.tsv / _regression.tsv   ← **不是**对照集来源，只用于记录判据差异
```

- **Fridén 2009**（J Med Chem 52:6233, PMID 19764786）原文**非开放获取**，
  ChEMBL assay `CHEMBL1798466` 是这批 `K(p,uu,brain)` 数据的公开途径。
  快照冻结后**无需 29 GB 的 ChEMBL 库也能复现**（脚本 `--db` 不存在时自动回落，两条路径产物逐字节一致）。
- **B3DB**（Meng et al. *Sci Data* 2021）是 **50 个已发表数据集的汇编**，
  `2026-08-02` 快照。上游会不定期追加数据，`main` 不是固定版本，结论对应这一版。

## 检索链路

```
Step1_01  UniProt P35557 → component_sequences → target_components → target_dictionary
Step1_02  tid → assays → 汇总实验类型/物种/可信度/来源/活性规模
Step1_03  assay → 分类（激活/抑制/结合/GKRP 相互作用/细胞表型/无法判断）
Step1_04  激活 assay → activity → 分子，分别给出效力/效能/证据（不合成总分、不排序）
Step1_05  分子层方向判定 → 排除打标 → 效力单轴分档排序 → 骨架去冗余 → 候选表
Step1_06  候选分子 → molecule_dictionary/compound_structures/compound_properties
          /molecule_hierarchy/molecule_synonyms → 理化性质主表
Step1_07  drug_mechanism / usan_stems / molecule_synonyms → 补回零活性的临床 GKA
Step1_整合 两条路径合表 → Step1_GKA_Candidates_with_Properties.csv（787 × 42）

Step2_01  标题正则 ∪ 实体锚定 → 同族扩展 → GCK 相关专利
Step2_02  专利分层 L1/L2/L3 → 权要化合物 → 噪声指标 → SureChEMBL 侧候选池

Step3_00  B3DB 分类对照 + Fridén 定量参考 → 入脑对照集（与 GKA 无关）
Step3_01  787 GKA 候选 ∪ 487 对照 → 统一输入表 → RDKit 标准化 + 描述符
Step3_02  → SwissADME 提交文件（7 批 × ≤200）
Step3_03  → ADMETlab 提交文件（16 批 × ≤98）
Step3_04  两个工具的结果 → 逐行结构校验 → 锚点合并 → 整合表（1,274 × 251）
Step3_05  排序与流程验收（**未做，需专门讨论**）
```

**Step3 是把「计算预测入脑」这条流程跑通，不是这一轮交出可入脑 GKA 名单。**
终点是一张表（候选与对照在同条件下跑完的预测值），**Step3_05 的阈值与验收标准留待讨论**——
先出数据再定判据，反过来就是踩过两次的那个坑。

**两个库各自独立检索，最后 Union。** ChEMBL 侧的分子**不进** SureChEMBL 侧的检索式，
只做事后验证（`val_*` 列）与最终比较——否则「SureChEMBL 独立贡献了多少」是循环论证。
同理，以 ChEMBL 分子为种子的**结构相似性检索（`fpsim2`）属于扩展臂，不是独立锚点**。

已锚定的 GCK 靶点（ChEMBL 37）：

| target | tid | type | assays | activities |
|---|---|---|---|---|
| `CHEMBL3820` Hexokinase-4 | 20095 | SINGLE PROTEIN | 227 | 3,222 |
| `CHEMBL3885579` Glucokinase/GKRP | 117123 | PROTEIN-PROTEIN INTERACTION | 1 | 40 |

各步骤的收敛规模（ChEMBL 37）：

| 步骤 | 产物 | 规模 |
|---|---|---|
| Step1_02 | GCK assay 清单 | 228 |
| Step1_03 | 其中「GCK 激活」且靶点身份无疑 | 142 |
| Step1_04 | 激活 assay 里出现过的分子 | 1,333（**不等于激活剂**） |
| Step1_05 | 过方向门与删失门后入选 | 1,222；进 5 个效力档 1,007 |
| Step1_05 | 后续实验清单 `Step1_05_Followup_Candidates.csv` | 782（P1 39 / P2 83 / P3 202 / P4 458） |
| Step1_06 | 理化性质主表 | 782 行 × 75 列，全部命中无缺失 |
| Step1_整合 | 加上 Step1_07 的 5 个零活性临床药 | **787 × 42**，`source` 分 activity / drug_annotation |

**ChEMBL 侧下游只读 `Step1_GKA_Candidates_with_Properties.csv`**（顶层）——
它带着 `priority` / `potency_band` / `pactivity_median` / `direction` / `curated_direction` /
骨架列 + 全部理化性质。

SureChEMBL 侧（`2026-07-17` 快照）：

| 步骤 | 产物 | 规模 |
|---|---|---|
| Step2_01 | GCK 相关专利（召回优先） | 35,793 篇 / **7,620 同族** |
| Step2_02 | L1 标题写 activator | 638 篇 / **143 同族** |
| Step2_02 | L2 标题/权要提靶点，或说明书提 ≥3 次，且权要有化合物 | 4,283 篇 / **1,388 同族** |
| Step2_02 | L1+L2 权要化合物 → 过噪声门 | 30,627 → **keep 21,488** |

**SureChEMBL 侧下游只读 `Step2_02_GKA_Compound_Pool.csv`**（`keep = TRUE`）。

Step3 侧（入脑预测）：

| 步骤 | 产物 | 规模 | 状态 |
|---|---|---|---|
| Step3_00 | 入脑对照集 | **487** = 351 BBB+ / 94 BBB− / 42 Fridén | ✅ |
| Step3_01 | 统一输入 + RDKit 标准化 | **1,274 × 67**（787 GKA + 445 B3DB + 42 Fridén） | ✅ |
| Step3_02 | SwissADME 提交 | 7 批，1,392 条（1,272 唯一结构 + 锚点重复 120） | ✅ |
| Step3_02 | SwissADME **结果** | 返回 1,391 条 → 合并 **1,271 个唯一分子** | ✅ 7/7 批 |
| Step3_03 | ADMETlab 提交 | 16 批，1,571 条（1,271 唯一结构 + 锚点重复 300） | ✅ |
| Step3_03 | ADMETlab **结果** | 1,571 条 → **1,271 个唯一分子**（每批 123 列） | ✅ 16/16 批 |
| Step3_04 | 结果合并 + CNS MPO | **1,274 × 261**，两个工具各覆盖 1,273 行，MPO 出分 1,273 | ✅ |
| Step3_05 | 排序与验收 | 只出了 4 张**描述性**图 + 说明，**判据仍未做** | ⬜ **需专门讨论** |

Step3 的两份讲解文档：

- `Step3_Methods_Explained_From_Scratch.md`（Step3 根目录）——
  RDKit / SwissADME / ADMETlab / CNS MPO 的原理，**面向零化学背景**，
  含血脑屏障、LogP vs LogD、TPSA、pKa 的基础解释与本项目踩过的坑。
- `Step3_05_Ranking_and_Validation/Step3_05_Figures_Explained.md`——
  四张图（工具比例 / 概率分布 / 化学空间 / MPO 六项拆解）的解读，
  **图里没有任何阈值线**，判定仍留给 Step3_05。

**Step3 下游只读 `Step3_04_Integrated_Brain_Penetration_Results.csv`**（1,274 × 251）。
以下数字全部出自它，**以此为准**——Step3_04 之前那次临时统计没剔除樟脑那条、
也没按锚点重复做中位数合并，个别格子差 1–2 个分子。

两个工具并排（**只是数据，判定属 Step3_05**）：

| 集合 | n | SwissADME BBB+ | Pgp 底物 | **BBB+ 且非 Pgp** | ADMETlab BBB 中位 | BBB > 0.5 | `pgp_sub` 中位 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B3DB 对照 | 444 | 63.5% | 47.1% | 32.7%（145） | 0.899 | 64.0% | 0.201 |
| Fridén | 42 | 47.6% | 50.0% | 31.0%（13） | 0.311 | 47.6% | 0.710 |
| **GKA 候选** | **787** | **3.7%**（29 个） | 51.8% | **3.2%**（25 个） | **0.014** | **10.8%** | **0.001** |

对照按实测值分层，两个工具走向一致：

| | SwissADME（BBB+ 个数） | ADMETlab（BBB 中位 / >0.5 占比） |
|---|---|---|
| Fridén `Kp,uu ≥ 0.3` | 18 个里 12 个 | 0.986 / 77.8% |
| Fridén `0.05–0.3` | 9 个里 6 个 | 0.058 / 33.3% |
| Fridén `Kp,uu ≤ 0.05` | 14 个里 12 个判 BBB− | 0.002 / 21.4% |
| B3DB `control_positive` | 351 个里 261 个 Yes | 0.961 / 75.2% |
| B3DB `control_negative` | 93 个里 72 个 No | 0.033 / 21.5% |

CNS MPO 已出分（1,273 行；拐点取自本地 PDF，见下）：

| 集合 | 总分中位 | ≥ 4 | ≥ 5 |
|---|---:|---:|---:|
| B3DB 对照 | 4.74 | 75.3% | 41.1% |
| Fridén | 4.80 | 69.0% | 33.3% |
| **GKA 候选** | **3.78** | **38.5%** | **6.9%** |

六项 T0 中位显示**差距只集中在两项**：MW（对照 1.00 / GKA **0.27**）与
TPSA（1.00 / **0.49**）；LogP、pKa 三组都是 1.00，HBD 上 GKA 反而更好（0.83）。
`logD` 中位：B3DB 2.67 / Fridén 1.72 / GKA 2.89；
`pka_basic` 中位：6.56 / 5.87 / **4.65**。

对照集的选取规则（Step3_00）：
BBB+ 取 B3DB regression `group B` 且 `logBB ≥ −0.5`；
BBB− 取 `group B`，不足再用 `group A` 补，`logBB ≤ −1.1`；
两侧同一 Murcko 骨架最多留 3 个，簇内按「离 −1 边界远 + 来源多 + 靠近 GKA 空间」综合排序。
**不设固定数量**，合格多少要多少。

**⚠ 787 个 GKA 候选只做骨架分组，不做删除式去冗余**——1,274 行进、1,274 行出，
去盐后结构相同的只打 `dup_group` / `is_dup_representative`，行全部保留。

⚠ 两侧重叠只有 **313 / 21,488（1.5%）**——两个库在看几乎不相交的化学空间。
但专利侧那 21,175 个**没有任何活性数据**，「是不是 GKA」尚未证实。

## 关键事实与陷阱

这些都是在真实数据上验证过的，不要凭直觉推翻：

- **必须按 UniProt accession 锚定，不能按名字。** GCK 在 ChEMBL 里的 `pref_name` 是
  **Hexokinase-4**；搜 "Glucokinase" 会漏掉主靶点，且会搜到 `CHEMBL1075152`——
  那是 GKRP（基因 GCKR），不是 GCK。
- **⚠ 挂在靶点下的 assay 未必真在测该靶点。** ChEMBL 37 把 **52 个 MAP4K2 /
  Germinal Center Kinase 的 assay 误挂在了 `CHEMBL3820` 下**——"GCK" 是歧义缩写
  （glucokinase 465 aa / MAP4K2 820 aa，后者是蛋白激酶）。
  判别规则：**描述只用缩写 "GCK" 而不写 "glucokinase" 的，全部是 MAP4K2**。
  下游筛 GKA 必须加 `target_identity_suspect == FALSE`。
  完整证据见 `Step1_Find_GKA_from_ChEMBL/Step1_03_GCK_Assay_Classification/
  Step1_03_Target_Mismapping_MAP4K2.md`。
  **靶点身份校验要作为独立一步**，方向分类（激活/抑制）回答不了这个问题。
  序列是最可靠的裁判：`component_sequences.sequence` 里可直接比对描述中的肽段与残基范围。
- **PPI 靶点不能与主靶点无脑合并**：`CHEMBL3885579` 测的是 GCK–GKRP 相互作用
  （GKRP 解离剂机制），与直接激活 GCK 是两类机制。
- **⚠ assay 有目的性，但 assay 里的分子未必都符合那个目的。**
  assay 与分子是一对多：142 个激活 assay → 1,333 个分子，中位 22 个/assay，
  最多 70 个（`CHEMBL1825593`）。assay 描述定义的是「测什么量、相对什么基线」，
  **不含任何分子的方向结论**。实测：24 个含比值读数的激活 assay 中，
  **11 个内部同时存在激活与降低的结果**——`CHEMBL1825590` 的 70 个分子里，
  52 个 >1.05（激活）、6 个 <0.95（酶活反而降低）。
  按 Step1_04 的数据估算，「激活 assay 里的 1,333 个分子」中约 **170 个没有正向
  激活证据、63 个方向冲突**。
  **当前采用的策略：以「assay 目的是测激活」作为筛选门槛。** 这是保险的一侧——
  宁可混进非激活分子，也不漏掉真 GKA。**分子层的方向判定尚未做**，
  下游把这批分子当 GKA 候选时必须知道这个口径，
  `Step1_04` 的分子数是「出现在激活 assay 里的分子」，不是「激活剂数」。
- **反方向的分子未必看得见**：单向量纲（`EC50`/`AC50`）装不下抑制——
  无效化合物不会得到负值，只会变成删失 `>` 或根本不录入
  （142 个激活 assay 里有 86 条 `>` EC50，分布在 30 个 assay）。
  **看不见不等于不存在**，删失值不能当等号值汇总，否则把弱化合物算强。
  双向量纲（`FC`/`Ratio`/`%max`/`Emax`）才会把两个方向都记下来。
- **比值型读数的「无效基线」不统一，必须逐 assay 从描述读，不能设全局阈值。**
  同一个 `0.69` 在不同 assay 里含义相反：
  - 分母是**未处理对照**（17 个 assay）→ 基线 1，`0.69` 是酶活降到 69%
  - 分母是**参比激活剂**（7 个 assay，如 `CHEMBL1825591/93` 的
    "ratio of ... in compound treated to **Ro-28-1675**"）→ 基线 0，
    `0.69` 是「达到参比效果的 69%」，是个不错的激活剂
  - percent 类同一个 `%max` 里两种基线都有：`CHEMBL2353000` min=0（基线 0）、
    `CHEMBL4014755` min=452（基线 100）
  Step1_05 已把 35 个 (assay, type, scale) 组逐组解析完，结论是
  `control` 21 组 / `reference` 8 组 / `kinetic` 3 组 / `max` 2 组 / `uncertain` 1 组，
  全表在 `Step1_05_..._Candidate_Selection.md` 里可直接复用，别再从头判一遍。
  按「fold ≤ 1.05 即无激活」粗判会误杀 8 个 reference 组的整组分子。
  还有极性问题：`S0.5`/`S50`/`Km` 是**值越小激活越强**，与 `FC` 方向相反；
  它们测的是酶对葡萄糖的半饱和常数（GKA 把它从野生型 ~7 mM 拉到 0.5–2 mM），
  是机制读数，**不是化合物效力**，当效力用会得出「EC50 = 0.6 mM」这种量级错误。
- **⚠ 整条链路挂在 `activities` 上，零活性记录的分子在结构上不可见。**
  `靶点 → assays → activity → 分子` 这条路走不到「有 ChEMBL ID 但没有活性数据」的分子。
  **ChEMBL 的分子有两条独立入口**：从文献活性数据抽取，和从药名/临床登记册收录
  （`compound_records.doc_type = DATASET`，来源 `USAN` / `INN` / `ATC` / `CANDIDATES`）。
  实测：**6 个带 `-gliatin` 词干的临床 GKA 里有 4 个不在 Step1_06 的 782 个候选中**，
  包括唯一的 III 期药**多格列艾汀**（`CHEMBL4297508`，中国已上市，`activities = 0`）。
  第三种漏法：`CHEMBL4297399` LY-2608204 有 6 条活性，但全来自 SARS_COV_2 / HDAC6
  筛选数据集，**没有一条打在 GCK 上**。
  补救见 `Step1_07`，两条独立锚定路径：
  - **`drug_mechanism`**：挂 `CHEMBL3820` 的 6 条，`action_type` 全是 `ACTIVATOR`、
    `mechanism_of_action` 全是 "Hexokinase type IV activator"。
    **这是人工审编的方向标注**——Step1_03/05 用规则 + LLM 从 assay 描述推的那个结论，
    这里现成且更可靠。**做方向判定前先查这张表。**
  - **`usan_stems`**：`-gliatin` 的官方注释就是 `glucokinase activator`，
    按药名后缀就能认出 GKA。`molecule_dictionary.usan_stem` 有的为空，
    要用 `molecule_synonyms LIKE '%gliatin%'` 兜底。
- **⚠ 同一个药以「游离碱 + 盐」两个 chembl_id 存在，而 InChIKey 归并不了它们，
  `molecule_hierarchy.parent_molregno` 才行。** 实测两对：
  - MK-0941：`CHEMBL3580737` 游离碱（有活性）+ `CHEMBL4297302` **甲磺酸盐**
    （SMILES 带 `.CS(=O)(=O)O`，0 活性），`parent_molregno` → 游离碱
  - Globalagliatin：`CHEMBL4297399` LY-2608204 + `CHEMBL5095182` **盐酸盐**
    （SMILES 带 `Cl.`），`parent_molregno` → 前者

  **盐与游离碱的 InChIKey 完全不同**（`KJSGTWFWVTYPFZ-…` vs `PIDNRTWDGDJKSQ-…`），
  所以**按 InChIKey 去重会把同一个药算成两个**。
  **药物层去重用 `parent_molregno`；跨库对齐（如与 SureChEMBL）才用 InChIKey，
  且要先归到母体再取 key。** 两个场景用不同的键，别混。
- **⚠ 任何依赖「被测了多少次 / 测了几个轴」的量，都不能进排序或分层——
  量到的会是数据可得性，不是分子质量。** 这条在 Step1_05 上栽过两次：
  - **证据强度**（重复 assay 数、独立文献数）：三个 phase 2 药 MK-0941、AZD-1656、
    PF-04991532 的证据都是「弱」（各只有 1 个 assay、1 篇文献），
    而参比化合物 Ro-28-1675 是「强」（15 个 assay、7 篇文献）——
    因为它是大家用来做对照的，不是因为它更好。折进排序会系统性地把最像药的分子排后面。
  - **效能佐证**（有没有 fold/% activation 读数）：**单 assay 的分子只有 14%
    有效能记录，多 assay 的有 67%**。初版把它当分层条件（A 层 = 效力 + 效能双证），
    结果 pAct ≥ 6.5 的 548 个分子里 343 个仅因 ChEMBL 没测效能就掉层，
    A/B 两层效力值域几乎完全重叠（A 6.50–8.78 / B 6.00–8.72），
    223 个 B 层分子的效力 ≥ A 层中位数。
  正确做法：**分档只用效力单轴，其余作独立标记列并列展示**。
  改完后 `has_efficacy_corroboration` 在前四个效力档稳定在 29–32%，
  与效力无关——这个「平」本身就是它不该进分档的证据。
- **阈值要用阳性对照标定，不能按整体分布拍。** Step1_05 初版按 pActivity 分布
  取 7.0 / 6.5，自检只过 4/6：Ro-28-1675 跨实验 EC50 是 127–690 nM
  （pAct 中位 6.39）、Piraglitin 是 364–6320 nM（6.145）——
  **已知临床 GKA 本来就在几百 nM 这一档**，阈值把临床药从中间劈开了。
  **阳性对照集（11 个，Step1_07 后的完整版）**，每步筛选都该拿它们自检，
  落不进候选说明规则错了，不是数据错了：

  | ChEMBL ID | 名称 | phase | 备注 |
  |---|---|---|---|
  | `CHEMBL4297508` | Dorzagliatin（多格列艾汀/Sinogliatin） | 3 | 中国已上市，**0 活性** |
  | `CHEMBL1783734` | Piraglitin | 2 | |
  | `CHEMBL2165615` | Neriglitin | 2 | |
  | `CHEMBL2165620` | PF-04991532 | 2 | |
  | `CHEMBL3219124` | AZD-1656 | 2 | |
  | `CHEMBL3580737` | MK-0941 游离碱 | 2 | 与下一条是同一个药 |
  | `CHEMBL4297302` | MK-0941 药物条目 | 2 | **0 活性** |
  | `CHEMBL4297399` | LY-2608204 / Globalagliatin | 2 | 活性全在别的靶点上 |
  | `CHEMBL5095182` | Globalagliatin 盐酸盐 | 1 | 与上一条同药，**0 活性** |
  | `CHEMBL5095262` | Cadisegliatin | 2 | **0 活性** |
  | `CHEMBL5072532` | BMS-820132 | — | 在 782 里 |

  另有 `CHEMBL1096435` Ro-28-1675——是文献常用**参比化合物**不是临床药，
  证据等级高但效力只有中位水平，适合当「规则会不会把参比化合物排下去」的探针。
- **收窄候选靠骨架，不靠阈值**：1,333 个分子只有 521 个 Murcko 骨架，
  最大一簇 46 个；效力 ≥7.0 的 264 个分子只归属 122 个骨架。
  直接按效力取 top-N 会拿到同一篇 SAR 论文的一串同系物——看着 50 个，其实 2 个化学起点。
- **`assay_type` 没有区分力**：226/228 都是 `B`(Binding)，包括描述明写
  "Activation of glucokinase" 的。不能用它判方向。
- **`bao_format` 不可信**：多个 E. coli 重组纯酶实验被标成 "tissue-based format"
  （如 `CHEMBL5048855`）。不能用它判细胞实验。
- **`confidence_score` 不是数据质量分**，而是靶点归因的粒度 + 直接/同源。
  9 = 直接指认到单一蛋白，8 = 同源映射。GCK 数据 222/228 为 9，很干净。
- **文本陷阱**：`expressed in CHO/Sf9 cells` 是表达宿主不是细胞实验；
  `enzyme affinity for glucose` 是底物亲和力不是化合物结合；
  `filter-binding assay` 里的 binding 只是方法名；
  `in absence of GKRP` 是实验条件不是在测 GKRP 相互作用。
- **`standard_type` 有大小写变体**（`Inhibition` vs `INHIBITION` 等），
  `GROUP BY` 会拆开，先归一。
- **只有部分 `standard_type` 有官方标准化规则**（见 `activity_stds_lookup`）。
  查不到的（`Activity`、`Ratio`、`Percent Effect`、`Z score`）是自由文本，
  含义随实验而异，**不能跨实验汇总**。GKA 数据里主力可定量指标是 `EC50`（1,378 条）。
- **实验条件没有结构化**：`assay_parameters` 里只有 EFO 疾病本体标注。
  葡萄糖浓度、孵育时间只写在 `assay_description` 自由文本里——
  低糖/高糖条件是区分 GKA 类型的关键，需要时得自己解析描述。

### 理化性质侧（Step1_06 验证）

- **ChEMBL 37 的 `compound_properties` 只有 15 列**，没有 `cx_logp` / `cx_logd` /
  `cx_most_apka` / `cx_most_bpka` / `molecular_species`。别按老版本的字段名去查。
  有的是：`mw_freebase`、`full_mwt`、`alogp`、`hba`、`hbd`、`psa`、`rtb`、`ro3_pass`、
  `num_ro5_violations`、`aromatic_rings`、`heavy_atoms`、`qed_weighted`、
  `full_molformula`、`np_likeness_score`。782 个候选上**一列不缺**。
  **`alogp` 是 logP 不是 logD**，782 个里 95 个含羧酸，生理 pH 带负电，
  拿 alogp 当 logD 用会高估膜通透性。
- **性质分布**（782 个候选）：PSA 中位 105.25、MW 中位 462.55、ALogP 中位 4.14、
  QED 中位 0.47、Lipinski 0 违规 467 个。这批分子是冲**肝和胰腺**做的。
- **别漏掉这两张表**，它们也是分子性质，且都有数据：
  - `compound_structural_alerts` → `structural_alerts` → `structural_alert_sets`：
    284/782 分子命中，473 条。MLSMR 268 / Dundee 200 / Glaxo 3 / BMS 2。
    最常见 Hetero_hetero 73、长脂链 47、Michael acceptor 59。
  - `ligand_eff`（**按 `activity_id` 关联，不是 molregno**）：726/782 分子、1,104 条
    LE/BEI/SEI/LLE。LLE 中位 2.81。**逐 activity 计算**，同一分子在不同葡萄糖浓度下
    EC50 不同、LLE 也跟着变，不能跨 assay 平均。
- **`molecule_dictionary` 的注解字段在这批分子上基本是空的**：
  `chirality` / `prodrug` / `first_in_class` / `inorganic_flag` 有 777/782 是 **`-1`
  （未标注，不是 0）**；`oral` / `parenteral` / `topical` / `therapeutic_flag` /
  `withdrawn_flag` 全是 0。`molecule_synonyms` 只覆盖 5 个分子。
  这是正常的——绝大多数是文献化合物，没进过开发流程。
- **回溯原文用 `compound_records.compound_key`**（论文/专利里的化合物编号），
  782 个全覆盖，比 ChEMBL ID 好对号。
- **结构标识**：`canonical_smiles` / `standard_inchi` / `standard_inchi_key` 782 全有，
  RDKit 全部可解析，**InChIKey 782 个全唯一**，可直接当跨库主键。两个坑：
  - **28 组分子的 InChIKey 前 14 位相同**（骨架层一致、立体层不同），
    按前 14 位去重会把立体异构体合并——GKA 的手性中心通常决定活性，别合。
    全表 282 个分子带 `@` 标记。
  - 782 个里只有 `CHEMBL1204008` 是多组分盐（SMILES 带 `.`，母体 `CHEMBL575092`）。
    它的 `full_mwt` 含盐、`mw_freebase` 不含，取值时别混。
- `compound_structures.molfile`（2D 坐标块）782 个都有但 Step1_06 没取——
  结构已由 SMILES/InChI 完整表达，需要时再补。
  `drug_mechanism` / `drug_indication` 各只覆盖 3 个分子，
  `drug_warning` / `formulations` / `molecule_atc_classification` 为 0。

### SureChEMBL 侧（Step2_01 / Step2_02 验证）

- **数据模型是两条互不相连的标注链，只在 `patents.id` 上碰面：**

  ```
  文本侧  biomedical_entities ──< biomedical_locations >── patents
  化学侧  compounds ──────────< patent_compound_map >──── patents
  ```

  **`compounds` 与 `biomedical_entities` 之间没有任何关联。** SureChEMBL 从不记录
  「这个化合物作用于这个蛋白」，只记录「谁出现在哪篇专利的哪个部分」。
  所以**方向（激活/抑制）在 bulk 数据里判不了**——`Mechanism` 类型那 8,202 条实体
  全是工业化学词（防锈剂、防霉剂、消泡剂），`resolved_form` 全空。
- **⚠ 标注管道会整篇缺失，这是单一锚点必然漏的根因。**
  实测四篇专利标题明写 "GLUCOKINASE ACTIVATOR"、化学侧有 72–238 个化合物，
  但 `biomedical_locations` **一条记录都没有**：
  `EP-4725482-A1`（GKA 用于认知障碍与神经退行，**本项目最该找到的一篇**）、
  `US-20260200881-A1`、`CN-118453592-A`、`US-12064416-B2`。
  只用实体锚定会把它们整个漏掉。**必须多锚点并集**：
  `patents.title` 正则（不经过标注管道）∪ 实体锚定，再用 `family_id` 展开
  （实测补 6,422 篇）。改完后权要含已知 GKA 的专利召回率 38% → **95.1%**。
- **⚠ `field_id` 决定语义，全库差 6.5 倍**（`desc` 12.18 亿关联 vs `clms` 1.87 亿）。
  说明书含背景技术、会大段引用**他人的**化合物；真正主张保护的只在权利要求。
  实测 glucokinase 在 `desc` 命中 28,272 篇、`clms` 只有 2,597 篇。
  锚定步骤两个都留（`hit_*` 列），抽分子时只取 `field_id = 2`。
- **⚠ 典型的化合物专利在权利要求里根本不提靶点。**
  标题是纯化学名、权要写通式与取代基，靶点只写在说明书。实测被漏的：
  `Therapeutic agents`（权要 256 个化合物、装着 52 个 Step1 已知分子，权要提靶点 **0 次**）、
  `2-(3,5-DISUBSTITUTEDPHENYL)PYRIMIDIN-4(3H)-ONE DERIVATIVES`（说明书提 37 次）。
  **筛专利不能只看标题和权要，要加「说明书提及 ≥3 次」这一支**——
  阈值 3 是量出来的：提 1–2 次的有 23,787 篇（顺带列举靶点），≥3 次的只有 2,364 篇。
- **⚠ 专利的层级不能继承给它的化合物**，与 Step1_04 那条同构。
  实测 `EP-4725482-A1` 权要里只有 3 个化合物，**其中 2 个是葡萄糖**（开链式 + 环式）——
  权要写「激活葡萄糖激酶」必然提到底物，抽取管道就把它注册成化合物了。
  化合物层要独立判别，最有效的指标是 **`n_global`（该结构在全库出现的专利数）**：
  水 `O` 出现 1,098 万篇，多格列艾汀 187 篇，**差 4 个数量级**。
- **⚠ 缩写歧义（`GCK` / `GCKs` / `GLK` / `GK`），是 ChEMBL 侧 MAP4K2 那个坑的翻版：**

  | 写法 | `corrected_text` | `resolved_form` | 收不收 |
  |---|---|---|---|
  | `GCK` | `GCK` | （未解析） | ❌ |
  | `GCKs` | **`germinal center kinase`** | （未解析） | ❌ **就是 MAP4K2** |
  | `GcK`（小写 c） | `glucokinase` | `HGNC:4195` | ✅ |
  | `GLK` | `GLK` | **`HGNC:6865`（MAP4K3）** | ❌ |
  | `GK` | `glucokinase` | `HGNC:4195` | ⚠ 收但打标 |

  `GK` 在糖尿病文献里更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型），1,343 篇，
  与本领域高度重叠，必须标 `risk_flags`。
  `GLK activator` 是 AstraZeneca 系列的写法（47 篇/10 同族），
  **要加进标题正则必须带上下文**（`GLK\s+activator`），单用会撞 MAP4K3 与植物 Golden2-like 基因。
- **术语提醒**：文档里写「提 GCK」指的是**提到葡萄糖激酶这个靶点**（35 个归一写法的任一），
  **不是字面的三个字母 `GCK`**——后者恰恰被排除。写文档时避免这个简写。
- **⚠ `family_id = -1` 是「未分配同族」哨兵**（全库 71,862 篇），
  `COUNT(DISTINCT family_id)` 必须 `FILTER (WHERE family_id > 0)`，
  否则这 7 万篇会被算成同一个发明。
  另：`family_id` 为空的 1,172,063 篇与 `publication_date` 为空的**完全是同一批**。
- **⚠ JP / CN 的标注基本失效**：JPO 不提供全文（只有著录项+英文标题摘要），
  CNIPA 只有英文机翻全文。实测命中率 JP 0.001%、CN 0.004%，US/EP/WO 都在 0.14–0.17%。
  **这批数据实质上是 US / EP / WO 的视图**，不能说「中国/日本没有 GKA 专利」——
  是看不见，不是没有。
- **其他**：`inchi_key` 在 SureChEMBL 里**不唯一**（3,099 万行 / 2,987 万唯一），
  join 前必须 `DISTINCT`；`biomedical_locations` 里有 96 个孤儿 `patent_id`
  （`patents` 表中不存在）；本地实测行数与官方宣传数字对不上
  （实测 4,491 万专利 / 3,099 万化合物，官方称 1.166 亿 / 4,770 万），
  **做规模陈述以本地实测为准**。

### Step3 侧（入脑预测，Step3_00–03 验证）

- **⚠⚠ 两套对照集回答两个不同的问题，绝不能混用。**

  | | GKA 身份对照（那 12 个） | 入脑对照（Step3_00 建的 487 个） |
  |---|---|---|
  | 回答 | 这个分子**是不是 GKA** | 这个分子**能不能进脑** |
  | 用途 | 校验 Step1 检索有没有漏药 | 校验 Step3 预测流程有没有分辨力 |

  **那 12 个 GKA 对照不参与 Step3 的阈值标定**——本项目没有它们任何一个的实测脑暴露数据，
  既不能当入脑阳性对照，**也不能当入脑阴性对照**。
  「它们是外周药所以应该不入脑」是未经验证的推断，不是事实。
  入脑对照**与 GKA、与葡萄糖激酶无关，也不需要有关**——正因为无关才能独立检验流程。
- **`Kp,uu` 与 `logBB` 测的不是同一个量，判据选择会实质改变对照集构成。**
  `logBB` 是脑/血浆**总**浓度比（含组织结合），`Kp,uu` 是**游离**浓度比。
  Step3 选 Kp,uu，理由是只有游离药物能作用于靶点。
  36 个共同化合物里 4 个分组不同，**三个 BBB+ 冲突全是 P-gp 底物**
  （Loperamide 0.007 / Paclitaxel 0.007 / Nelfinavir 0.019），
  另一个 Propranolol 在 B3DB 中根本无数值、标签继承自上游文献。
  ⚠ **这是判据差异，不是谁对谁错**；换个问题（组织总蓄积、训 ML 模型）logBB 才是对的。
- **B3DB 的结构实测**：7,807 个分类条目里**只有 1,058（13.6%）有数值 `logBB`**，
  其余仅有继承自上游文献的分类标签。`group` 列（regression 表）按 reference 数量实测是
  A=恰好 1 个来源 / B=2–33 / C=恰好 2 / D=3–35；**官方文档未逐个定义 A/B/C/D**，
  原文在 Nature 登录墙，这是从数据量出来的理解，不是官方口径。
- **⚠ RDKit 的 TPSA 默认不计 S/P 的极性贡献，SwissADME 计入。**
  实测 `Descriptors.TPSA(m)` 只有 106/199 与 SwissADME 一致，
  `includeSandP=True` 则 **193/199** 一致。GKA 候选含 73 个磺酰胺 + 大量砜，
  **TPSA 中位从 105.2 变成 122.0**，「TPSA ≤ 90」的个数从 181 掉到 85。
  **BOILED-Egg 用的是含 S/P 那个**，CNS MPO 的 TPSA 项选错口径会直接改分数。
  同理旋转键：RDKit 默认 Strict 只对上 80/199，`NonStrict` 才 **199/199**。
  Step3_01 两种口径都留（`tpsa` / `tpsa_sandp`、`rtb` / `rtb_nonstrict`）。
- **⚠ 送外部工具前必须中性化，这不是可选项。** SwissADME FAQ 原文：
  "Is it preferable to input the neutral form of the molecule? **Yes, definitively.**
  The SMILES entry is taken as given and **not neutralized**"——模型基本在中性化合物上训练。
  Step3_01 的流程是 `Cleanup → LargestFragment(去盐) → Uncharger(中性化)`。
  1,274 个里去盐 8 个、中性化 8 个；**4 个季铵/肟盐是永久正电荷，化学上无法中性化**，已标出。
- **CNS MPO 是 Pfizer 的评分函数，RDKit 里没有任何实现。**
  原文（2026-08-03 在 PubMed 逐条核对过题名与卷期页）：
  - Wager TT, Hou X, Verhoest PR, Villalobos A. *"Moving beyond Rules: The Development of a
    Central Nervous System Multiparameter Optimization (CNS MPO) Approach To Enable Alignment
    of Druglike Properties."* **ACS Chem Neurosci 2010;1(6):435–449**，
    doi:10.1021/cn100008c（PMID 22778837）—— **拐点数值在这一篇**
  - Wager TT, Hou X, Verhoest PR, Villalobos A. *"Central Nervous System Multiparameter
    Optimization Desirability: Application in Drug Discovery."*
    **ACS Chem Neurosci 2016;7(6):767–775**，doi:10.1021/acschemneuro.6b00029（PMID 26991242）
  - 同期配套文（性质空间的数据来源，119 个上市 CNS 药 + 108 个辉瑞候选）：
    Wager TT, Chandrasekaran RY, Hou X, Troutman MD, Verhoest PR, Villalobos A, Will Y.
    *"Defining Desirable Central Nervous System Drug Space through the Alignment of Molecular
    Properties, in Vitro ADME, and Safety Attributes."*
    **ACS Chem Neurosci 2010;1(6):420–434**，doi:10.1021/cn100007x

  六项 = MW / cLogP / **cLogD7.4** / TPSA / HBD / **最碱性 pKa**，
  各自过一条 desirability 曲线映射到 0–1 再相加，得 0–6；TPSA 是唯一的「驼峰」函数。
  **RDKit 只给得出四项**，logD7.4 与 pKa 得靠 ADMETlab（`logD` / `pka_basic` 列）。

  ✅ **拐点已从原文 PDF 取得**（`Step3_GKA_Brain_Penetration_Prediction/cn100008c.pdf`，
  Table 1 + Figure 4；变换公式见 Methods eq 1/2：拐点之间**线性**插值、六项**等权**求和）。
  实现在 `Step3_04_Result_Integration/Step3_04_CNS_MPO.py`：

  | 项 | 形状 | T0 = 1.0 | T0 = 0.0 |
  |---|---|---|---|
  | cLogP | 单调下降 | ≤ 3 | > 5 |
  | cLogD7.4 | 单调下降 | ≤ 2 | > 4 |
  | MW | 单调下降 | ≤ 360 | > 500 |
  | TPSA | **驼峰** | 40 < TPSA ≤ 90 | ≤ 20 或 > 120 |
  | HBD | 单调下降 | ≤ 0.5 | > 3.5 |
  | 最碱性 pKa | 单调下降 | ≤ 8 | > 10 |

  **自检写死在脚本里**：原文 **Table 4（p.446）的算例输入未被四舍五入**，
  六项 T0 与总分 4.0 能对到小数点后两位；再加 Table 3 三个候选（容差 0.03）。
  四组算例任一项对不上，脚本直接退出、拒绝出分。
  ⚠ 三处**必须显式选择**的口径（原文用 BioByte ClogP / ACD logD 与 pKa，本项目都没有）：
  ClogP 取 `admetlab_logp`（与同源的 logD/pKa 一致）、TPSA 取 **N/O 口径 `tpsa`**
  （原文 ref 9 是 Ertl 2000，**与 BOILED-Egg 用含 S/P 的相反**）、MW 取平均分子量。
  换口径的三个敏感性版本一并算出（`cnsmpo_score_tpsa_sandp` / `_logp_rdkit` / `_logp_swiss`），
  逐分子给 `cnsmpo_score_variant_spread`：GKA 候选中位 **0.49**，对照只有 0.16——
  **候选这一侧对口径明显更敏感**，Step3_05 用分数排序前要先看这一列。
  ⚠⚠ 原文 p.446 明写 *"the algorithm is **not intended to be used purely as a predictor
  of CNS penetration**"*——它是**成药性对齐**工具，不是入脑预测器。
  实测佐证：本项目 487 个入脑对照上，CNS MPO 中位**阳性 4.8 / 阴性 4.4，几乎分不开**
  （而两个工具的 BBB 项分得很开）。**拿它当入脑判据前必须先面对这件事。**
- **⚠ SwissADME 的 `BBB permeant` 是 (WLOGP, TPSA) 的确定性几何规则，不是独立模型。**
  实测 141 个 "No" 中 **0 个**落在 58 个 "Yes" 点的凸包内部，完全可分。它就是 BOILED-Egg 的卵黄椭圆。
  对比之下 `Pgp substrate` 有 81 个 "No" 落在 "Yes" 凸包内——**它才携带独立信息**。
  **下游不能单用 `BBB permeant`，必须 `BBB+ 且非 Pgp 底物` 组合**，否则看不见外排。
- **⚠ SwissADME 现状（2026-08 实测）**：`iLOGP` 与 **5 个 CYP 抑制预测整列 `n/d`**（199/199 全空），
  应是其首页公告的「ChemAxon 停止对学术网站支持、5 月 15 日改版」的后果。
  副作用：`Consensus Log P` 变成**四者均值**（XLOGP3/WLOGP/MLOGP/Silicos-IT），不再是五者。
  另：**ChemAxon 的学术免费路径已关**，别再按老经验指望它。
- **两个工具的输入格式差异是链路风险的根源**：

  | | SwissADME | ADMETlab 3.0 |
  |---|---|---|
  | 格式 | `SMILES<空格>名称`，一行一个 | **单列 CSV（表头 `SMILES`）或裸 TXT** |
  | 能否带 ID | ✅ 能 | ❌ **不能** |
  | 批量上限 | FAQ 明写 **≤200**，且要求**串行** | 页面写 MAX 1000，**但 646/批实测 504 超时**，98/批可行 |
  | 回填 | 按名称 | **只能按行序** |

  **名称一律用 `mol_id`，不能用化合物名**——1,274 个里 93 个名字含空格，会被 SwissADME 截断。
  ADMETlab 无 ID 字段，行序一旦错位后面全错且不报错；**好在它返回 `raw_smiles` 原样**，
  可按结构做独立于行序的校验（实测 99/99 逐行一致、InChIKey 全同）。
- **ADMETlab 的 API 存在但主计算端点是坏的**（2026-08-02 实测）：
  `/api/single/admet` → `KeyError: "['BSEP'] not in index"`（服务端 bug，必然 500）；
  `/api/admetCSV` → `FileNotFoundError`；`/api/uploadfile` 返回 `null`（试过 6 种字段名）。
  **只有 `/api/washmol` 可用**（SMILES 数组 → 标准化去盐）。走网页上传。
- **ADMETlab 返回 123 列，CNS MPO 缺的两个值都在里面**（batch1 实测，99/99 无缺失）：

  | 列 | 用途 |
  |---|---|
  | **`logD`** / **`pka_basic`** / `pka_acidic` / `logP` | **补齐 CNS MPO 的第 3、6 项** |
  | `BBB` | ⚠ **0–1 概率值，不是 Yes/No 标签**，与 SwissADME 的二分类不能直接对齐 |
  | `pgp_sub` / `pgp_inh` / `BCRP` / `MRP1` / `OATP1B1` / `OATP1B3` / `BSEP` | 外排/摄取转运体，比 SwissADME 只有一个 `pgp_sub` 丰富得多 |
  | `caco2` / `MDCK` / `PAMPA` / `PPB` / `Fu` / `logVDss` | 通透性与结合 |

  **有了它就不必装 OPERA**（NIEHS 那个免费本地 QSAR，也给 logD7.4 + pKa，
  且带适用域评估）——但 ADMETlab 若长期不可用，OPERA 是首选替补，
  网页替补是 pkCSM（给 logBB / logPS / P-gp，但**不给 logD 和 pKa**）。
- **`raw_smiles` 列让 ADMETlab 的行序风险降了一档**：它把输入原样返回，
  所以除了行序还能按结构独立校验。**16 批全量实测：行序 16/16 逐字一致、
  工具标准化后 InChIKey 1,571/1,571 全同、锚点漂移只有 3 处且差在小数点后第 8 位**
  （`0.1683771461` vs `0.1683771610`，浮点噪声不是模型漂移）。
  **ADMETlab 的可复现性明显好于 SwissADME**，后者出过返回错误分子的事故（见下）。
- **⚠ 两个工具在 P-gp 上给出几乎相反的图景**，而这恰恰是判断
  「像 CNS 药但被外排挡住」的关键项：

  | | SwissADME `Pgp substrate` | ADMETlab `pgp_sub` 中位 |
  |---|---|---|
  | GKA 候选 | **51.7% 判为底物** | **0.001** |
  | B3DB 对照 | 47.1% | 0.192 |
  | Fridén | 50.0% | 0.710 |

  不同模型、不同训练集、不同输出形式（二分类 vs 概率），分歧本身不算错误。
  但 **Step3_05 讨论「多工具结论不一致怎么办」时，P-gp 会是最突出的分歧点**，
  不能简单取交集或多数——BBB 项两个工具高度一致，P-gp 项则近乎相反。
- **每批都要放同一组锚点分子重复提交。** 网页模型会静默更新，多批不是同时跑的；
  同一锚点在不同批给出不同结果就是漂移的直接证据。且**各批的集合构成必须均匀混合**——
  若把候选与对照分装不同批，漂移与「候选 vs 对照」的真实差异就分不开了。
  两个工具用**同一组 20 个锚点**，结果才好互相印证。
- **⚠⚠ 网页工具会返回你从未提交过的分子——实测抓到过一次。**
  SwissADME batch6 的 `B3D_0441`（trifluopromazine，C18H19F3N2S，MW 352.42）
  返回的是 `O=C1[C@@]2(C)CC[C@@H](C1(C)C)C2` = **樟脑**（C10H16O，MW 152.23），
  连 `BBB permeant` 都从 No 翻成 Yes。**樟脑根本不在本项目的 1,274 个分子里**，
  推测是浏览器会话残留。名称是对的，只有结构是错的——**光核对名称发现不了。**

  1,391 条里就这 1 条（其余 1,390 条结构逐一相符），是孤立单行不是整批错位。
  **两条独立路径同时抓到它**：锚点的 7 次重复（另外 6 次全对）＋ 提交/返回的全量结构核对。
  **两道校验都不能省**——锚点只覆盖 20 个分子，非锚点分子出事只有结构核对能发现。

  ⇒ **凡是网页工具的结果，回填前必须逐行比对「提交的结构」与「返回的结构」的 InChIKey，
  不符的整行剔除并记录，绝不能只按名称或行序对。**

  Step3_04 把这道校验写进了脚本并独立重现了这一条（`Step3_04_Verification_Failures.csv`）。
  **剔除它有实质影响**：`B3D_0441` 是 `control_positive`，错的那次给 `BBB permeant = Yes`，
  其余 6 次全是 `No`——不剔除就会被众数合并带偏。
  ADMETlab 侧 1,571 行全部通过（行序 16/16 批对上、`raw_smiles` 逐字一致）。
- **⚠ 两个工具的 `MW` 口径不同，别混着比阈值。**（Step3_04 实测）
  **SwissADME 的 `MW` 是平均分子量**（与 RDKit `MolWt` 一致，最大差 0.02 Da），
  **ADMETlab 的 `MW` 是单同位素质量**（与 `ExactMolWt` 一致，最大差 0.005 Da）。
  1,274 个里 **73 个两者差 >1 Da，最大 3.07 Da，全部落在含 Cl（45）/ Br（32）的分子上**。
  差值对 CNS MPO 的 MW 项影响极小，但**拿 `admetlab_mw` 去比 Lipinski 500 这类硬边界，
  含卤分子会系统性偏低**。整合表里四个口径都在：`mw` / `mw_exact` / `swissadme_mw` / `admetlab_mw`。
  ⚠ 这也是自检的一个假警报来源——初版拿 `admetlab_mw` 比 `mw`，报了「73 行回填可能错位」，
  实际是口径差异。**交叉校验要比同一个量。**
- **GKA 候选与入脑对照的化学空间几乎不重叠，且这不是去冗余造成的**（已双向验证）：

  | | n | MW 中位 | TPSA 中位（不含 S/P） |
  |---|---:|---:|---:|
  | BBB+ 对照 | 351 | 319 | 44 |
  | **BBB− 对照** | 94 | **366** | **89** |
  | GKA 候选 | 787 | **463** | **105** |

  **GKA 候选的 MW 与 TPSA 超过了阴性对照**，不是落在两类之间。
  验证过两头：选取只把 B3DB 往 GKA 方向拉近（合格池 292.5 → 选中 319.4），
  而未经任何筛选的 B3DB 全表（1,058）是 MW 309、Fridén 42 个是 MW 292；
  GKA 侧换五种口径（逐分子 / 每骨架 1 个 / 簇内最小 MW / 只看 226 个单例骨架 / 只看最大 10 簇）
  MW 全在 460–475、TPSA 全在 104–115。**差异是固有的。**
  ⚠ 由此引出 Step3_05 必须面对的**外推问题**：
  对照集再怎么挑也够不到 GKA 那个区间（B3DB 全池 487 个 BBB+ 里四项同时落进 GKA
  5–95 分位区间的只有 37 个），**流程在对照集上分得开 ≠ 在 GKA 上同样可靠**。

## 方法论：两次踩同一个坑

- **⚠⚠ 阈值必须用阳性对照标定，且自检要写死在脚本里。**
  这条在 Step1_05 踩过一次（初版 7.0/6.5 只过 4/6，改 6.5/6.0 后 6/6），
  **Step2_02 又踩了第二次**：看着分布拍了 5 个阈值、一个对照没验，
  结果 **12 个已知 GKA 只有 3 个通过**。被误杀的两条规则及其错因：

  | 规则 | 误杀 | 为什么错 |
  |---|---:|---|
  | `specificity >= 0.01` | 2 | **越重要的药被后续专利引用越多**，全库出现数越大、特异性反而越低（AZD-1656 0.0063、Piraglitin 0.0072）。对成名已久的药有系统性偏见 |
  | `n_in_sel >= 2` | 5 | 理由本身就错：**一个药通常只在它自己那一族专利的权要里被具体画出来**，别处提它是当现有技术写在说明书。「只出现在 1 篇权要里」恰恰是原研化合物的特征 |

  两条已降级为标注列。**真正的修补不是调阈值，是把自检做成脚本的一部分**——
  Step1_05 的自检写死在脚本里所以守住了，Step2_02 初版没有所以没守住。
  现在 `selfcheck_controls()` 误杀会直接报警。
- **召回优先的步骤，宁可宽也不能把已知的药筛掉。** 收窄留给后面有方向证据的步骤。
- **⚠ 不对外部数据源下「对错」判决，只陈述差异与成因，判断权留给人。**
  别人的库可能采用不同标准、不同测量情景、不同适用场景。
  B3DB 用 `logBB > −1` 判 BBB±，这个定义本身自洽；它与 Kp,uu 在 P-gp 底物上分歧，
  是**两个量测的东西不同**，不是谁标错了。
  写法上：「判据 A 判为 X，判据 B 判为 Y，差异集中在 Z 类化合物，机制是……」，
  不写「❌ / 标反了 / 错在……」。选用某个来源要表述成**带理由的选择**，
  而不是对落选者的否决。差异如实落进产物（列名用 `differ` 而非 `conflict`），让人能自己复核。
- **别拿硬指标凑软目标。** 固定数量（"选够 100 个"）会逼着在「凑数」与「保多样性」之间二选一：
  Step3_00 阴性侧池子只有 105 个分子 / 77 个骨架，硬凑 100 就等于放弃骨架去冗余。
  正确做法是**放弃固定 N，合格多少要多少**，把约束写成规则而不是数字。
- **骨架去冗余要软，不要硬。** 「一骨架一分子」会把证据充分的同系物一并砍掉；
  同一 Murcko 骨架保留 2–3 个代表，簇内按证据质量排。
  ⚠ **无环分子的 Murcko 骨架是空串**——全当成一簇会把彼此毫不相干的分子挤在一起
  （B3DB 阳性池 487 个里有 70 个无环）。是否合并要显式选择，不能默认。
- **⚠ 去重的键必须在标准化之后取。** Step3_00 用标准化前的 InChIKey 去重，
  漏掉了 6 对同一化合物的不同写法（4 对立体标注详略不同、2 对结构本身差一个 CH₂）。
  跨库比对至少要看**标准化后的完整 InChIKey**，必要时再看前 14 位（骨架层）与名称。

## 各步骤约定

- 每步一个目录：`Readme.md`（任务定义）+ 脚本 + `.csv`（主产物）+ `.md`（报告）。
- 每步显式读取上一步的 CSV，保持链路可追溯；也留 `--tid` 之类的参数可覆盖。
- 报告开头必须记录出处：ChEMBL 版本、数据库路径、运行时间、输入来源。
- 多值字段聚合成 JSON 列表写进单元格，**保持一行一个实体**。
- 空值也是事实，如实记录（如 `variant_id` 0/228），不要省略字段。
- **过滤一律「加列不删行」**：算指标 → 打 `keep` / `exclude_reason`，被排除的行留在 CSV 里，
  可复核可反悔、改阈值不用重跑上游。
- **凡是有阳性对照的步骤，自检必须写进脚本**，误杀要报警（见「方法论」一节）。
- **人工干预要留痕。** 手工删掉的分子记进 `Step3_02_Manual_Exclusions.csv` 那种表
  （mol_id / 批次 / 时间 / 操作人 / 原因 / 当时的数值依据），下游脚本读它，
  把「人工排除」与「真正的缺失」区分开，否则对账会误报。
  已记：`B3D_0012`（MW 1802.7、127 重原子，全表唯一 MW>1000，超出 SwissADME 适用范围）。
- **外部工具的提交与回填要拆成两半**：导出提交文件（带逐行清单）+ `parse_<tool>()` 导入并校验。
  网页工具常改列名、丢行、改大小写，**对不上要报错而不是静默丢行**。
  这样上游可以先本地跑通，外部结果到位后直接并入，不必重跑。

## 规则 + LLM 分类的分工

`Step1_03` 拆成两个脚本，规则在前，LLM 只处理规则判不明确的记录：

- 规则给出 high 置信度的结论，**LLM 不得覆盖**。
- LLM 必须输出**逐字引自原描述的证据句**；脚本校验该句是否真实存在于原文，
  核验不通过则不采纳标签、直接转人工。这是防编造的硬约束。
- LLM 返回"无法判断"或低置信度时保持 `review_required`，**不允许强行分类**。
- 模型 `deepseek-v4-pro`，key 读环境变量 `DEEPSEEK_API_KEY`。
  它是推理模型，`max_tokens` 给小了正文会全空（推理 token 吃光预算），默认 3000。
- 响应按输入哈希缓存到 `.llm_cache.json`（已 gitignore），重跑不重复计费。

## 环境

micromamba 环境 **`GKA_in_Brain`**。micromamba 不在 PATH 上（`.bashrc` 里只有 alias，
非交互 shell 取不到），脚本和命令里一律写全路径：

```
/home/sgao30/micromamba/bin/micromamba run -n GKA_in_Brain python xxx.py
```

| | |
|---|---|
| python | 3.11.15 |
| rdkit | 2026.03.4 |
| numpy | 2.4.6 |
| pandas | 3.0.5 |
| pyarrow | 25.0.0（读 parquet） |
| duckdb | 1.5.5（直查 SureChEMBL 大表） |
| requests | 2.34.2（Step1_03 调 LLM 用） |
| sqlite3 模块 | 3.53.4 |

没有 `sqlite3` 命令行，用 python 的 `sqlite3` 模块。

- **不要用系统 python（3.8.10）。** 那里的 rdkit 缺 numpy，
  `Descriptors` / `Crippen` / `MurckoScaffold` 全部 import 失败，
  只有 `rdMolDescriptors` 这类纯 C++ 接口能用。
- Step1_01–04 的脚本在 3.11 下重跑，产物与 3.8 逐字节一致（只有报告时间戳变）。
- **f-string 表达式内仍不能含反斜杠**——这条限制是 3.12（PEP 701）才放开的，
  3.11 下照旧报 `SyntaxError`。用 `chr(92)` 绕开。
