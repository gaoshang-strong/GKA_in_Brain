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
```

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

**下游要挑分子，直接读 `Step1_06_GKA_Physicochemical_Properties.csv`**——
它带着 Step1_05 的 `priority` / `potency_band` / `pactivity_median` / `direction` /
骨架列，是链路末端唯一需要读的表。

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
- **⚠ 同一个药在 ChEMBL 里有多个 chembl_id，`molecule_hierarchy` 归并不全。**
  MK-0941 = `CHEMBL3580737`（游离碱，有活性）+ `CHEMBL4297302`（药物条目，0 活性）；
  Globalagliatin = `CHEMBL4297399`（LY-2608204，研发代号）+ `CHEMBL5095182`（盐酸盐）。
  这几对之间**没有 parent 关系**。**去重要按结构（InChIKey）或人工归并，不能按 ID 计数。**
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

## 各步骤约定

- 每步一个目录：`Readme.md`（任务定义）+ 脚本 + `.csv`（主产物）+ `.md`（报告）。
- 每步显式读取上一步的 CSV，保持链路可追溯；也留 `--tid` 之类的参数可覆盖。
- 报告开头必须记录出处：ChEMBL 版本、数据库路径、运行时间、输入来源。
- 多值字段聚合成 JSON 列表写进单元格，**保持一行一个实体**。
- 空值也是事实，如实记录（如 `variant_id` 0/228），不要省略字段。

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
| requests | 2.34.2（Step1_03 调 LLM 用） |
| sqlite3 模块 | 3.53.4 |

没有 `sqlite3` 命令行，用 python 的 `sqlite3` 模块。

- **不要用系统 python（3.8.10）。** 那里的 rdkit 缺 numpy，
  `Descriptors` / `Crippen` / `MurckoScaffold` 全部 import 失败，
  只有 `rdMolDescriptors` 这类纯 C++ 接口能用。
- Step1_01–04 的脚本在 3.11 下重跑，产物与 3.8 逐字节一致（只有报告时间戳变）。
- **f-string 表达式内仍不能含反斜杠**——这条限制是 3.12（PEP 701）才放开的，
  3.11 下照旧报 `SyntaxError`。用 `chr(92)` 绕开。
