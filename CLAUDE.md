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
```

已锚定的 GCK 靶点（ChEMBL 37）：

| target | tid | type | assays | activities |
|---|---|---|---|---|
| `CHEMBL3820` Hexokinase-4 | 20095 | SINGLE PROTEIN | 227 | 3,222 |
| `CHEMBL3885579` Glucokinase/GKRP | 117123 | PROTEIN-PROTEIN INTERACTION | 1 | 40 |

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
  还有极性问题：`S0.5`/`S50`/`Km` 是**值越小激活越强**，与 `FC` 方向相反；
  它们测的是酶对葡萄糖的半饱和常数（GKA 把它从野生型 ~7 mM 拉到 0.5–2 mM），
  是机制读数，**不是化合物效力**，当效力用会得出「EC50 = 0.6 mM」这种量级错误。
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
