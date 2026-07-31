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
- **PPI 靶点不能与主靶点无脑合并**：`CHEMBL3885579` 测的是 GCK–GKRP 相互作用
  （GKRP 解离剂机制），与直接激活 GCK 是两类机制。
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

Python 3.8.10，标准库 + `requests`。没有 `sqlite3` 命令行，用 python 的 `sqlite3` 模块。
注意 3.8 下 f-string 表达式内不能含反斜杠。
