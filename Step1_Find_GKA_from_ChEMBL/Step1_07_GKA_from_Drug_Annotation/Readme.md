补 Step1 的召回缺口：把**因为没有活性数据而被整条链路漏掉**的 GKA 找回来。

## 为什么会漏

Step1_01–06 的链路是：

```
靶点 → assays → 激活 assay → activities → 分子
```

**每一步都挂在 `activities` 上。零活性记录的分子在结构上不可见。**

ChEMBL 的分子有两条独立入口：

1. **从文献活性数据抽取** —— Step1 走的这条
2. **从药名/临床登记册收录**（USAN / INN / ATC / Clinical Candidates）—— 完全没走

多格列艾汀（Dorzagliatin，`CHEMBL4297508`，III 期，中国已上市）就是典型：
`compound_records` 里 4 条来源**全是 `doc_type = DATASET`（USAN / INN / ATC / CANDIDATES），
没有一篇论文**，`activities = 0`。它在 ChEMBL 里，但我们看不见它。

## 从哪几张表找

三条独立路径，各自记录证据，**产物里逐分子写明是被哪条路径命中的**。

### 路径 A：`drug_mechanism` —— 人工审编的方向标注 ⭐

```
drug_mechanism.molregno → molecule_dictionary
drug_mechanism.tid      → target_dictionary   （筛 CHEMBL3820 Hexokinase-4）
```

取 `mechanism_of_action`、`action_type`、`direct_interaction`、`molecular_mechanism`。

**这张表直接给出方向**：挂在 `CHEMBL3820` 上的 6 条，`action_type` 全部是 **`ACTIVATOR`**，
`mechanism_of_action` 全部是 **"Hexokinase type IV activator"**。

这正是 Step1_03/05 用规则 + LLM 从 assay 描述里辛苦推的那个结论，
ChEMBL 已经人工标好了。**比我们推的可靠，且免费。**

**安全网**：同时不限 `tid` 扫一遍
`mechanism_of_action LIKE '%glucokinase%' OR '%hexokinase type IV%'`，
防止有 GKA 被挂到别的靶点上。实测结果与按 `tid` 筛完全一致（都是那 6 条），
但这条检查要留在脚本里，换版本时能自动发现变化。

### 路径 B：`usan_stems` —— 按药名词干识别 ⭐

```
usan_stems.stem = '-gliatin'  →  annotation = 'glucokinase activator'
molecule_dictionary.usan_stem = '-gliatin'
```

WHO/USAN 的命名规则里，**后缀 `-gliatin` 的定义就是「glucokinase activator」**
（`usan_stems` 表原文，`molecule_dictionary.usan_stem_definition` 同）。
**按名字就能认出 GKA，零成本、零歧义。**

### 路径 C：`molecule_synonyms` —— 兜住 `usan_stem` 字段没填的

```
molecule_synonyms.synonyms LIKE '%gliatin%'
```

有些分子 `usan_stem` 字段是空的，但同义词里有 `-gliatin` 名。
例：`CHEMBL4297399` 主名是研发代号 `LY-2608204`，同义词才是 `Globalagliatin`。

## ⚠ 同一个药有多个 ChEMBL ID

必须用 `molecule_hierarchy` 归并，否则会重复计数：

| 药 | 条目 1 | 条目 2 |
|---|---|---|
| MK-0941 | `CHEMBL3580737` MK-0941 FREE BASE（有活性，**在 782 里**） | `CHEMBL4297302` MK-0941（药物条目，0 活性） |
| Globalagliatin | `CHEMBL4297399` LY-2608204（研发代号，6 条活性） | `CHEMBL5095182` GLOBALAGLIATIN HYDROCHLORIDE（盐酸盐，0 活性） |

产物保留全部条目（一行一个 ChEMBL ID），另给 `parent_chembl_id` 与 `dedup_group` 列供归并。

## 口径：不合并进 Step1_05/06 的主候选表

本步骤的分子**来源与 782 个候选不同**——它们没有活性数值，无法参与效力分档与排序。
产物**单独成表**，带 `source = drug_annotation` 标记。下游要合并时必须知道：

- 782 个候选有 `pactivity`，可排序；这批**没有**，只有方向标注
- 这批的方向是**人工审编的**（`action_type = ACTIVATOR`），比 Step1_05 从读数推的更可靠

## 本步骤不做什么

| 不做 | 为什么 |
|---|---|
| 不返工 Step1_03/04/05 | 分类与排序方法是对的，缺的是入口，不是逻辑 |
| 不并入「GCK–GKRP 相互作用」assay 的 115 个分子 | 另一类机制，与单列 `CHEMBL3885579` 的处理一致 |
| 不并入「GCK 抑制」assay 的 50 个分子 | 方向相反 |
| 不并入「无法判断」/「结合」assay 的 18 个分子 | 方向不明，收益抵不上把口径搞浑 |

## 输出

- `Step1_07_GKA_from_Drug_Annotation.csv` —— 一行一个 ChEMBL ID，含三条路径的命中证据、
  结构、理化性质（字段与 Step1_06 对齐）、活性计数、是否已在 782 里
- `Step1_07_GKA_from_Drug_Annotation.md` —— 报告
- 顺带更新阳性对照集：**8 → 11**（Step2 自检用）

代码：`Step1_07_GKA_from_Drug_Annotation.py`
