# Step1_07 从药物注释补回被漏掉的 GKA

- ChEMBL 版本：**CHEMBL_37（2026-05-01）**
- 数据库文件：`/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db`
- 运行时间：2026-07-31 18:11:31
- 命中分子：**9** 个 ChEMBL ID，其中 **5** 个不在 Step1_06 的 782 个候选里

> Step1_01–06 的链路每一步都挂在 `activities` 上，**零活性记录的分子在结构上不可见**。本步骤从药物注释侧独立锚定，把它们找回来。

## 一、从哪几张表找的

| 路径 | 表 | 判据 | 命中 |
| --- | --- | --- | ---: |
| **A** | `drug_mechanism` × `target_dictionary` | `tid` → `CHEMBL3820` 且 `action_type` 属激活类 | 6 |
| **B** | `usan_stems` × `molecule_dictionary` | `usan_stem = '-gliatin'` | 6 |
| **C** | `molecule_synonyms` | `synonyms LIKE '%gliatin%'` | 6 |

路径 B 的依据是 `usan_stems` 表原文：**`-gliatin` → 「glucokinase activator」**。WHO/USAN 的命名规则里这个后缀就是给 GKA 用的，**按药名就能认出来**。

路径 A 的 `action_type` 实测分布：`ACTIVATOR` × 6。

**路径 A 给的是人工审编的方向标注**——`mechanism_of_action` 写着 「Hexokinase type IV activator」。这正是 Step1_03/05 用规则 + LLM 从 assay 描述里推的结论，ChEMBL 已经标好了，比推的可靠。

**安全网**：另外不限 `tid` 扫了一遍 `mechanism_of_action` 文本（含 glucokinase / hexokinase type IV），命中 0 条，与按 `tid` 筛的结果完全一致，没有 GKA 被挂到别的靶点上。

## 二、命中清单

| 分子 | 名称 | phase | 命中路径 | GCK 活性数 | 总活性数 | 在 782 里 |
| --- | --- | ---: | --- | ---: | ---: | :---: |
| `CHEMBL4297508` | DORZAGLIATIN | 3 | drug_mechanism+usan_stem+synonym | 0 | 0 | **❌ 新增** |
| `CHEMBL1783734` | PIRAGLIATIN | 2 | usan_stem+synonym | 19 | 140 | ✅ |
| `CHEMBL2165615` | NERIGLIATIN | 2 | drug_mechanism+usan_stem+synonym | 13 | 72 | ✅ |
| `CHEMBL2165620` | PF-04991532 | 2 | drug_mechanism | 1 | 101 | ✅ |
| `CHEMBL3219124` | AZD-1656 | 2 | drug_mechanism | 6 | 23 | ✅ |
| `CHEMBL4297302` | MK-0941 | 2 | drug_mechanism | 0 | 0 | **❌ 新增** |
| `CHEMBL4297399` | LY-2608204 | 2 | drug_mechanism+usan_stem+synonym | 0 | 6 | **❌ 新增** |
| `CHEMBL5095262` | CADISEGLIATIN | 2 | usan_stem+synonym | 0 | 0 | **❌ 新增** |
| `CHEMBL5095182` | GLOBALAGLIATIN HYDROCHLORIDE | 1 | usan_stem+synonym | 0 | 0 | **❌ 新增** |

## 三、为什么以前会漏

**4 / 9** 个分子的 `activities` 为 0。看它们的 `compound_records` 来源就明白了：

| 分子 | 来源 | 有文献记录 |
| --- | --- | :---: |
| `CHEMBL4297302` MK-0941 | `CANDIDATES`(DATASET) | ❌ |
| `CHEMBL4297508` DORZAGLIATIN | `USAN`(DATASET)、`CANDIDATES`(DATASET)、`ATC`(DATASET)、`INN`(DATASET) | ❌ |
| `CHEMBL5095182` GLOBALAGLIATIN HYDROCHLORIDE | `CANDIDATES`(DATASET) | ❌ |
| `CHEMBL5095262` CADISEGLIATIN | `USAN`(DATASET)、`INN`(DATASET)、`CANDIDATES`(DATASET) | ❌ |
| `CHEMBL4297399` LY-2608204 | `CANDIDATES`(DATASET)、`SARS_COV_2`(DATASET)、`HDAC6`(DATASET) | ❌ |
| `CHEMBL3219124` AZD-1656 | `LITERATURE`(PUBLICATION)、`CANDIDATES`(DATASET) | ✅ |
| `CHEMBL2165615` NERIGLIATIN | `LITERATURE`(PUBLICATION)、`ASTRAZENECA`(DATASET)、`CANDIDATES`(DATASET)、`SARS_COV_2`(DATASE | ✅ |
| `CHEMBL2165620` PF-04991532 | `LITERATURE`(PUBLICATION)、`CANDIDATES`(DATASET) | ✅ |
| `CHEMBL1783734` PIRAGLIATIN | `LITERATURE`(PUBLICATION)、`USAN`(DATASET)、`PUBCHEM_BIOASSAY`(DATASET)、`ASTRAZENECA`(DATASE | ✅ |

`doc_type = DATASET` 的来源（USAN / INN / ATC / CANDIDATES）是**药名与临床登记册**，不带活性数值。分子经这条路进 ChEMBL，`activities` 就是 0，Step1 那条 `靶点 → assay → activity → 分子` 的链路自然看不见它。

## 四、⚠ 同一个药有多个 ChEMBL ID

`molecule_hierarchy` 的 `parent_molregno` 能归并一部分，但**不是全部**——同一个药以「游离碱」「盐」「药物条目」等多种身份注册时，彼此之间未必有 parent 关系。实测到的：

| 药 | 条目 | 说明 |
| --- | --- | --- |
| MK-0941 | `CHEMBL3580737`、`CHEMBL4297302` | 本表命中 2 个 |
| Globalagliatin / LY-2608204 | `CHEMBL4297399`、`CHEMBL5095182` | 本表命中 2 个 |

**做去重统计时必须按结构（InChIKey）或人工归并，不能只按 ChEMBL ID 计数。**

## 五、口径：不并入 782 个候选

本表的分子**没有活性数值**，无法参与 Step1_05 的效力分档与排序，因此单独成表。下游合并时必须知道两点：

1. 782 个候选有 `pactivity`、可排序；**这批没有**，只有方向标注
2. 这批的方向是 ChEMBL **人工审编**的（`action_type = ACTIVATOR`），比 Step1_05 从读数推出来的更可靠

## 六、更新后的阳性对照集

Step2 自检用。原 8 个 → 现 **11** 个。

| 分子 | 名称 | phase | InChIKey | 来源 |
| --- | --- | ---: | --- | --- |
| `CHEMBL4297508` | DORZAGLIATIN | 3 | `HMUMWSORCUWQJO-QAPCUYQASA-N` | Step1_07 |
| `CHEMBL1783734` | PIRAGLIATIN | 2 | `XEANIURBPHCHMG-SWLSCSKDSA-N` | Step1_07 |
| `CHEMBL2165615` | NERIGLIATIN | 2 | `MASKQITXHVYVFL-UHFFFAOYSA-N` | Step1_07 |
| `CHEMBL2165620` | PF-04991532 | 2 | `GKMLFBRLRVQVJO-ZDUSSCGKSA-N` | Step1_07 |
| `CHEMBL3219124` | AZD-1656 | 2 | `FJEJHJINOKKDCW-INIZCTEOSA-N` | Step1_07 |
| `CHEMBL4297302` | MK-0941 | 2 | `PIDNRTWDGDJKSQ-UQKRIMTDSA-N` | Step1_07 |
| `CHEMBL4297399` | LY-2608204 | 2 | `QIIVJLHCZUTGSD-CUBQBAPOSA-N` | Step1_07 |
| `CHEMBL5095262` | CADISEGLIATIN | 2 | `HPGJSAAUJGAMLV-QAQDUYKDSA-N` | Step1_07 |
| `CHEMBL5095182` | GLOBALAGLIATIN HYDROCHLORIDE | 1 | `FRUQQNDJVRDIRH-JOFLZTHPSA-N` | Step1_07 |

另加 Step1_06 里已有的 `CHEMBL1096435` Ro-28-1675（参比化合物，非临床药）与 `CHEMBL5072532` BMS-820132。

