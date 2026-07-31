# Step1 候选分子与理化性质整合表

- 运行时间：2026-07-31 18:16:42
- 输入：`Step1_06_GKA_Physicochemical_Property_Extraction/Step1_06_GKA_Physicochemical_Properties.csv`、`Step1_07_GKA_from_Drug_Annotation/Step1_07_GKA_from_Drug_Annotation.csv`
- 合计 **787** 个 ChEMBL ID

> 本表是 **Step1 全链路的最终产物**，Step2 直接读它。两条来源不同的路径合并在一起，`source` 列标明每一行从哪来。

## 一、两条来源，口径不同

| source | 分子数 | 有效力值 | 方向来自 | 能否排序 |
| --- | ---: | :---: | --- | :---: |
| `activity` | 778 | ✅ | 从 assay 读数推（Step1_05） | ✅ |
| `activity+drug_annotation` | 4 | ✅ | **两者都有**，人工审编可复核推断结果 | ✅ |
| `drug_annotation` | 5 | ❌ | **ChEMBL 人工审编**（`action_type`） | ❌ |

**`drug_annotation` 那批没有 `pactivity_median` / `priority` / `potency_band`——这是事实不是缺失。** 它们在 ChEMBL 里一条活性记录都没有（或没有打在 GCK 上），无法参与效力分档与排序。但它们的方向是**人工审编**的，比推断的可靠。

只靠药物注释路径进来的（Step1_07 补回的）：

| 分子 | 名称 | phase | 命中路径 | 审编方向 |
| --- | --- | ---: | --- | --- |
| `CHEMBL4297508` | DORZAGLIATIN | 3 | drug_mechanism+usan_stem+synonym | ACTIVATOR |
| `CHEMBL4297302` | MK-0941 | 2 | drug_mechanism | ACTIVATOR |
| `CHEMBL4297399` | LY-2608204 | 2 | drug_mechanism+usan_stem+synonym | ACTIVATOR |
| `CHEMBL5095262` | CADISEGLIATIN | 2 | usan_stem+synonym | — |
| `CHEMBL5095182` | GLOBALAGLIATIN HYDROCHLORIDE | 1 | usan_stem+synonym | — |

## 二、⚠ 去重键：`parent_chembl_id`，不是 InChIKey

**787 个 chembl_id → 785 个去重组**，其中 2 组含多个 ID。

| 去重组 | 成员 | 说明 |
| --- | --- | --- |
| `CHEMBL3580737` | `CHEMBL3580737`、`CHEMBL4297302`（盐） | MK-0941 (free base)、MK-0941 (mesylate) |
| `CHEMBL4297399` | `CHEMBL4297399`、`CHEMBL5095182`（盐） | LY-2608204 / Globalagliatin、Globalagliatin HCl |

**同一个药以「游离碱 + 盐」两个 chembl_id 存在时，两者 InChIKey 完全不同**——

```
MK-0941        CHEMBL3580737 游离碱   KJSGTWFWVTYPFZ-AWEZNQCLSA-N
               CHEMBL4297302 甲磺酸盐 PIDNRTWDGDJKSQ-UQKRIMTDSA-N  ← SMILES 带 .CS(=O)(=O)O
Globalagliatin CHEMBL4297399 游离碱   QIIVJLHCZUTGSD-CUBQBAPOSA-N
               CHEMBL5095182 盐酸盐   FRUQQNDJVRDIRH-JOFLZTHPSA-N  ← SMILES 带 Cl.
```

按 InChIKey 去重会把同一个药算成两个。`molecule_hierarchy.parent_molregno` 能正确归并，所以**药物层去重用 `parent_chembl_id`**；跨库对齐结构（如与 SureChEMBL）才用 InChIKey，且要先归到母体。**两个场景用不同的键。**

产物保留一行一个 chembl_id，`is_dedup_representative` 标出每组代表（优先非盐、有活性证据），**不擅自删行**。

## 三、理化性质覆盖

| 字段 | 有值 | 缺失 |
| --- | ---: | ---: |
| `canonical_smiles` | 787 | 0 |
| `standard_inchi_key` | 787 | 0 |
| `mw_freebase` | 787 | 0 |
| `alogp` | 787 | 0 |
| `psa` | 787 | 0 |
| `hbd` | 787 | 0 |
| `hba` | 787 | 0 |
| `rtb` | 787 | 0 |
| `qed_weighted` | 787 | 0 |
| `murcko_scaffold` | 787 | 0 |

- **MW**：中位 462.55，落在常用 CNS 窗口（≤ 450）的 327 / 787
- **ALogP**：中位 4.15，落在常用 CNS 窗口（1–5）的 617 / 787
- **PSA**：中位 105.25，落在常用 CNS 窗口（≤ 90）的 181 / 787

> 性质窗口只作描述，**本表不做任何脑暴露筛选**——留给 Step2。另注意 ChEMBL 37 没有 `cx_logd` / `molecular_species`，`alogp` 不是 logD。

## 四、阳性对照在表内的位置

**12 / 12** 个对照在表内。

| 分子 | 名称 | phase | source | 效力档 | pAct 中位 |
| --- | --- | ---: | --- | --- | ---: |
| `CHEMBL4297508` | Dorzagliatin | 3 | drug_annotation | — | — |
| `CHEMBL3219124` | AZD-1656 | 2 | activity+drug_annotation | ≥7.5（EC50 ≤ 32 nM） | 7.55 |
| `CHEMBL2165620` | PF-04991532 | 2 | activity+drug_annotation | 7.0–7.5（32–100 nM） | 7.05 |
| `CHEMBL2165615` | Neriglitin | 2 | activity+drug_annotation | 6.5–7.0（0.1–0.32 µM） | 6.725 |
| `CHEMBL1783734` | Piraglitin | 2 | activity+drug_annotation | 6.0–6.5（0.32–1 µM） | 6.145 |
| `CHEMBL3580737` | MK-0941 (free base) | 2 | activity | ≥7.5（EC50 ≤ 32 nM） | 8.02 |
| `CHEMBL4297302` | MK-0941 (mesylate) | 2 | drug_annotation | — | — |
| `CHEMBL4297399` | LY-2608204 / Globalagliatin | 2 | drug_annotation | — | — |
| `CHEMBL5095262` | Cadisegliatin | 2 | drug_annotation | — | — |
| `CHEMBL5095182` | Globalagliatin HCl | 1 | drug_annotation | — | — |
| `CHEMBL1096435` | Ro-28-1675 (参比化合物) | — | activity | 6.0–6.5（0.32–1 µM） | 6.39 |
| `CHEMBL5072532` | BMS-820132 | — | activity | 7.0–7.5（32–100 nM） | 7.28 |

## 五、下游怎么用

| 想做的事 | 该怎么筛 |
| --- | --- |
| 按效力挑候选 | `source` 含 `activity` 且 `priority` / `potency_band` 非空 |
| 要方向最可靠的 | `curated_direction = 'ACTIVATOR'`（ChEMBL 人工审编） |
| 结构去冗余 | `is_scaffold_representative = TRUE` |
| 药物层计数 | 按 `dedup_group` 或 `is_dedup_representative = TRUE` |
| 与 SureChEMBL 对齐 | `standard_inchi_key`，**先归母体再取 key** |

