# Step1_06 GKA 候选分子理化性质提取

- ChEMBL 版本：**CHEMBL_37（2026-05-01）**
- 数据库文件：`/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db`
- 运行时间：2026-07-31 14:42:09
- 输入：`Step1_05_Followup_Candidates.csv`（Step1_05 后续实验清单）
- 提取分子：**782** 个，全部命中

> **本步骤只取数，不做筛选、不做判断。** 性质窗口（能不能进脑）留给 Step2。

## 一、字段覆盖

| 字段 | 有值 | 缺失 | 说明 |
| --- | ---: | ---: | --- |
| `canonical_smiles` | 782 | 0 | 结构 |
| `standard_inchi_key` | 782 | 0 | 结构哈希，跨库对齐用 |
| `mw_freebase` | 782 | 0 | 游离碱分子量 |
| `full_mwt` | 782 | 0 | 含盐分子量 |
| `alogp` | 782 | 0 | 计算 logP，**不是 logD** |
| `psa` | 782 | 0 | 极性表面积（TPSA） |
| `hbd` | 782 | 0 | 氢键供体 |
| `hba` | 782 | 0 | 氢键受体 |
| `rtb` | 782 | 0 | 可旋转键 |
| `aromatic_rings` | 782 | 0 | 芳环数 |
| `heavy_atoms` | 782 | 0 | 重原子数 |
| `qed_weighted` | 782 | 0 | 类药性 QED (0–1) |
| `num_ro5_violations` | 782 | 0 | Lipinski 五规则违规数 |
| `ro3_pass` | 782 | 0 | 是否通过 Rule of 3（片段筛选用） |
| `np_likeness_score` | 782 | 0 | 天然产物相似度 |
| `full_molformula` | 782 | 0 | 分子式 |

## 二、ChEMBL 37 缺的性质（Step2 必须自己算）

这一版的 `compound_properties` 只有 15 列，**没有**下面这些。缺的恰好是判断脑暴露最关键的：

| 缺失字段 | 是什么 | 为什么这里要命 |
| --- | --- | --- |
| `cx_logd` | pH 7.4 下的分配系数 | GKA 里有相当一部分是羧酸，生理 pH 下带负电，logD 比 logP 低几个数量级。**只看 `alogp` 会高估膜通透性** |
| `molecular_species` | 酸 / 碱 / 中性 | 酸性化合物入脑普遍差，这是 CNS 项目的第一刀 |
| `cx_most_apka` / `cx_most_bpka` | 最强酸/碱解离常数 | 没有 pKa 就算不出 logD |

粗查 SMILES，**95 / 782** 个分子含羧酸基团（只是字符串匹配，Step2 应改用 SMARTS 正式判定）。本步骤把 `alogp` 如实取出并标注它**不是 logD**，不代替。

## 三、性质分布

**描述性统计，不构成筛选。** 阈值列只是让你先看到量级。

| 性质 | 中位 | 最小 | 最大 | 常用 CNS 窗口 | 落在窗口内 |
| --- | ---: | ---: | ---: | --- | ---: |
| MW (`mw_freebase`) | 462.55 | 255.30 | 658.80 | ≤ 450 | 327 |
| ALogP (`alogp`) | 4.14 | 1.46 | 7.67 | 1 – 5 | 615 |
| PSA (`psa`) | 105.25 | 37.91 | 171.65 | ≤ 90 | 179 |
| HBD (`hbd`) | 1.00 | 1.00 | 3.00 | ≤ 2 | 736 |
| HBA (`hba`) | 7.00 | 2.00 | 12.00 | ≤ 7 | 436 |
| RotB (`rtb`) | 7.00 | 2.00 | 15.00 | ≤ 8 | 534 |
| 芳环数 (`aromatic_rings`) | 3.00 | 1.00 | 5.00 | ≤ 3 | 509 |
| QED (`qed_weighted`) | 0.47 | 0.20 | 0.91 | ≥ 0.5 | 355 |

| Lipinski 违规数 | 分子数 |
| ---: | ---: |
| 0 | 467 |
| 1 | 193 |
| 2 | 119 |
| 3 | 3 |

六项一起看（MW / ALogP / PSA / HBD / HBA / RotB 全部落在窗口内）：**85 / 782**。

> 这批分子是冲着**肝和胰腺**做的，不是冲着脑做的——PSA 中位数已经超过常用 CNS 上限，是意料之中的结果。**这不是筛选结论**，Step2 会用正式的判据重做。

## 四、注解字段（大多为空，如实记录）

`-1` 是 ChEMBL 的「未标注」，不是 0。这批分子绝大多数是文献化合物，没进过开发流程，本来就不会有这些注解。

| 字段 | 取值分布 |
| --- | --- |
| `max_phase` | `(空)`×777、`2`×5 |
| `chirality` | `-1`×777、`1`×4、`2`×1 |
| `prodrug` | `-1`×777、`0`×5 |
| `natural_product` | `0`×780、`1`×2 |
| `first_in_class` | `-1`×777、`0`×5 |
| `inorganic_flag` | `-1`×777、`0`×5 |
| `therapeutic_flag` | `0`×782 |
| `availability_type` | `-1`×777、`(空)`×5 |
| `chemical_probe` | `0`×781、`1`×1 |
| `orphan` | `-1`×777、`0`×5 |

## 五、盐型归并与同义词

`molecule_hierarchy` 里 **781** 个分子本身就是母体，**1** 个是盐型/衍生记录（`is_parent = FALSE`，`parent_chembl_id` 给出母体）。

| 分子 | 母体 |
| --- | --- |
| `CHEMBL1204008` | `CHEMBL575092` |

`molecule_synonyms` 只覆盖 **5 / 782** 个分子——有名字的基本就是进过临床的那几个，其余是文献化合物只有 ChEMBL ID。

| 分子 | 名称 / 研发代号 |
| --- | --- |
| `CHEMBL2165620` | Pf-04991532（OTHER） |
| `CHEMBL3219124` | Azd 1656（OTHER）、Azd-1656（OTHER）、Azd1656（OTHER）、AZD1656（RESEARCH_CODE） |
| `CHEMBL1783734` | Piragliatin（INN）、Piragliatine（INN_FRENCH）、Piragliatina（INN_SPANISH）、R04389620（RESEARCH_CODE）、R04389620\（RESEARCH_CODE）、R |
| `CHEMBL2165615` | Nerigliatin（INN）、Nerigliatine（INN_FRENCH）、Nerigliatina（INN_SPANISH）、Nerigliatin（OTHER）、Pf-04937319（OTHER）、PF-04937319（RE |
| `CHEMBL5072532` | BMS-820132（PND） |

## 附：按 priority 的性质概览

| priority | 分子数 | MW 中位 | ALogP 中位 | PSA 中位 | 六项全落窗口 |
| --- | ---: | ---: | ---: | ---: | ---: |
| P1 | 39 | 490.5 | 4.54 | 105.2 | 3 |
| P2 | 83 | 458.5 | 4.35 | 97.1 | 10 |
| P3 | 202 | 470.3 | 3.98 | 109.8 | 21 |
| P4 | 458 | 458.5 | 4.16 | 105.4 | 51 |

