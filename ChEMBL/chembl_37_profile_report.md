# ChEMBL 数据库结构与内容概览报告

> 自动生成于 2026-07-30 16:14　|　数据库文件：`/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db`
>
> 本报告面向**有生物学 / 生信背景、但没有药物化学背景**的读者。每一节先解释「这是什么」，再给出这个数据库里的实际统计。


## 0. 三十秒速览


**ChEMBL 是什么？** 一个由 EMBL-EBI 维护的人工整理（manually curated）的生物活性数据库。它把散落在几万篇论文和专利里的一句句话——「化合物 X 抑制蛋白 Y，IC50 = 12 nM」——抽取成结构化的表格。你可以把它理解成**化合物-靶点相互作用领域的 UniProt + GEO**：

- 像 UniProt 一样，它有稳定 ID（`CHEMBL25` 就是阿司匹林）、人工审编、版本化发布；
- 像 GEO 一样，它的一条记录必须绑定「在什么实验条件下测的」，脱离实验条件的数值没有意义。

**一句话的数据模型**：

```
某个化合物  在某个实验里  针对某个靶点  测到了某个数值  出处是某篇文献
(compound)    (assay)       (target)      (activity)      (document)
```

这五个实体就是 ChEMBL 的骨架，其余 60 多张表都是围绕它们的注释、字典和映射。

| 条目 | 值 |
| --- | --- |
| 数据库版本 | CHEMBL_37 |
| 版本创建日期 | 2026-05-01 00:00:00.000000 |
| 版本说明 | ChEMBL Release 37 (https://www.ebi.ac.uk/chembl) |
| 文件大小 | 30.48 GB |
| 表数量 | 72 |
| 化合物数（molecule_dictionary） | 2,921,148 |
| 活性数据点数（activities） | 24,527,044 |
| 实验数（assays） | 1,970,438 |
| 靶点数（target_dictionary） | 18,552 |
| 来源文献数（docs） | 101,100 |
| SQLite 引擎（运行本脚本的） | 3.31.1 |


### 0.1 本版所依赖的外部资源版本


ChEMBL 的注释是**站在其他数据库肩膀上**的。下表来自 `version` 表，记录了本次发布所用的外部资源版本。做可重复性分析或跨库整合时，这决定了你的 UniProt / GO / MeSH 注释应该对齐到哪个版本。

| 资源 | 用途 |
| --- | --- |
| **Bioassay Ontology 2.0** | BAO version used for assays (http://bioassayontology.org/) |
| **COCONUT 2025-07** | COCONUT version used for natural product flagging (https://coconut.naturalproducts.net/) |
| **ChEMBL_Structure_Pipeline 1.2.0** | ChEMBL_Structure_Pipeline version used for chemical structure curation (https://github.com/chembl/ChEMBL_Structure_Pipeline) |
| **EFO 3.74.0** | EFO version used for indication and warning data (https://www.ebi.ac.uk/efo/) |
| **Gene Ontology 2024-02-22** | GO version used for genes (https://www.ebi.ac.uk/QuickGO/) |
| **InChI v1.06** | InChI version used for compound registration (https://github.com/IUPAC-InChI/InChI) |
| **MeSH 2025** | MeSH version used for indication data (https://www.nlm.nih.gov/mesh/meshhome.html) |
| **RDKit 2022.09.4** | RDKit version used for standardisation and salt stripping of chemical structures (https://www.rdkit.org/) |
| **Swiss-Prot 2025_03** | UniProtKB version used for component_sequences data (https://www.expasy.org/resources/uniprotkb) |
| **UBERON 2024-03-22** | UBERON version used for tissues (http://obophenotype.github.io/uberon/) |


### 0.2 最近的发布历史


> ChEMBL 大约每 6–12 个月发一版。**每一版都会回溯性地修订旧数据**（合并重复化合物、重新指认靶点、调整临床阶段标注），所以任何分析都必须写明用的是哪一版。

| 版本 | 创建日期 |
| --- | --- |
| CHEMBL_37 | 2026-05-01 |
| CHEMBL_36 | 2025-07-28 |
| CHEMBL_35 | 2024-12-01 |
| CHEMBL_34 | 2024-03-28 |
| CHEMBL_33 | 2023-05-31 |
| CHEMBL_32 | 2023-01-26 |
| CHEMBL_31 | 2022-07-12 |
| CHEMBL_30 | 2022-02-22 |


### 0.3 所有 CHEMBLxxxxx 标识符的构成


拿到一个陌生的 `CHEMBLxxxxx`，可以在 `chembl_id_lookup` 里查它到底是什么实体。

| entity_type | 数量 | 是什么 |
| --- | ---: | --- |
| `COMPOUND` | 3,082,236 | 化合物 |
| `ASSAY` | 2,271,683 | 实验 |
| `DOCUMENT` | 101,436 | 文献/数据集 |
| `TARGET` | 20,483 | 靶点 |
| `CELL` | 2,312 | 细胞系 |
| `TISSUE` | 802 | 组织 |


## 1. 先搞懂这些词


如果你是从生物学过来的，下面这几个概念是读懂 ChEMBL 的全部前提。


### 1.1 化合物（compound / molecule）


一个小分子药物或研究用化合物。它在数据库里的身份证有两套：

- **`molregno`** — 内部整数主键，所有表之间用它连接（相当于 NCBI 的 GeneID）；
- **`chembl_id`** — 对外稳定 ID，形如 `CHEMBL25`（相当于基因的 Ensembl ID）。

化学结构本身以文本形式存放，最常用的是 **SMILES**：一行 ASCII 字符串就编码了一个分子的原子连接关系，例如阿司匹林是 `CC(=O)Oc1ccccc1C(=O)O`。你不需要会读 SMILES —— RDKit 之类的库可以把它解析成分子对象、算描述符、算相似度。另外还有 **InChIKey**，一个 27 位的定长哈希（如 `BSYNRYMUTXBXSQ-UHFFFAOYSA-N`），作用完全等同于序列的 MD5：用来跨数据库精确匹配同一个分子。


### 1.2 靶点（target）


化合物作用的对象。**注意它不总是一个蛋白**：可能是单一蛋白、蛋白复合物、一条蛋白家族、一个细胞系，甚至一整个物种（抗菌实验里靶点就是「大肠杆菌」）。`target_dictionary.target_type` 记录了是哪一种。

要从**基因名**找到靶点，路径是：`component_synonyms`（基因名/蛋白别名）→ `component_sequences`（UniProt accession）→ `target_components` → `target_dictionary`。`component_sequences.accession` 就是 UniProt 号，这是 ChEMBL 与你熟悉的生物学世界之间最重要的一座桥。


### 1.3 实验（assay）


一次具体的测量设定：用什么体系、测什么读数、针对哪个靶点。**同一个化合物-靶点对，在不同实验里测出的数值可能差一个数量级**，这不是数据错误，而是实验条件不同（ATP 浓度、细胞类型、孵育时间…）。所以 ChEMBL 从不把活性直接挂在「化合物-靶点」上，而是必须经过 assay 这一层。这一点和 GEO 里样本必须绑定实验设计是一个道理。

两个字段决定了这条数据能不能用：

- **`assay_type`** — 实验大类（见 §5）；
- **`confidence_score`（0–9）** — 「靶点指认可信度」：ChEMBL 的编者对『这个实验真的测的是这个靶点吗』的信心打分。做定量建模一般要求 **≥ 8**（单一蛋白，明确）。


### 1.4 活性（activity）——最容易踩坑的一张表


`activities` 是核心事实表，一行 = 一个数值。关键在于它有**两套值**：

| 字段 | 含义 |
| --- | --- |
| `type` / `value` / `units` / `relation` | **原文照抄**：论文里怎么写就怎么存（可能是 `Ic-50`、单位 µM） |
| `standard_type` / `standard_value` / `standard_units` / `standard_relation` | **标准化后**：类型名统一（`IC50`）、单位统一（浓度一律 nM） |

**永远用 `standard_*` 这一套做分析**，`type/value` 只在追溯原文时用。

`standard_relation` 是个容易被忽略的坑：它可能是 `=`，也可能是 `>` 或 `<`。`> 10000 nM` 的意思是「测到最高浓度也没打到 50%，实际值未知，只知道比这大」——这是**删失数据（censored data）**，直接当成 10000 会污染回归模型。生存分析里的删失是同一个概念。


### 1.5 pChEMBL —— 你真正应该用的那个数


`pchembl_value = −log10(标准化浓度，单位 M)`，只对浓度型指标（IC50/EC50/XC50/AC50/Ki/Kd/Potency）且 `relation = '='` 时计算。

- IC50 = 1 nM  → pChEMBL = 9
- IC50 = 100 nM → pChEMBL = 7
- IC50 = 10 µM  → pChEMBL = 5

为什么要取负对数？因为**结合自由能与浓度的对数成正比**（ΔG = −RT·ln K），所以 pChEMBL 才是那个物理上线性、分布近似正态、适合做回归的量。这跟你对表达量取 log2 是完全一样的动机。**数值越大 = 活性越强**（方向和 IC50 相反，注意别搞反）。

经验尺度：pChEMBL ≥ 6（IC50 ≤ 1 µM）算「有活性」，≥ 8（≤ 10 nM）算「强效」。


### 1.6 max_phase —— 这个化合物走到哪一步了


`molecule_dictionary.max_phase` 记录化合物达到过的最高临床阶段：4 = 已上市药物，3/2/1 = 对应临床期，0.5 = 早期，0 或空 = 纯研究化合物。想快速捞出「所有已知药物」，条件就是 `max_phase = 4`。

> 注意：不同 ChEMBL 版本之间 `max_phase` 会变动（官方在持续复核「什么算获批」），跨版本比较药物数量时务必注明版本号。


## 2. 数据模型：表是怎么连起来的


```
                    compound_properties   compound_structures
                          │ molregno            │ molregno
                          └──────┬──────────────┘
                                 │
                        molecule_dictionary          docs ──────┐
                                 │ molregno       (文献/专利)   │ doc_id
                                 │                              │
        compound_records ────────┤                              │
         (文献级记录)  record_id │                              │
                                 │                              │
                          ┌──────┴──────────────────────────────┴──┐
                          │            ACTIVITIES                  │  ← 核心事实表
                          │  molregno / assay_id / doc_id / src_id │
                          └──────────────────┬─────────────────────┘
                                             │ assay_id
                                          ASSAYS
                                     (实验设定 + 可信度)
                                             │ tid
                                     target_dictionary
                                             │ tid
                                      target_components
                                             │ component_id
                              component_sequences  (UniProt accession)
                                             │
                        component_synonyms / component_go / component_class
```

**读法**：从 `activities` 出发，向左连到「测的是什么分子」，向下连到「怎么测的、测的什么靶点」，向上连到「出处是哪篇文献」。几乎所有实用查询都是这个骨架的变形。


### 2.1 核心表的外键（从数据库实际读取）


| 表 | 外键列 | 指向表 | 指向列 |
| --- | --- | --- | --- |
| `activities` | `data_validity_comment` | `data_validity_lookup` | `data_validity_comment` |
| `activities` | `src_id` | `source` | `src_id` |
| `activities` | `record_id` | `compound_records` | `record_id` |
| `activities` | `molregno` | `molecule_dictionary` | `molregno` |
| `activities` | `doc_id` | `docs` | `doc_id` |
| `activities` | `bao_endpoint` | `bioassay_ontology` | `bao_id` |
| `activities` | `assay_id` | `assays` | `assay_id` |
| `activities` | `action_type` | `action_type` | `action_type` |
| `assays` | `bao_format` | `bioassay_ontology` | `bao_id` |
| `assays` | `variant_id` | `variant_sequences` | `variant_id` |
| `assays` | `tissue_id` | `tissue_dictionary` | `tissue_id` |
| `assays` | `tid` | `target_dictionary` | `tid` |
| `assays` | `src_id` | `source` | `src_id` |
| `assays` | `relationship_type` | `relationship_type` | `relationship_type` |
| `assays` | `doc_id` | `docs` | `doc_id` |
| `assays` | `confidence_score` | `confidence_score_lookup` | `confidence_score` |
| `assays` | `chembl_id` | `chembl_id_lookup` | `chembl_id` |
| `assays` | `cell_id` | `cell_dictionary` | `cell_id` |
| `assays` | `assay_type` | `assay_type` | `assay_type` |
| `target_dictionary` | `target_type` | `target_type` | `target_type` |
| `target_dictionary` | `chembl_id` | `chembl_id_lookup` | `chembl_id` |
| `target_components` | `tid` | `target_dictionary` | `tid` |
| `target_components` | `component_id` | `component_sequences` | `component_id` |
| `molecule_dictionary` | `chembl_id` | `chembl_id_lookup` | `chembl_id` |
| `compound_records` | `src_id` | `source` | `src_id` |
| `compound_records` | `molregno` | `molecule_dictionary` | `molregno` |
| `compound_records` | `doc_id` | `docs` | `doc_id` |
| `compound_structures` | `molregno` | `molecule_dictionary` | `molregno` |


### 2.2 主键命名约定


| 主键 | 属于 | 说明 |
| --- | --- | --- |
| `molregno` | 化合物 | molecule registration number |
| `tid` | 靶点 | target id |
| `assay_id` | 实验 |  |
| `activity_id` | 活性数据点 |  |
| `doc_id` | 文献 |  |
| `record_id` | 文献级化合物记录 | 同一化合物出现在 N 篇文献 → N 个 record_id，但只有 1 个 molregno |
| `component_id` | 靶点组件（蛋白） |  |
| `src_id` | 数据来源 |  |

> `chembl_id_lookup` 表可以反查任意 `CHEMBLxxxxx` 属于哪一类实体——拿到一个陌生 ChEMBL ID 时先查它。


## 3. 表清单与规模


共 **72** 张表，按行数排序（精确计数）。「说明」一栏是本脚本内置的中文注解。

| 表名 | 行数 | 说明 |
| --- | ---: | --- |
| `activities` | 24,527,044 | 核心事实表：一次测量得到的一个活性数值（如某化合物对某靶点的 IC50） |
| `activity_properties` | 12,213,211 | 某条活性记录的附加参数（如测定时的底物浓度、pH） |
| `chembl_id_lookup` | 5,478,952 | 全局 ChEMBL ID 索引：给一个 CHEMBLxxxx，告诉你它是化合物/靶点/文献/实验 |
| `compound_structural_alerts` | 5,020,133 | 结构警示：命中已知易假阳性/易毒性的子结构 |
| `compound_records` | 3,824,604 | 文献级记录：同一化合物在不同文献里出现一次就有一行 |
| `molecule_dictionary` | 2,921,148 | 化合物主表：一个 molregno = 一个唯一化合物（含名称、研发阶段等） |
| `compound_properties` | 2,901,464 | 计算得到的理化性质（分子量、logP、氢键供受体数等） |
| `compound_structures` | 2,897,819 | 化学结构：SMILES / InChI / InChIKey / molfile |
| `molecule_hierarchy` | 2,828,129 | 化合物家族关系：盐/前药 → 母体活性分子 |
| `ligand_eff` | 2,223,169 | 配体效率指标（活性按分子大小归一化） |
| `activity_supp_map` | 2,010,125 | activities 与 activity_supp 的映射 |
| `assays` | 1,970,438 | 实验（assay）描述：这次测量是怎么做的、测的哪个靶点、可信度多高 |
| `activity_supp` | 1,776,415 | 存放补充/原始的多维测定数据（深度数据，通常不需要） |
| `activity_smid` | 1,732,478 |  |
| `predicted_binding_domains` | 822,313 | 预测的结合结构域 |
| `assay_parameters` | 460,048 | 实验的参数（如给药剂量、给药途径、动物品系） |
| `assay_class_map` | 244,490 | assays 与 assay_classification 的多对多映射 |
| `target_relations` | 155,208 | 靶点之间的关系（如亚基属于复合物、复合物属于家族） |
| `component_go` | 153,709 | 组件 → GO 注释的映射 |
| `molecule_synonyms` | 136,061 | 化合物别名（商品名、代号、INN 名等） |
| `component_synonyms` | 121,157 | 组件别名（基因名、蛋白名等）— 用基因名找靶点时从这里入手 |
| `docs` | 101,100 | 数据来源文献/专利/数据集（PubMed ID、DOI、期刊、年份） |
| `indication_refs` | 93,733 | drug_indication 的文献出处 |
| `drug_indication` | 60,055 | 药物适应症（映射到 MeSH / EFO 疾病本体） |
| `formulations` | 53,430 | 药品制剂 → 所含化合物的映射 |
| `products` | 45,752 | FDA 批准的药品（制剂层面，不是分子层面） |
| `component_domains` | 28,547 | 组件 → Pfam 结构域的映射 |
| `biotherapeutics` | 23,849 | 生物药（抗体、多肽等）的额外信息 |
| `product_patents` | 19,705 | 药品对应的专利信息 |
| `target_dictionary` | 18,552 | 靶点主表：一个 tid = 一个作用对象（蛋白、蛋白复合物、细胞系、整个生物体…） |
| `target_components` | 17,284 | 靶点 → 组成它的分子组件（蛋白）的映射 |
| `mechanism_refs` | 13,600 | drug_mechanism 的文献出处 |
| `component_class` | 13,180 | 组件 → 蛋白家族分类的映射 |
| `component_sequences` | 12,986 | 靶点组件的序列信息（UniProt accession、序列、物种） |
| `drug_mechanism` | 7,561 | 已知药物的作用机制：药物 × 靶点 × 作用类型（激动/抑制…） |
| `site_components` | 5,908 | 结合位点由哪些组件构成 |
| `atc_classification` | 5,579 | WHO ATC 药物分类字典（按解剖-治疗-化学分级） |
| `warning_refs` | 4,946 | drug_warning 的文献出处 |
| `molecule_atc_classification` | 4,567 | 化合物 → WHO ATC 药物分类的映射 |
| `binding_sites` | 4,545 | 结合位点定义 |
| `biotherapeutic_components` | 4,415 | 生物药 → 组成序列的映射 |
| `organism_class` | 4,280 | 物种分类字典（对齐 NCBI taxonomy） |
| `patent_use_codes` | 4,043 | 专利用途代码字典 |
| `domains` | 3,953 | Pfam 结构域字典 |
| `bio_component_sequences` | 3,478 | 生物药组分的序列 |
| `metabolism_refs` | 3,296 | metabolism 的文献出处 |
| `variant_sequences` | 2,836 | 突变体序列（如激酶耐药突变） |
| `defined_daily_dose` | 2,721 | WHO 定义的日剂量 |
| `drug_warning` | 2,304 | 药物安全性警示（黑框警告、撤市等） |
| `cell_dictionary` | 2,238 | 细胞系字典 |
| `metabolism` | 2,147 | 药物代谢途径（底物 → 代谢产物 → 催化酶） |
| `structural_alerts` | 936 | 结构警示子结构定义 |
| `protein_classification` | 905 | ChEMBL 自建的蛋白家族树（激酶/GPCR/离子通道…） |
| `usan_stems` | 834 | USAN 药物命名词干（如 -tinib 表示激酶抑制剂） |
| `tissue_dictionary` | 791 | 组织字典（对齐 UBERON 本体） |
| `pesticide_classification` | 595 | 农药分类字典 |
| `pesticide_class_mapping` | 593 | 化合物 → 农药分类映射 |
| `assay_classification` | 584 | 体内实验的疾病/治疗领域分类 |
| `bioassay_ontology` | 311 | BioAssay Ontology (BAO) 术语字典 |
| `go_classification` | 309 | GO 分类字典 |
| `sqlite_stat1` | 210 | SQLite 内部统计表（查询优化器用，不是 ChEMBL 数据） |
| `activity_stds_lookup` | 151 | 标准化字典：哪些活性类型允许被标准化成哪些单位 |
| `source` | 67 | 数据源字典：这条数据来自文献、PubChem、BindingDB 还是某个捐赠数据集 |
| `chembl_release` | 37 | 本数据库的版本与发布日期 |
| `action_type` | 35 | 作用类型字典（激动剂、抑制剂…）及其上位归类 |
| `target_type` | 28 | 靶点类型字典（SINGLE PROTEIN / PROTEIN COMPLEX / ORGANISM / CELL-LINE…） |
| `version` | 11 | 版本信息（旧字段） |
| `confidence_score_lookup` | 10 | 靶点指认可信度 0–9 的含义字典 |
| `data_validity_lookup` | 7 | 数据可疑标记的含义字典 |
| `assay_type` | 6 | 实验大类字典：B=结合、F=功能、A=ADME、T=毒性、P=理化、U=未分类 |
| `relationship_type` | 6 | assay 与靶点关系类型字典 |
| `structural_alert_sets` | 5 | 结构警示规则集（PAINS、Dundee 等） |


## 4. 数据从哪来


ChEMBL 不是单一来源。`source` 表列出了所有数据源，`activities.src_id` 指向它。了解来源很重要，因为**不同来源的数据质量与稠密度差别很大**：文献数据经过人工审编但稀疏，高通量筛选数据量大但多为单浓度、阴性居多。

| src_id | 简称 | 描述 | activities 数 |  |
| ---: | --- | --- | ---: | --- |
| 1 | `LITERATURE` | Scientific Literature | 9,554,640 | ████████████████████ |
| 7 | `PUBCHEM_BIOASSAY` | PubChem BioAssays | 7,434,992 | ████████████████ |
| 37 | `BINDINGDB` | BindingDB Patent Bioactivity Data | 2,682,137 | ██████ |
| 33 | `GATES_LIBRARY` | Gates Foundation Compound Collection | 1,482,491 | ███ |
| 55 | `EUBOPEN_CGL` | EUbOPEN Chemogenomic Library | 755,560 | ██ |
| 15 | `DRUGMATRIX` | DrugMatrix | 494,046 | █ |
| 51 | `WINZ_PLASMO` | Winzeler Lab Plasmodium Screening | 399,067 | █ |
| 60 | `MMV_MALARIA_HGL` | MMV Malaria Hit Generation Library | 308,377 | █ |
| 11 | `TG_GATES` | Open TG-GATEs | 210,708 |  |
| 38 | `PATENT` | SureChEMBL Patent Bioactivity Data | 180,540 |  |
| 16 | `GSK_PKIS` | GSK Published Kinase Inhibitor Set | 169,451 |  |
| 40 | `COADD` | CO-ADD Antimicrobial Screening | 99,793 |  |
| 54 | `DONATED_PROBES` | SGC Frankfurt - Donated Chemical Probes | 83,959 |  |
| 2 | `GSK_TCMDC` | GSK Malaria Screening | 81,198 |  |
| 72 | `LIT_CHEM_PROBES` | Chemical Probe data from Scientific Literature | 76,924 |  |
| 5 | `SANGER` | Sanger Institute Genomics of Drug Sensitivity in Cancer | 73,169 |  |
| 48 | `TUM_PROTEOMIC_KUSTER` | Kuster Lab Chemical Proteomics Drug Profiling | 70,505 |  |
| 69 | `ZIMM_BT_12_23` | EMBL Heidelberg Gut Microbiome Host Interactions | 47,351 |  |
| 17 | `MMV_MBOX` | Medicines for Malaria Venture (MMV) Malaria Box | 45,158 |  |
| 32 | `ST_JUDE_LEISH` | St. Jude Children’s Research Hospital Leishmania Screening | 42,105 |  |
| 52 | `SARS_COV_2` | SARS-CoV-2 Screening Data | 37,209 |  |
| 71 | `ASAP` | AI-driven Structure-enabled Antiviral Platform (ASAP) | 30,781 |  |
| 3 | `NOVARTIS` | Novartis Malaria Screening | 27,888 |  |
| 65 | `LIT_EUBOPEN_CGL` | EUbOPEN Chemogenomic Library Literature Data | 23,222 |  |
| 61 | `KI_EUBOPEN` | Karolinska Institute dNTPase SAMHD1 screening | 17,834 |  |
| 14 | `DNDI` | Drugs for Neglected Diseases Initiative (DNDi) | 14,452 |  |
| 27 | `ASTRAZENECA` | AstraZeneca DMPK/physicochemical | 11,687 |  |
| 59 | `HDAC6` | Fraunhofer Institute HDAC6 screening | 11,680 |  |
| 57 | `CARE` | IMI-CARE SARS-CoV-2 Data | 9,646 |  |
| 29 | `GSK_TCAKS` | GSK Kinetoplastid Screening | 7,235 |  |
| 18 | `TP_TRANSPORTER` | TP-search Transporter Database | 6,765 |  |
| 20 | `WHO_TDR` | WHO Tropical Disease Research (TDR) Malaria Screening | 5,853 |  |
| 4 | `ST_JUDE` | St. Jude Children’s Research Hospital Malaria Screening | 5,456 |  |
| 34 | `MMV_PBOX` | Medicines for Malaria Venture (MMV) Pathogen Box | 5,056 |  |
| 21 | `SUPPLEMENTARY` | Deposited Supplementary Bioactivity Data | 4,817 |  |
| 67 | `DUNDEE_T_CRUZI` | University of Dundee T. Cruzi | 3,328 |  |
| 30 | `K4DD` | Kinetics for Drug Discovery (K4DD) | 2,064 |  |
| 22 | `GSK_TB` | GSK Tuberculosis Screening | 1,814 |  |
| 68 | `EU-OPENSCREEN` | EU-OPENSCREEN | 1,813 |  |
| 28 | `FDA_APPROVAL` | FDA approval pharmacokinetics/metabolism  | 1,387 |  |
| 56 | `SALVENSIS_LSHTM` | London School of Hygiene and Tropical Medicine and Salvensis Schistosomiasis screening | 1,222 |  |
| 39 | `DRUG_PK` | Curated Drug Pharmacokinetic Data | 1,163 |  |
| 49 | `HESI` | HESI Cardiac Safety Committee Myocyte Subteam dataset | 986 |  |
| 43 | `PKIS2` | GSK Published Kinase Inhibitor Set 2 | 491 |  |
| 23 | `OSM` | Open Source Malaria Screening | 344 |  |
| 70 | `ABER_NTD` | Aberystwyth University Schistosomiasis | 268 |  |
| 64 | `CSD23` | Cardiff University Schistosomiasis | 194 |  |
| 19 | `HARVARD` | Harvard University Malaria Screening | 111 |  |
| 58 | `RESOLUTE` | Research Empowerment on Solute Carriers (RESOLUTE) | 96 |  |
| 31 | `METABOLISM` | Curated Drug Metabolism Pathways | 11 |  |
| 0 | `UNDEFINED` | Undefined | 0 |  |
| 6 | `PDBE` | PDBe Ligands (DEPRECATED) | 0 |  |
| 8 | `CANDIDATES` | Clinical Candidate Compounds | 0 |  |
| 9 | `FDA_ORANGE_BOOK` | FDA Orange Book Drugs | 0 |  |
| 10 | `GRAC` | Guide to Receptors and Channels (DEPRECATED) | 0 |  |
| 12 | `FDA_NEW_DRUGS` | FDA Novel Drugs and Biotherapeutics | 0 |  |
| 13 | `USAN` | United States Adopted Names (USAN) | 0 |  |
| 24 | `MILLIPORE` | Millipore Kinase Screening (DEPRECATED - MERGED WITH SRC_ID = 1) | 0 |  |
| 25 | `EXT. PROJECT CPDS` | External Project Compounds | 0 |  |
| 26 | `ATLAS` | Gene Expression Atlas Compounds (EMBL-EBI) | 0 |  |
| 35 | `HECATOS` | Hepatic and Cardiac Toxicity Systems modelling (HeCaToS) Compounds | 0 |  |
| 36 | `WITHDRAWN` | Withdrawn Drugs | 0 |  |
| 41 | `ATC` | WHO Anatomical Therapeutic Chemical (ATC) Classification of Drugs | 0 |  |
| 42 | `BNF` | British National Formulary (BNF) | 0 |  |
| 53 | `PRODRUG_ACTIVE` | Active Ingredient of a Prodrug | 0 |  |
| 63 | `INN` | International Nonproprietary Names (INN) for Pharmaceutical Substances | 0 |  |
| 66 | `EMA` | European Medicines Agency (EMA) | 0 |  |


## 5. 化合物层面


- 唯一化合物：**2,921,148**
- 其中有化学结构（SMILES 等）：**2,897,819**（99.2%）—— 没有结构的多为抗体、细胞疗法等大分子
- 有计算理化性质：**2,901,464**（99.3%）


### 5.1 化合物类型

| molecule_type | 数量 | 占比 |  | 说明 |
| --- | ---: | ---: | --- | --- |
| Small molecule | 1,920,259 | 65.7% | ██████████████████ | 小分子（绝大多数传统药物与研究化合物） |
| (空) | 571,492 | 19.6% | █████ |  |
| Unknown | 404,621 | 13.9% | ████ | 未标注 |
| Protein | 22,799 | 0.8% |  | 蛋白类（含酶、融合蛋白） |
| Antibody | 1,032 | 0.0% |  | 抗体 |
| Oligonucleotide | 260 | 0.0% |  | 寡核苷酸（ASO、siRNA 等） |
| Gene | 191 | 0.0% |  | 基因疗法 |
| Enzyme | 129 | 0.0% |  | 酶 |
| Antibody drug conjugate | 109 | 0.0% |  | 抗体偶联药物（抗体连上小分子毒素） |
| Vaccine component | 90 | 0.0% |  | 疫苗组分 |
| Cell | 85 | 0.0% |  | 细胞疗法 |
| Oligosaccharide | 81 | 0.0% |  | 寡糖 |


### 5.2 研发阶段 max_phase

| max_phase | 化合物数 | 含义 |
| --- | ---: | --- |
| None | 2,901,840 | 未标注 |
| 2 | 9,054 | II 期临床 |
| 4 | 4,225 | 已获批上市药物 |
| -1 | 2,499 | 曾出现在临床相关来源但阶段未知 |
| 3 | 1,892 | 处于/已完成 III 期临床 |
| 1 | 1,613 | I 期临床 |
| 0.5 | 25 | 临床前 / 早期临床（ChEMBL 新增的中间档） |

> 想拿「已上市药物」集合：`SELECT * FROM molecule_dictionary WHERE max_phase = 4`，本库共 **4,225** 个。


### 5.3 化合物标记位（值为 1 的数量）

| 字段 | 数量 | 含义 |
| --- | ---: | --- |
| `therapeutic_flag` | 3,890 | 被标记为治疗用药 |
| `natural_product` | 98,989 | 天然产物来源 |
| `oral` | 1,979 | 可口服 |
| `parenteral` | 1,573 | 注射给药 |
| `topical` | 601 | 外用 |
| `black_box_warning` | 936 | 有黑框警告（严重安全性问题） |
| `withdrawn_flag` | 364 | 已撤市 |
| `prodrug` | 487 | 前药（体内代谢后才有活性） |
| `chemical_probe` | 807 | 化学探针（工具分子，选择性经过验证） |
| `inorganic_flag` | 204 | 无机物 |
| `polymer_flag` | 261 | 聚合物 |


### 5.4 计算理化性质的分布


这些是**用软件从结构算出来的**（不是实验测的），常用于快速过滤化合物。对非化学背景读者，最需要知道的是所谓 **Lipinski 类药五规则（Rule of Five）**：分子量 ≤ 500、logP ≤ 5、氢键供体 ≤ 5、氢键受体 ≤ 10 —— 满足这些的分子更可能有良好口服吸收。`num_ro5_violations` 就是违反了几条。

| 字段 | 有值数 | 最小 | 中位 | 均值 | 最大 | 单位 | 含义 |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `mw_freebase` | 2,901,461 | 4.00 | 402.45 | 441.35 | 13,170.65 | Da | 分子量（游离碱形式） |
| `alogp` | 2,824,312 | -14.26 | 3.50 | 3.53 | 22.57 |  | 计算 logP，亲脂性；越大越亲脂 |
| `hba` | 2,824,312 | 0.00 | 5.00 | 5.45 | 32.00 |  | 氢键受体数 |
| `hbd` | 2,824,312 | 0.00 | 1.00 | 1.63 | 25.00 |  | 氢键供体数 |
| `psa` | 2,824,312 | 0.00 | 78.27 | 84.47 | 595.22 | Å² | 极性表面积；与膜通透性、口服吸收相关 |
| `rtb` | 2,824,312 | 0.00 | 5.00 | 5.81 | 67.00 |  | 可旋转键数；反映分子柔性 |
| `num_ro5_violations` | 2,824,312 | 0.00 | 0.00 | 0.44 | 4.00 |  | 违反 Lipinski 五规则的条数 |
| `aromatic_rings` | 2,824,312 | 0.00 | 3.00 | 2.58 | 30.00 |  | 芳香环数 |
| `heavy_atoms` | 2,824,312 | 1.00 | 28.00 | 29.24 | 79.00 |  | 重原子（非氢原子）数 |
| `qed_weighted` | 2,824,312 | 0.01 | 0.55 | 0.54 | 0.95 |  | QED 类药性综合评分，0–1，越大越「像药」 |
| `np_likeness_score` | 2,824,312 | -4.13 | -1.05 | -0.88 | 4.13 |  | 天然产物相似度评分 |


## 6. 靶点层面


靶点总数：**18,552**。再强调一次：靶点不等于蛋白，见下表分布。


### 6.1 靶点类型

| target_type | 数量 |  | 说明 |
| --- | ---: | --- | --- |
| SINGLE PROTEIN | 11,055 | ██████████████████ | 单一蛋白 — 可直接对应一个 UniProt / 基因 |
| ORGANISM | 2,784 | █████ | 整个生物体（如某种细菌、寄生虫） |
| CELL-LINE | 2,000 | ███ | 细胞系整体（表型筛选） |
| PROTEIN-PROTEIN INTERACTION | 804 | █ | 蛋白-蛋白相互作用界面（如分子胶/降解剂的作用对象） |
| PROTEIN COMPLEX | 647 | █ | 蛋白复合物（多亚基共同构成作用对象） |
| PROTEIN FAMILY | 428 | █ | 蛋白家族（未细分到具体成员） |
| TISSUE | 293 |  | 组织 |
| NUCLEIC-ACID | 129 |  | 核酸 |
| SELECTIVITY GROUP | 123 |  | 选择性分组（用于比较同一化合物对一组相关靶点的选择性） |
| PROTEIN COMPLEX GROUP | 67 |  | 复合物家族 |
| SMALL MOLECULE | 49 |  | 小分子作为作用对象（如螯合剂） |
| CHIMERIC PROTEIN | 35 |  | 嵌合蛋白 |
| UNKNOWN | 24 |  | 未知 |
| OLIGOSACCHARIDE | 22 |  | 寡糖 |
| SUBCELLULAR | 20 |  | 亚细胞结构 |
| MACROMOLECULE | 19 |  | 大分子 |
| PROTEIN NUCLEIC-ACID COMPLEX | 15 |  | 蛋白-核酸复合物 |
| LIPID | 11 |  | 脂质 |
| 3D CELL CULTURE | 11 |  | 三维细胞培养体系（类器官等） |
| METAL | 10 |  | 金属离子 |
| PHENOTYPE | 2 |  | 表型 |
| UNCHECKED | 1 |  | 未审核 |
| NON-MOLECULAR | 1 |  | 非分子实体 |
| NO TARGET | 1 |  | 无靶点信息 |
| ADMET | 1 |  | ADMET 性质 |


### 6.2 靶点物种分布（Top 20）

| 物种 | 靶点数 |
| --- | ---: |
| Homo sapiens | 9,118 |
| Mus musculus | 1,710 |
| Rattus norvegicus | 1,347 |
| (空) | 275 |
| Bos taurus | 248 |
| Sus scrofa | 134 |
| Oryctolagus cuniculus | 109 |
| Canis lupus familiaris | 101 |
| Escherichia coli (strain K12) | 93 |
| Mycobacterium tuberculosis | 89 |
| Cavia porcellus | 73 |
| Saccharomyces cerevisiae S288c | 72 |
| Escherichia coli K-12 | 71 |
| Escherichia coli | 69 |
| Staphylococcus aureus | 64 |
| Plasmodium falciparum | 56 |
| Gallus gallus | 48 |
| Pseudomonas aeruginosa | 45 |
| Macaca mulatta | 45 |
| Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv) | 43 |

> **人源单一蛋白靶点**共 **5,869** 个 —— 这通常是做人类靶点分析时的起点子集。


### 6.3 与 UniProt / 基因的对接

- `component_sequences` 共 12,986 行，其中 **12,959 个不同的 UniProt accession**。
- 组件类型：PROTEIN 12,912，RNA 69，DNA 5

从基因名查靶点的标准写法：

```sql
SELECT DISTINCT td.chembl_id, td.pref_name, td.organism, td.target_type
FROM component_synonyms cs
JOIN component_sequences seq ON cs.component_id = seq.component_id
JOIN target_components tc    ON tc.component_id = seq.component_id
JOIN target_dictionary td    ON td.tid = tc.tid
WHERE cs.component_synonym = 'GCK'        -- 基因名
  AND cs.syn_type = 'GENE_SYMBOL';
```


### 6.4 蛋白家族分类（Top 15，按组件数）

| 蛋白家族 | 组件数 |
| --- | ---: |
| Unclassified protein | 2,497 |
| Hydrolase | 1,343 |
| Transferase | 1,268 |
| Oxidoreductase | 899 |
| Other cytosolic protein | 631 |
| Enzyme | 496 |
| Transcription factor | 253 |
| Secreted protein | 250 |
| Lyase | 230 |
| Membrane receptor | 207 |
| Ligase | 193 |
| Isomerase | 185 |
| Structural protein | 147 |
| Serine protease S1A subfamily | 96 |
| Voltage-gated potassium channel | 93 |


## 7. 实验（assay）层面


实验总数：**1,970,438**。


### 7.1 实验类型 assay_type

| 类型 | 实验数 |  | 含义 |
| --- | ---: | --- | --- |
| F | 921,809 | ██████████████████ | Functional：测功能后果（酶活、细胞信号、表型），靶点归属可能间接 |
| B | 615,457 | ████████████ | Binding：直接测化合物与靶点的结合（Ki/Kd/IC50），机制最明确，做 SAR/建模首选 |
| A | 325,436 | ██████ | ADME：吸收/分布/代谢/排泄性质（溶解度、渗透性、肝微粒体稳定性…） |
| T | 76,319 | █ | Toxicity：毒性相关 |
| P | 27,714 | █ | Physicochemical：纯理化性质，不涉及生物体系 |
| U | 3,703 |  | Unassigned：未分类 |


### 7.2 靶点指认可信度 confidence_score ⭐


**这是筛选数据时最重要的字段之一。** 它回答：「这个实验真的能归因到这个靶点吗？」

| 分数 | 实验数 |  | 含义 |
| ---: | ---: | --- | --- |
| 9 | 436,802 | ████████ | 同源蛋白复合物，靶点指认最明确 |
| 8 | 90,636 | ██ | 单一蛋白靶点，明确 — **建模常用的最低门槛** |
| 7 | 17,611 |  | 同源蛋白复合物（多亚基） |
| 6 | 2,846 |  | 蛋白复合物 |
| 5 | 30,162 | █ | 蛋白家族/未明确到具体成员 |
| 4 | 8,991 |  | 多个蛋白（如整条通路） |
| 3 | 20,956 |  | 细胞系或亚细胞组分层面 |
| 2 | 5,809 |  | 细胞系/亚细胞，靶点为推测 |
| 1 | 992,848 | ██████████████████ | 靶点仅为推测 |
| 0 | 363,777 | ███████ | 未指认靶点（如整体动物表型实验） |

> `confidence_score >= 8` 的实验共 **527,438**（26.8%）。做 QSAR / 机器学习建模时，这几乎是标配过滤条件。


### 7.3 实验体系物种（Top 20）


> 注意区分 `assay_organism`（实验体系来自哪个物种，例如用大鼠肝微粒体）与 `target_dictionary.organism`（靶点蛋白本身的物种）—— 两者可以不同。

| assay_organism | 实验数 |
| --- | ---: |
| Homo sapiens | 894,401 |
| Mus musculus | 269,412 |
| Rattus norvegicus | 206,361 |
| (空) | 167,981 |
| Staphylococcus aureus | 38,786 |
| Canis lupus familiaris | 30,333 |
| Escherichia coli | 25,908 |
| Pseudomonas aeruginosa | 15,721 |
| Human immunodeficiency virus 1 | 14,370 |
| Cavia porcellus | 14,190 |
| Candida albicans | 12,925 |
| Mycobacterium tuberculosis | 11,868 |
| Plasmodium falciparum | 10,717 |
| Oryctolagus cuniculus | 8,320 |
| Bos taurus | 7,999 |
| Macaca fascicularis | 7,737 |
| Macaca mulatta | 6,698 |
| Klebsiella pneumoniae | 6,295 |
| Salmonella enterica subsp. enterica serovar Typhimurium | 5,251 |
| Sus scrofa | 4,751 |


### 7.4 实验数最多的靶点（Top 20）

| ChEMBL ID | 靶点名 | 物种 | 实验数 |
| --- | --- | --- | ---: |
| `CHEMBL612558` | ADMET |  | 171,614 |
| `CHEMBL612545` | Unchecked |  | 160,540 |
| `CHEMBL375` | Mus musculus | Mus musculus | 113,921 |
| `CHEMBL376` | Rattus norvegicus | Rattus norvegicus | 106,633 |
| `CHEMBL352` | Staphylococcus aureus | Staphylococcus aureus | 35,254 |
| `CHEMBL3879801` | NON-PROTEIN TARGET |  | 34,765 |
| `CHEMBL2362975` | No relevant target |  | 31,623 |
| `CHEMBL387` | MCF7 | Homo sapiens | 19,772 |
| `CHEMBL354` | Escherichia coli | Escherichia coli | 19,608 |
| `CHEMBL392` | A549 | Homo sapiens | 19,101 |
| `CHEMBL373` | Canis familiaris | Canis lupus familiaris | 16,398 |
| `CHEMBL394` | HCT-116 | Homo sapiens | 14,777 |
| `CHEMBL400` | MDA-MB-231 | Homo sapiens | 14,208 |
| `CHEMBL348` | Pseudomonas aeruginosa | Pseudomonas aeruginosa | 14,060 |
| `CHEMBL399` | HeLa | Homo sapiens | 13,601 |
| `CHEMBL372` | Homo sapiens | Homo sapiens | 12,767 |
| `CHEMBL395` | HepG2 | Homo sapiens | 12,077 |
| `CHEMBL366` | Candida albicans | Candida albicans | 11,832 |
| `CHEMBL360` | Mycobacterium tuberculosis | Mycobacterium tuberculosis | 10,049 |
| `CHEMBL612546` | Molecular identity unknown |  | 9,868 |


## 8. 活性数据层面 ⭐


活性数据点总数：**24,527,044**。这是整个数据库的重心，也是最需要小心处理的部分。


### 8.1 数据完整性与质量标记

| 指标 | 数量 | 占比 | 备注 |
| --- | ---: | ---: | --- |
| 有标准化数值 `standard_value` | 21,129,848 | 86.1% | 没有数值的多为定性结论（Active/Inactive） |
| 有 `pchembl_value` | 4,969,278 | 20.3% | **做定量分析的可用子集** |
| `standard_relation = '='` | 15,047,664 | 61.4% | 精确值；其余为 > / < 的删失数据 |
| 被标记数据可疑 `data_validity_comment` | 353,204 | 1.4% | 建议剔除 |
| 疑似重复引用 `potential_duplicate = 1` | 355,508 | 1.4% | 同一数值被多篇文献转述，建议剔除 |
| 已人工标准化 `standard_flag = 1` | 15,880,656 | 64.7% |  |
| 有 `modality` 标注 | 29,083 | 0.1% | ChEMBL 37 起新增，目前主要标注靶向蛋白降解 |


### 8.2 测量指标 standard_type（Top 20）


每一行是一种「测的是什么量」。下面的中文解释是本脚本内置的，供非药学背景读者参考。

| standard_type | 数量 |  | 这是什么 |
| --- | ---: | --- | --- |
| `Potency` | 4,473,542 | ████████████████ | 效价，PubChem 高通量筛选里对 IC50/EC50 的统称 |
| `IC50` | 3,623,879 | █████████████ | 半数抑制浓度：让目标活性下降 50% 所需的化合物浓度。数值越小 = 越强 |
| `GI50` | 2,631,731 | █████████ | 抑制 50% 细胞生长所需浓度（细胞水平） |
| `Inhibition` | 1,624,513 | ██████ | 在某个固定浓度下的抑制百分比（%），不是浓度值 |
| `Activity` | 1,415,489 | █████ | 笼统的“活性”，单位/含义随实验而异，通常需要看 assay 描述 |
| `Percent Effect` | 1,328,366 | █████ | 在固定浓度下的效应百分比 |
| `Ki` | 887,151 | ███ | 抑制常数：抑制剂与靶点的结合亲和力（热力学量）。越小 = 结合越紧 |
| `k_off` | 826,525 | ███ | 解离速率常数：配体从靶点上脱落的快慢；k_off 越小停留时间越长 |
| `kon` | 826,356 | ███ | 结合速率常数：配体与靶点结合的快慢（结合动力学） |
| `MIC` | 798,852 | ███ | 最低抑菌浓度：抑制细菌生长的最低浓度（抗菌实验） |
| `EC50` | 613,608 | ██ | 半数效应浓度：产生 50% 最大效应所需浓度。数值越小 = 越强 |
| `INHIBITION` | 339,133 | █ |  |
| `AC50` | 286,628 | █ | 半数活性浓度（IC50/EC50 的中性叫法，常见于高通量筛选） |
| `Kd` | 213,575 | █ | 解离常数：配体-靶点结合亲和力。越小 = 结合越紧 |
| `Z score` | 147,592 | █ | 高通量筛选的标准化打分（相对对照组的偏离程度），不是浓度 |
| `Ratio IC50` | 147,260 | █ | 两个 IC50 的比值，常用于表示选择性 |
| `GI` | 132,461 |  | 生长抑制（百分比或定性） |
| `Tissue Severity Score` | 128,999 |  | 组织病理学评分（毒理实验中对病变严重程度的分级） |
| `Ratio` | 126,526 |  | 两个测量值的比值 |
| `CC50` | 107,782 |  | 细胞毒性浓度：杀死 50% 细胞所需浓度 |

> ⚠️ **陷阱：同一指标存在大小写/写法不同的变体**，它们在数据库里是不同的字符串，`GROUP BY standard_type` 会把它们拆开。本库中最主要的几组：
>
> - `Inhibition`（1,624,513）、`INHIBITION`（339,133）
> - `Activity`（1,415,489）、`activity`（373）
> - `kon`（826,356）、`Kon`（19）
> - `Kd`（213,575）、`KD`（3）
> - `T1/2`（95,426）、`t1/2`（3,565）
> - `CL`（82,022）、`Cl`（4,799）
> - `Residual Activity`（73,944）、`Residual activity`（2,326）
> - `Solubility`（62,193）、`solubility`（8,506）
>
> 分析前建议先 `UPPER(standard_type)` 归一，或明确列出你要的写法。


### 8.3 标准化单位 standard_units（Top 15）


> 标准化后，浓度一律为 **nM**。看到 `%` 说明是百分比读数（如抑制率），看到 `(空)` 多半是无量纲比值或定性结论。**不同单位的数值绝不能混在一起做统计。**

| 单位 | 数量 |
| --- | ---: |
| `nM` | 12,701,404 |
| `%` | 5,257,341 |
| `(空)` | 3,419,100 |
| `ug.mL-1` | 970,562 |
| `s-1` | 827,745 |
| `uM` | 249,613 |
| `hr` | 136,149 |
| `mm` | 96,101 |
| `ug ml-1` | 79,336 |
| `mg.kg-1` | 65,941 |
| `mL.min-1.kg-1` | 48,673 |
| `cells.uL-1` | 42,368 |
| `ng.hr.mL-1` | 40,941 |
| `degrees C` | 33,571 |
| `g` | 31,476 |


### 8.4 关系符 standard_relation

| 关系 | 数量 | 含义 |
| --- | ---: | --- |
| `=` | 15,047,664 | 精确值 |
| `(空)` | 7,458,386 |  |
| `>` | 1,603,288 | 大于（未达到效应，实际值更大 → 活性更弱） |
| `<` | 346,366 | 小于 |
| `<=` | 33,610 | 小于等于 |
| `>=` | 32,494 | 大于等于 |
| `~` | 4,777 | 约等于 |
| `>>` | 435 |  |
| `<<` | 24 |  |


### 8.5 pChEMBL 分布


**读法**：pChEMBL 9 = 1 nM（很强），7 = 100 nM（不错），5 = 10 µM（弱）。分布通常在 5–8 之间呈钟形，这既反映真实的活性分布，也反映发表偏倚（太弱的化合物没人报道）。

| pChEMBL 区间 | 对应浓度 | 数量 |  |
| --- | --- | ---: | --- |
| [0, 1) | — | 1 |  |
| [1, 2) | — | 29 |  |
| [2, 3) | 1–10 mM | 645 |  |
| [3, 4) | 0.1–1 mM | 5,926 |  |
| [4, 5) | 10–100 µM | 1,367,408 | ██████████████████████████████ |
| [5, 6) | 1–10 µM | 1,267,198 | ████████████████████████████ |
| [6, 7) | 0.1–1 µM | 947,217 | █████████████████████ |
| [7, 8) | 10–100 nM | 753,259 | █████████████████ |
| [8, 9) | 1–10 nM | 447,595 | ██████████ |
| [9, 10) | 0.1–1 nM | 150,578 | ███ |
| [10, 11) | 10–100 pM | 28,296 | █ |
| [11, 12) | ≤ 10 pM | 1,081 |  |
| [12, 13) | — | 18 |  |
| [13, 14) | — | 20 |  |
| [14, 15) | — | 7 |  |


### 8.6 被标记为可疑的数据

| data_validity_comment | 数量 |
| --- | ---: |
| Outside typical range | 345,635 |
| Potential transcription error | 4,398 |
| Potential missing data | 2,509 |
| Manually validated | 505 |
| Potential author error | 155 |
| Author confirmed error | 2 |


### 8.7 modality（作用模态）


ChEMBL 37 新增字段，标注化合物的设计模态。目前主要值是「靶向蛋白降解」（PROTAC、分子胶这类不是抑制靶点、而是诱导其被降解的分子）。注意它与 `action_type` 不同：标了 modality 的化合物**未必**有活性。

| modality | 数量 |
| --- | ---: |
| Targeted Protein Degradation | 29,083 |


### 8.8 活性数据点最多的靶点（Top 20）


**这张表本身就是一堂课。** 排在最前面的往往不是什么热门蛋白，而是 `Unchecked`、`ADMET`、细胞系（HepG2、MCF7）和整个物种（疟原虫、金黄色葡萄球菌）——因为大规模高通量筛选和表型筛选贡献了海量数据点，但它们的靶点归属是模糊的。**数据量大 ≠ 数据可用**。要找某个具体蛋白的可建模数据，必须叠加 `target_type = 'SINGLE PROTEIN'` 和 `confidence_score >= 8`。

| ChEMBL ID | 靶点名 | 物种 | 活性数 |
| --- | --- | --- | ---: |
| `CHEMBL612545` | Unchecked |  | 2,317,536 |
| `CHEMBL364` | Plasmodium falciparum | Plasmodium falciparum | 974,632 |
| `CHEMBL376` | Rattus norvegicus | Rattus norvegicus | 780,950 |
| `CHEMBL612558` | ADMET |  | 538,391 |
| `CHEMBL1075138` | Tyrosyl-DNA phosphodiesterase 1 | Homo sapiens | 345,639 |
| `CHEMBL375` | Mus musculus | Mus musculus | 299,086 |
| `CHEMBL3879801` | NON-PROTEIN TARGET |  | 264,560 |
| `CHEMBL352` | Staphylococcus aureus | Staphylococcus aureus | 241,275 |
| `CHEMBL3706568` | HEK-293T | Homo sapiens | 234,606 |
| `CHEMBL615023` | U2OS | Homo sapiens | 234,604 |
| `CHEMBL5314315` | Fibroblast | Homo sapiens | 229,743 |
| `CHEMBL2362975` | No relevant target |  | 220,860 |
| `CHEMBL360` | Mycobacterium tuberculosis | Mycobacterium tuberculosis | 209,493 |
| `CHEMBL395` | HepG2 | Homo sapiens | 204,643 |
| `CHEMBL612653` | Plasmodium berghei | Plasmodium berghei | 193,322 |
| `CHEMBL354` | Escherichia coli | Escherichia coli | 147,762 |
| `CHEMBL392` | A549 | Homo sapiens | 139,460 |
| `CHEMBL387` | MCF7 | Homo sapiens | 137,820 |
| `CHEMBL348` | Pseudomonas aeruginosa | Pseudomonas aeruginosa | 135,437 |
| `CHEMBL1293278` | Geminin | Homo sapiens | 128,009 |


### 8.9 加上质量过滤后的排行（Top 20）


过滤条件：`target_type = 'SINGLE PROTEIN'` + `confidence_score ≥ 8` + 有 pChEMBL + `standard_relation = '='` + 排除可疑与重复。

| ChEMBL ID | 靶点名 | 物种 | 可用活性数 |
| --- | --- | --- | ---: |
| `CHEMBL1293224` | Microtubule-associated protein tau | Homo sapiens | 95,340 |
| `CHEMBL1293231` | Nuclear receptor ROR-gamma | Mus musculus | 90,262 |
| `CHEMBL3577` | Aldehyde dehydrogenase 1A1 | Homo sapiens | 76,392 |
| `CHEMBL2026` | Beta-lactamase | Escherichia coli K-12 | 62,061 |
| `CHEMBL1293254` | Ferritin light chain | Equus caballus | 42,642 |
| `CHEMBL1293226` | Lysine-specific demethylase 4E | Homo sapiens | 40,061 |
| `CHEMBL1293235` | Prelamin-A/C | Homo sapiens | 36,678 |
| `CHEMBL2608` | Lysosomal alpha-glucosidase | Homo sapiens | 35,148 |
| `CHEMBL1293232` | Survival motor neuron protein | Homo sapiens | 33,735 |
| `CHEMBL1293303` | Nonstructural protein 1 | Influenza A virus | 32,314 |
| `CHEMBL3563` | Cruzipain | Trypanosoma cruzi | 31,453 |
| `CHEMBL340` | Cytochrome P450 3A4 | Homo sapiens | 29,092 |
| `CHEMBL1963` | Thyrotropin receptor | Homo sapiens | 28,877 |
| `CHEMBL6110` | Thioredoxin glutathione reductase | Schistosoma mansoni | 28,538 |
| `CHEMBL4096` | Cellular tumor antigen p53 | Homo sapiens | 26,222 |
| `CHEMBL1293255` | 15-hydroxyprostaglandin dehydrogenase [NAD(+)] | Homo sapiens | 24,850 |
| `CHEMBL1293248` | 4'-phosphopantetheinyl transferase ffp | Bacillus subtilis | 24,681 |
| `CHEMBL2392` | DNA polymerase beta | Homo sapiens | 23,156 |
| `CHEMBL1293294` | Ras-related protein Rab-9A | Homo sapiens | 21,980 |
| `CHEMBL4159` | 3-hydroxyacyl-CoA dehydrogenase type-2 | Homo sapiens | 20,628 |

**注意榜单可能还是不太对劲**：占据前列的（Microtubule-associated protein tau、Nuclear receptor ROR-gamma、Aldehyde dehydrogenase 1A1 ……）未必是最重要的药物靶点，往往只是因为 PubChem 上的大规模高通量筛选（一次几十万化合物）恰好打了这些靶点。质量标记只能保证「这条数据本身可信」，**保证不了「这批数据代表了该靶点的研究现状」**。


### 8.10 再叠加「仅科学文献来源」（`src_id = 1`，Top 20）


再限定到人工审编的文献数据后，榜首变成了 Epidermal growth factor receptor、D(2) dopamine receptor、Voltage-gated inwardly rectifying potassium channel KCNH2 等被药物化学界长期反复研究的靶点。与上一张表对比即可看出，**「哪个靶点数据最多」这个问题的答案，取决于你是否把高通量筛选数据算进来**。

| ChEMBL ID | 靶点名 | 物种 | 可用活性数 |
| --- | --- | --- | ---: |
| `CHEMBL203` | Epidermal growth factor receptor | Homo sapiens | 13,652 |
| `CHEMBL217` | D(2) dopamine receptor | Homo sapiens | 11,342 |
| `CHEMBL240` | Voltage-gated inwardly rectifying potassium channel KCNH2 | Homo sapiens | 11,273 |
| `CHEMBL247` | Human immunodeficiency virus type 1 reverse transcriptase | Human immunodeficiency virus 1 | 9,808 |
| `CHEMBL205` | Carbonic anhydrase 2 | Homo sapiens | 9,785 |
| `CHEMBL261` | Carbonic anhydrase 1 | Homo sapiens | 8,243 |
| `CHEMBL279` | Vascular endothelial growth factor receptor 2 | Homo sapiens | 8,105 |
| `CHEMBL253` | Cannabinoid receptor 2 | Homo sapiens | 7,509 |
| `CHEMBL233` | Mu-type opioid receptor | Homo sapiens | 7,265 |
| `CHEMBL220` | Acetylcholinesterase | Homo sapiens | 7,263 |
| `CHEMBL325` | Histone deacetylase 1 | Homo sapiens | 7,223 |
| `CHEMBL1163125` | Bromodomain-containing protein 4 | Homo sapiens | 7,007 |
| `CHEMBL214` | 5-hydroxytryptamine receptor 1A | Homo sapiens | 6,959 |
| `CHEMBL218` | Cannabinoid receptor 1 | Homo sapiens | 6,933 |
| `CHEMBL4822` | Beta-secretase 1 | Homo sapiens | 6,747 |
| `CHEMBL234` | D(3) dopamine receptor | Homo sapiens | 6,706 |
| `CHEMBL3594` | Carbonic anhydrase 9 | Homo sapiens | 6,605 |
| `CHEMBL4078` | Acetylcholinesterase | Electrophorus electricus | 6,603 |
| `CHEMBL224` | 5-hydroxytryptamine receptor 2A | Homo sapiens | 6,338 |
| `CHEMBL251` | Adenosine receptor A2a | Homo sapiens | 6,224 |

> **这三张表（8.8 → 8.9 → 8.10）是本报告最重要的一节。** 同一个数据库，只是换了过滤条件，「最热门靶点」的答案就完全不同。用 ChEMBL 时，**先想清楚你的科学问题决定了哪种过滤**，再写 SQL。


## 9. 文献与药物注释


- 文献/数据集条目：**101,100**；有 PubMed ID：88,134；有 DOI：93,787

- 文献类型：PUBLICATION 93,488；PATENT 7,145；DATASET 465；BOOK 2


### 9.1 主要期刊（Top 20）

| 期刊 | 文献数 |
| --- | ---: |
| J Med Chem | 26,839 |
| Bioorg Med Chem Lett | 24,365 |
| Eur J Med Chem | 11,799 |
| Bioorg Med Chem | 9,438 |
| J Nat Prod | 8,945 |
| ACS Med Chem Lett | 3,547 |
| Antimicrob Agents Chemother | 2,127 |
| Medchemcomm | 1,374 |
| Med Chem Res | 1,309 |
| RSC Med Chem | 692 |
| J Agric Food Chem | 422 |
| Drug Metab Dispos | 320 |
| J Pestic Sci | 245 |
| J Biol Chem | 185 |
| Nat Chem Biol | 177 |
| Crop Prot | 129 |
| Pest Manag Sci | 126 |
| Proc Natl Acad Sci U S A | 116 |
| J Pharmacol Exp Ther | 105 |
| Biosci Biotechnol Biochem | 64 |


### 9.2 文献年份分布（1990 起）

| 年份 | 文献数 |  |
| ---: | ---: | --- |
| 1990 | 556 | ███ |
| 1991 | 662 | ███ |
| 1992 | 958 | █████ |
| 1993 | 1,000 | █████ |
| 1994 | 1,099 | ██████ |
| 1995 | 1,170 | ██████ |
| 1996 | 1,245 | ██████ |
| 1997 | 1,160 | ██████ |
| 1998 | 1,372 | ███████ |
| 1999 | 1,411 | ███████ |
| 2000 | 1,428 | ███████ |
| 2001 | 1,446 | ███████ |
| 2002 | 1,672 | █████████ |
| 2003 | 1,790 | █████████ |
| 2004 | 2,117 | ███████████ |
| 2005 | 2,021 | ██████████ |
| 2006 | 2,145 | ███████████ |
| 2007 | 3,847 | ████████████████████ |
| 2008 | 4,087 | █████████████████████ |
| 2009 | 4,471 | ███████████████████████ |
| 2010 | 4,658 | ████████████████████████ |
| 2011 | 4,262 | ██████████████████████ |
| 2012 | 4,484 | ███████████████████████ |
| 2013 | 4,697 | ████████████████████████ |
| 2014 | 3,923 | ████████████████████ |
| 2015 | 3,893 | ████████████████████ |
| 2016 | 4,225 | ██████████████████████ |
| 2017 | 4,369 | ██████████████████████ |
| 2018 | 3,884 | ████████████████████ |
| 2019 | 3,899 | ████████████████████ |
| 2020 | 3,931 | ████████████████████ |
| 2021 | 3,862 | ████████████████████ |
| 2022 | 3,283 | █████████████████ |
| 2023 | 2,865 | ███████████████ |
| 2024 | 2,933 | ███████████████ |
| 2025 | 1,281 | ███████ |


### 9.3 药物层面的注释表


这几张表只覆盖**已知药物**（不是全部化合物），但信息密度很高，适合做药物重定位、靶点-疾病关联一类的分析。

| 表 | 行数 | 说明 |
| --- | ---: | --- |
| `drug_mechanism` | 7,561 | 药物 × 靶点 × 作用类型（抑制剂/激动剂…），人工审编 |
| `drug_indication` | 60,055 | 药物 × 适应症，疾病映射到 MeSH / EFO |
| `drug_warning` | 2,304 | 安全性警示（黑框警告、撤市） |
| `metabolism` | 2,147 | 代谢途径：底物 → 产物 → 催化酶 |
| `molecule_atc_classification` | 4,567 | 药物 × WHO ATC 分类 |
| `products` | 45,752 | FDA 批准的药品（制剂级） |
| `formulations` | 53,430 | 药品 → 成分化合物 |


### 9.4 药物作用类型 action_type（Top 15）


> `INHIBITOR`（抑制剂）降低靶点活性，`AGONIST`（激动剂）模拟天然配体激活靶点，`ANTAGONIST`（拮抗剂）阻断天然配体，`BLOCKER` 多用于离子通道。`action_type` 表还给出了上位归类（正向/负向调节）。

| action_type | 数量 |
| --- | ---: |
| `INHIBITOR` | 3,586 |
| `ANTAGONIST` | 980 |
| `AGONIST` | 949 |
| `(空)` | 577 |
| `BINDING AGENT` | 284 |
| `BLOCKER` | 179 |
| `MODULATOR` | 111 |
| `POSITIVE ALLOSTERIC MODULATOR` | 82 |
| `HYDROLYTIC ENZYME` | 77 |
| `ACTIVATOR` | 73 |
| `PARTIAL AGONIST` | 64 |
| `DISRUPTING AGENT` | 59 |
| `VACCINE ANTIGEN` | 56 |
| `EXOGENOUS PROTEIN` | 50 |
| `SEQUESTERING AGENT` | 46 |


## 10. 上手：几个可直接运行的查询


以下 SQL 可直接在本库上运行（`python3 -c` 或 `sqlite3` 均可）。


### 10.1 已知某个基因，取它的高质量活性数据


这是最常用的一条流水线：基因名 → 靶点 → 高可信实验 → 精确的定量活性。

```sql
SELECT md.chembl_id                AS compound,
       cs.canonical_smiles         AS smiles,
       act.standard_type,
       act.standard_value, act.standard_units,
       act.pchembl_value,
       td.pref_name                AS target,
       d.year, d.pubmed_id
FROM activities act
JOIN assays a               ON act.assay_id = a.assay_id
JOIN target_dictionary td   ON a.tid        = td.tid
JOIN target_components tc   ON tc.tid       = td.tid
JOIN component_sequences seq ON seq.component_id = tc.component_id
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON cs.molregno  = md.molregno
LEFT JOIN docs d            ON act.doc_id   = d.doc_id
WHERE seq.accession = 'P35557'          -- UniProt：人葡萄糖激酶 GCK
  AND td.target_type = 'SINGLE PROTEIN'
  AND a.confidence_score >= 8           -- 靶点指认可信
  AND act.pchembl_value IS NOT NULL     -- 只要能定量的
  AND act.standard_relation = '='       -- 排除删失数据
  AND act.data_validity_comment IS NULL -- 排除可疑数据
  AND act.potential_duplicate = 0       -- 排除重复引用
ORDER BY act.pchembl_value DESC;
```


### 10.2 同一化合物-靶点对有多条记录时怎么办


很常见：不同实验室、不同实验条件重复测过。标准做法是**取中位数**，并丢弃离散度过大的对（例如 max−min > 1 个 log 单位说明数据打架）。

```sql
SELECT md.chembl_id, COUNT(*) n,
       ROUND(AVG(act.pchembl_value), 2)  AS mean_p,
       ROUND(MAX(act.pchembl_value) - MIN(act.pchembl_value), 2) AS spread
FROM activities act
JOIN assays a ON act.assay_id = a.assay_id
JOIN molecule_dictionary md ON act.molregno = md.molregno
WHERE a.tid = ? AND act.pchembl_value IS NOT NULL
GROUP BY md.chembl_id
HAVING n > 1
ORDER BY spread DESC;
```


### 10.3 查一个陌生的 CHEMBL ID 是什么


```sql
SELECT entity_type, entity_id FROM chembl_id_lookup WHERE chembl_id = 'CHEMBL25';
```


### 10.4 Python + pandas 的读法


```python
import sqlite3, pandas as pd
con = sqlite3.connect('file:/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db?mode=ro', uri=True)   # 只读，避免误写
df = pd.read_sql_query(open('query.sql').read(), con)
```

> 30 GB 的库不要 `SELECT *` 全表读进内存。先用 SQL 过滤好，再交给 pandas。


## 11. 使用这份数据前必须知道的坑


| 坑 | 后果 | 对策 |
| --- | --- | --- |
| 用了 `value` 而不是 `standard_value` | 单位混杂（µM 与 nM 混在一起），结论全错 | 永远用 `standard_*` 列 |
| 忽略 `standard_relation` | 把 `>10 µM` 当成 10 µM，把「无活性」当成「弱活性」 | 定量建模时限定 `= '='`；或按删失数据处理 |
| 不过滤 `confidence_score` | 把细胞表型、整体动物数据当成直接的靶点活性 | 建模用 `>= 8` |
| 混用不同 `standard_type` | IC50 与 Ki 与 %抑制率不可比 | 分开处理；至少分开 IC50/EC50 与 Ki/Kd |
| 忽略 `potential_duplicate` / `data_validity_comment` | 同一数值被重复计数；纳入已知错误值 | 两者都加进过滤条件 |
| 把化合物记录数当成化合物数 | `compound_records` 是文献级的，数量远大于唯一化合物 | 唯一化合物看 `molecule_dictionary` |
| 把盐和母体当成两个分子 | 同一药物的不同盐型被算作不同化合物 | 用 `molecule_hierarchy` 归并到 parent_molregno |
| 假设「没有数据 = 没有活性」 | ChEMBL 存在强烈发表偏倚，阴性结果严重缺失 | 做机器学习时需谨慎构造负样本 |
| 跨版本直接比较数字 | 官方会回溯性地修订（合并重复、重新指认靶点、改 max_phase） | 论文里注明具体版本号 |


## 12. 延伸资料


- 官方 schema 文档：https://chembl.gitbook.io/chembl-interface-documentation/db-schema-description
- 同目录下的 `schema_documentation.txt`（逐字段说明）与 `chembl_*_schema.pdf`（ER 图）
- 官方 release notes：`chembl_*_release_notes.txt`
- Web 界面（适合抽查单个化合物/靶点）：https://www.ebi.ac.uk/chembl/
- 数据许可 CC BY-SA；发表时需引用 ChEMBL 论文并注明 release 号（见 `REQUIRED.ATTRIBUTION`）

---

*本报告由 `chembl_profile.py` 自动生成。表中的中文解释为脚本内置的领域注解，统计数字均实时查询自上述数据库文件。*
