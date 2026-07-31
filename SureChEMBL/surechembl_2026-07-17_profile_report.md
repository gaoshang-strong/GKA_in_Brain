# SureChEMBL 专利化学数据库结构与内容概览报告

> 自动生成于 2026-07-31 16:53　|　快照目录：`/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17`
>
> 本报告面向**有生物学 / 生信背景、但没有专利检索经验**的读者。每一节先解释「这是什么」，再给出这份快照里的实际统计。


## 0. 三十秒速览

**SureChEMBL 是什么？** EMBL-EBI 维护的**专利化学数据库**。它把专利全文（含扫描件 OCR、结构图识别）里的化学结构和生物医学名词自动抽取出来，变成可查询的表。与 ChEMBL 的关键区别：

|  | ChEMBL | SureChEMBL |
| --- | --- | --- |
| 数据来源 | 论文 + 部分专利 | **专利全文** |
| 加工方式 | **人工审编** | **全自动抽取**（NER + OSR + OCR） |
| 核心事实 | 化合物 × 靶点 × 活性数值 | 化合物 × 专利 × **出现位置** |
| 有没有活性数据 | 有（IC50/EC50…） | **没有**。只知道「这个化合物出现在这篇专利里」 |
| 规模 | ~292 万化合物 | ~3,099 万化合物 |

**最重要的一条**：SureChEMBL **不告诉你化合物有没有活性**，只告诉你它出现在哪篇专利的哪个部分。活性要回 ChEMBL 或读专利原文。

**一句话的数据模型**：
```
某个化合物  出现在某篇专利  的某个部分
(compounds)   (patents)       (fields)
         └────────┬────────────────┘
            patent_compound_map

某个生物医学实体  出现在某篇专利  的某个部分  出现了几次
(biomedical_entities) (patents)   (fields)    (count)
         └────────────┬─────────────────────────┘
                biomedical_locations
```

| 条目 | 值 |
| --- | --- |
| 快照目录 | `SureChEMBL_2026-07-17` |
| 文件总大小 | 16.94 GB |
| 专利数（`patents`） | 44,912,542 |
| 化合物数（`compounds`） | 30,990,818 |
| 专利-化合物关联数（`patent_compound_map`） | 1,537,106,020 |
| 实体位置记录数（`biomedical_locations`） | 453,327,904 |
| 生物医学实体数（`biomedical_entities`） | 1,059,724 |
| duckdb 版本 | 1.5.5 |

> ⚠ **这份快照是某一期的全量镜像，不是增量。** SureChEMBL 每两周发一版并覆盖 `latest/`，因此本项目固定在日期目录上。任何数字都对应这一版。


## 1. 先搞懂这些词


### 1.1 专利文档（patent document）≠ 发明

同一项发明会在多个国家/多个阶段分别公开，每次公开都是**一篇独立的专利文档**，有各自的专利号。所以「专利数」远大于「发明数」。

把同一发明的所有文档串起来的字段是 **`family_id`（专利同族）**。**做去重统计必须按 `family_id`，不能按 `patent_number`**，否则一个化合物会因为在 5 个国家申请而被数 5 次。

### 1.2 `WO` 不是一个国家

`country` 里的 `WO` 是 WIPO 的 **PCT 国际申请**——申请人先交一份国际申请，之后再决定进入哪些国家。它和 `US`/`EP`/`JP`/`CN` 不是并列的地理关系，**同一发明常常既有 WO 文档又有各国文档**。

### 1.3 专利被切成 6 个部分（`fields`）⭐

这是全库**最容易踩坑**的设计。一篇专利被切成标题/摘要/说明书/权利要求/图片/MOL 附件，`patent_compound_map.field_id` 指出化合物是在哪个部分被发现的。

| id | field_name | 说明 |
| --- | --- | --- |
| 1 | `desc` | 说明书正文。**包含背景技术**，所以这里出现的化合物很多是他人的、被引用的 |
| 2 | `clms` | **权利要求**。这是专利真正主张保护的范围——判断「这篇专利要保护什么化合物」只能看这里 |
| 3 | `abst` | 摘要 |
| 4 | `ttl` | 标题。信息密度最高但覆盖最少 |
| 5 | `image` | 从图片里识别出的结构（化学结构图 OSR） |
| 6 | `molattachment` | 专利附带的 MOL 文件（序列表/结构文件） |

**`desc`（说明书）里出现 ≠ 这篇专利要保护它。** 说明书含背景技术，会大段引用他人的化合物做对比。真正主张保护的范围只在 **`clms`（权利要求）**。

### 1.4 实体解析（entity resolution）：`resolved_form` 才是锚点

`biomedical_entities` 里，`original_text` 是原文里出现的字面写法，`resolved_form` 是归一后的标识符。因为专利全文来自扫描件 OCR，**同一个概念会有几十种字面写法**（断字、错字、缩写）。

**按 `original_text` 做字符串匹配一定会漏，必须按 `resolved_form` 锚定。**

| id | type_name | 说明 | 本快照实体数 |
| ---: | --- | --- | ---: |
| 1 | `GeneOrProtein` | 基因/蛋白名，归一到 UniProt、HGNC、Entrez Gene | 913,773 |
| 2 | `Disease` | 疾病名，归一到 MeSH、Disease Ontology、Wikipedia | 137,749 |
| 3 | `Mechanism` | 化学物质的作用或作用机制 | 8,202 |
| 4 | `Physquant` | 各类物理量 | 0 |

**只有 497,974 / 1,059,724（47.0%）的实体有 `resolved_form`**，其余是抽出来了但没能归一到任何标识符的词。未归一的实体**无法可靠地用于检索**——你不知道它到底指什么。

⚠ `resolved_form` **不是单一命名空间**：同一列里混着 `HGNC:4195`（带前缀）、`Q14397`（裸 UniProt accession）、MeSH ID 等。解析前必须先判断命名空间。

### 1.5 分类号：CPC / IPC / IPCR / ECLA

| 字段 | 说明 |
| --- | --- |
| `cpc` | Cooperative Patent Classification，EPO 与 USPTO 联合分类体系，最细（约 25 万个条目） |
| `ipc` | International Patent Classification，WIPO 的国际分类，粗一些 |
| `ipcr` | IPC 的 reformed 版本（2006 年后的写法） |
| `ecla` | European Classification，EPO 的旧体系，**已被 CPC 取代**，新专利基本没有 |

这四个字段在 `patents.parquet` 里都是 **list 类型**（一篇专利可以有多个分类号），**`LIKE` 匹配不上**，要用 `list_contains()` / `UNNEST()`。

## 2. 数据模型与表清单

Schema 镜像 SureChEMBL 内部关系库，**不做冗余展开**——这是 2.0 相对旧版 MAP files 的主要改进，文件体积因此小很多。

```
compounds.id ──────┐
                   ├──< patent_compound_map >── fields.id
patents.id ────────┘                    (field_id)

biomedical_entities.id ──┐
                         ├──< biomedical_locations >── fields.id
patents.id ──────────────┘                    (field_id)
        │
        └── biomedical_entities.type_id ──> biomedical_types.id
```

**没有外键约束**（parquet 不支持），完整性靠上游保证——本报告 §7 有实测的完整性检查。

| 表 | 行数 | row group | 文件大小 | 说明 |
| --- | ---: | ---: | ---: | --- |
| `patents` | 44,912,542 | 1,000 | 5.51 GB | 专利主表：一行一篇专利文档，含专利号、国别、公开日、同族、分类号、申请人、标题 |
| `compounds` | 30,990,818 | 2,385 | 3.91 GB | 化合物主表：一行一个唯一化学结构（SMILES / InChI / InChIKey / 分子量） |
| `patent_compound_map` | 1,537,106,020 | 5,006 | 4.64 GB | 核心事实表：某化合物出现在某专利的某个部分。**全库最大，15.4 亿行** |
| `fields` | 6 | 1 | 2 KB | 字典：专利被切成哪几个部分（标题/摘要/说明书/权利要求/图片/MOL 附件） |
| `biomedical_entities` | 1,059,724 | 2 | 31.66 MB | 生物医学实体词表：文中出现的基因/蛋白/疾病/机制词，及其归一化 ID |
| `biomedical_types` | 4 | 1 | 3 KB | 字典：实体的四种类型 |
| `biomedical_locations` | 453,327,904 | 433 | 1.60 GB | 实体位置表：某实体出现在某专利某部分，及出现次数 |
| `fpsim2_fingerprints.h5` | — | — | 1.26 GB | FPSim2 指纹库（HDF5，非 parquet），可直接做全库相似性检索 |


### 2.1 各表字段（从 parquet schema 实读）

- **`patents`** — `id`:INT64, `patent_number`:BYTE_ARRAY, `country`:BYTE_ARRAY, `publication_date`:INT32, `family_id`:INT64, `element`:BYTE_ARRAY, `element`:BYTE_ARRAY, `element`:BYTE_ARRAY, `element`:BYTE_ARRAY, `element`:BYTE_ARRAY, `title`:BYTE_ARRAY
- **`compounds`** — `id`:INT64, `smiles`:BYTE_ARRAY, `inchi`:BYTE_ARRAY, `inchi_key`:BYTE_ARRAY, `mol_weight`:DOUBLE
- **`patent_compound_map`** — `patent_id`:INT64, `compound_id`:INT64, `field_id`:INT64
- **`fields`** — `id`:INT64, `field_name`:BYTE_ARRAY
- **`biomedical_entities`** — `id`:INT64, `type_id`:INT64, `corrected_text`:BYTE_ARRAY, `original_text`:BYTE_ARRAY, `resolved_form`:BYTE_ARRAY
- **`biomedical_types`** — `id`:INT64, `type_name`:BYTE_ARRAY, `description`:BYTE_ARRAY
- **`biomedical_locations`** — `entity_id`:INT64, `patent_id`:INT64, `field_id`:INT64, `count`:INT64

注意 `patents` 的 `cpc` / `ipcr` / `ipc` / `ecla` / `assignee` 是 **LIST 类型**，`compounds` 与 `patents` 的字符串列是 `large_string`（BYTE_ARRAY）。

## 3. 专利层面


### 3.1 专利局分布

| country | 专利数 | 占比 |  | 说明 |
| --- | ---: | ---: | --- | --- |
| `CN` | 23,884,165 | 53.2% | ██████████████████████████ | 中国国家知识产权局（CNIPA）。SureChEMBL 2.0 新增，专利数最多但化学信息密度最低 |
| `US` | 9,691,977 | 21.6% | ███████████ | 美国专利商标局（USPTO）。化合物贡献最大 |
| `EP` | 5,265,735 | 11.7% | ██████ | 欧洲专利局（EPO） |
| `JP` | 3,062,582 | 6.8% | ███ | 日本特许厅（JPO） |
| `WO` | 3,008,081 | 6.7% | ███ | 世界知识产权组织（WIPO）的 PCT 国际申请。**不是某国专利**，是进入各国前的国际阶段申请 |
| `GB` | 2 | 0.0% | █ | **异常值，数量极少，来源不明** |

⚠ 出现了不在官方 5 家专利局之列的 `country` 值：GB。数量极少，做统计时应显式排除或单独核查。


### 3.2 公开年份分布

| 年份 | 专利数 |  |
| ---: | ---: | --- |
| 2026 | 1,006,583 | ███████ |
| 2025 | 2,358,189 | █████████████████ |
| 2024 | 3,543,402 | █████████████████████████ |
| 2023 | 3,543,471 | █████████████████████████ |
| 2022 | 3,640,429 | ██████████████████████████ |
| 2021 | 2,650,719 | ███████████████████ |
| 2020 | 2,325,791 | █████████████████ |
| 2019 | 2,084,398 | ███████████████ |
| 2018 | 2,073,236 | ███████████████ |
| 2017 | 1,788,517 | █████████████ |
| 2016 | 1,659,983 | ████████████ |
| 2015 | 1,507,571 | ███████████ |
| 2014 | 1,313,507 | █████████ |
| 2013 | 1,253,427 | █████████ |
| 2012 | 1,110,252 | ████████ |
| 2011 | 940,445 | ███████ |
| 2010 | 879,852 | ██████ |
| 2009 | 807,432 | ██████ |
| 2008 | 762,499 | █████ |
| 2007 | 711,521 | █████ |
| 2006 | 677,805 | █████ |
| 2005 | 633,920 | █████ |
| 2004 | 599,823 | ████ |
| 2003 | 559,521 | ████ |
| 2002 | 504,263 | ████ |

`publication_date` 为空的有 **1,172,063** 篇（2.6%）。最新一年通常不完整（快照日期之后的还没收录）。


### 3.3 专利同族 family_id ⭐

**43,668,617 篇有真实同族的专利文档只对应 26,481,620 个同族**（平均每族 1.65 篇）。

**这是去重的关键**：直接数专利篇数会把同一发明重复计入。做「有多少个 GKA 发明」这类统计，一律 `COUNT(DISTINCT family_id)`。

⚠ **但 `family_id` 有两类无效值，必须先排除，否则 `COUNT(DISTINCT)` 会被污染**：

| 值 | 文档数 | 含义 |
| --- | ---: | --- |
| `-1`（哨兵） | 71,862 | **不是同族编号**，是「未分配同族」的占位。不排除的话这 71,862 篇会被当成同一个发明 |
| `NULL` | 1,172,063 | 字段缺失 |

```sql
-- 正确写法
COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)
```

另外实测：**`family_id` 为空的 1,172,063 篇，与 `publication_date` 为空的完全是同一批**。说明这是一整块元数据缺失的记录，不是随机缺失——做时间趋势或同族分析时它们会整体消失，要意识到这个盲区。

最大的几个真实同族（一项发明在全球公开了多少次）：

| family_id | 文档数 |
| --- | ---: |
| 21841414 | 616 |
| 27290628 | 398 |
| 25646068 | 356 |
| 21978554 | 339 |
| 34749937 | 323 |


### 3.4 申请人 assignee

`assignee` 是 LIST 类型，要 `UNNEST` 后再统计。

| 申请人 | 专利数 |
| --- | ---: |
| SAMSUNG ELECTRONICS CO., LTD. (KR) | 205,122 |
| IBM (US) | 186,000 |
| SAMSUNG ELECTRONICS CO LTD (KR) | 182,822 |
| INTERNATIONAL BUSINESS MACHINES CORPORATION | 165,996 |
| CANON KABUSHIKI KAISHA (JP) | 121,949 |
| MICRON TECHNOLOGY, INC. | 109,826 |
| INTERNATIONAL BUSINESS MACHINES CORPORATION (US) | 108,942 |
| CANON KK (JP) | 102,994 |
| BASF SE (DE) | 101,359 |
| SIEMENS AG (DE) | 100,343 |
| KABUSHIKI KAISHA TOSHIBA (JP) | 96,601 |
| MORGAN STANLEY SENIOR FUNDING, INC. | 92,307 |
| BOSCH GMBH ROBERT (DE) | 87,269 |
| TOYOTA JIDOSHA KABUSHIKI KAISHA (JP) | 85,517 |
| 中国石油化工股份有限公司 | 84,437 |

申请人名称**没有做机构消歧**——同一家公司会有多种写法（含子公司、不同语言转写、OCR 噪声）。要按公司统计必须自己归一。


## 4. 化合物层面

**30,990,818** 个唯一化学结构。每行给 `smiles` / `inchi` / `inchi_key` / `mol_weight`，**没有活性、没有名称、没有理化性质**——要这些得回 ChEMBL 或自己算。


### 4.1 结构标识的完整性

| 字段 | 有值 | 缺失 | 缺失率 |
| --- | ---: | ---: | ---: |
| `smiles` | 30,990,818 | 0 | 0.0% |
| `inchi` | 30,990,818 | 0 | 0.0% |
| `inchi_key` | 30,990,818 | 0 | 0.0% |

`inchi_key` 唯一值 **29,874,136** / 30,990,818——有 1,116,682 个重复，跨库对齐前要先去重。


### 4.2 分子量分布

| 统计量 | 值 |
| --- | ---: |
| 有值 | 30,990,818 |
| 最小 | 0 |
| 25% | 288.23 |
| 中位 | 402.43 |
| 75% | 516.60 |
| 99% | 1,193.43 |
| 最大 | 141,073 |

⚠ **两端都有脏数据**：`mol_weight = 0` 的有 **2,470** 条，`> 2000`（超出小分子范围，多为聚合物/多肽/OSR 误识别）有 **46,963** 条，最大值 141,073。**做小分子筛选务必加分子量下限与上限**。

> 抽取是全自动的（含化学结构图 OSR 与 OCR），**必然存在错误结构**。SureChEMBL 的化合物没有经过人工审编，这与 ChEMBL 是本质区别——拿来做训练集或候选池前要自己过一遍 RDKit 合法性与合理性检查。


## 5. 专利 ↔ 化合物关联（核心表）⭐

**1,537,106,020 行**，全库最大。三列：`patent_id`、`compound_id`、`field_id`。


### 5.1 按出现部分 field_id 拆开 ⭐

| field_id | field_name | 关联数 | 占比 |  |
| ---: | --- | ---: | ---: | --- |
| 1 | `desc` | 1,218,172,303 | 79.3% | ██████████████████████████ |
| 2 | `clms` | 186,683,469 | 12.1% | ████ |
| 3 | `abst` | 52,098,610 | 3.4% | █ |
| 5 | `image` | 46,057,688 | 3.0% | █ |
| 6 | `molattachment` | 24,041,451 | 1.6% | █ |
| 4 | `ttl` | 10,052,499 | 0.7% | █ |

**说明书（`desc`）占了 79.3%，权利要求（`clms`）只有 12.1%，相差 6.5 倍。**

这就是本库最大的坑：**如果不按 `field_id` 过滤，你拿到的绝大部分是「说明书里被提到过」的化合物**——包括背景技术里引用的他人化合物、对比例、乃至试剂。判断「这篇专利要保护什么」必须 `field_id = 2`。


### 5.2 参与关联的实体数

|  | 数量 | 占该表总数 |
| --- | ---: | ---: |
| 出现在关联表里的专利 | 44,912,568 | 100.0% |
| 出现在关联表里的化合物 | 30,990,818 | 100.0% |


### 5.3 每篇专利含多少化合物

| 最小 | 中位 | 90% | 99% | 最大 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3 | 59 | 452 | 24,225 |

**分布极度右偏**：中位只有 3 个，但最大一篇有 **24,225** 个化合物——那是马库什结构（Markush）被枚举展开的组合库专利。**这类专利会主导任何按化合物计数的统计**，做分析时要么按专利归一，要么把超大专利单列。


## 6. 生物医学标注

`biomedical_entities` **1,059,724** 个实体，`biomedical_locations` **453,327,904** 条位置记录。


### 6.1 各类型实体的出现规模

| 类型 | 实体数 | 说明 |
| --- | ---: | --- |
| `GeneOrProtein` | 913,773 | 基因/蛋白名，归一到 UniProt、HGNC、Entrez Gene |
| `Disease` | 137,749 | 疾病名，归一到 MeSH、Disease Ontology、Wikipedia |
| `Mechanism` | 8,202 | 化学物质的作用或作用机制 |

⚠ 字典里定义了但**本快照一条实体都没有**的类型：Physquant。字典有定义不等于有数据。


### 6.2 归一化率

| 类型 | 实体数 | 有 resolved_form | 归一化率 |
| --- | ---: | ---: | ---: |
| `GeneOrProtein` | 913,773 | 361,589 | 39.6% |
| `Disease` | 137,749 | 136,385 | 99.0% |
| `Mechanism` | 8,202 | 0 | 0.0% |

**归一化率低意味着大量实体只能当自由文本用。** 检索时按 `resolved_form` 锚定虽然会漏掉未归一的那部分，但按 `original_text` 匹配又会引入大量歧义——两头都有代价，正确做法是**先按 `resolved_form` 取，再人工看未归一实体里有没有该捞的**。


### 6.3 GCK / glucokinase 的实测（本项目锚点）⭐

`resolved_form = 'HGNC:4195'`（人 GCK 基因）对应 **33 个不同的 `original_text`**：

```
glucokinase、GK、Hexokinase 4、gluco- kinase、GlcK、gluco kinase、giucokinase、gluco-kinase、glucokmase、Glucoki- nase、glk、GlkA、glu- cokinase、Glucoki nase、GKAs、glucokina se、GKA、GIu- cokinase、glu- co kinase、gl uc ok i na se、Glu cokinase、Glu-cokinase、Hexokinase-4、glucok-inase、4、Hk4、gki、Glucokin ase、Glucokinas e、GcK、gukA、GluK、glucok inase
```

里面有 OCR 破碎形（`glucokmase`、`gl uc ok i na se`）、ChEMBL 那边的名字（`Hexokinase 4`）、以及一批缩写。**按字符串搜 `glucokinase` 会漏掉其中大部分。**

| field | 专利数 | 提及次数 |
| --- | ---: | ---: |
| `desc` | 28,272 | 123,965 |
| `clms` | 2,597 | 10,442 |
| `abst` | 1,462 | 3,189 |
| `ttl` | 1,085 | 1,925 |

各写法贡献的专利数：

| original_text | 专利数 | 备注 |
| --- | ---: | --- |
| `glucokinase` | 29,191 |  |
| `GK` | 1,343 | ⚠ 糖尿病文献里 `GK` 更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型），与本领域高度混淆 |
| `GKAs` | 131 |  |
| `Hexokinase 4` | 103 |  |
| `glucokmase` | 80 |  |
| `giucokinase` | 51 |  |
| `GlcK` | 48 | 细菌 glucokinase 基因名，不是人 GCK |
| `glu- cokinase` | 37 |  |
| `Hexokinase-4` | 25 |  |
| `gluco- kinase` | 20 |  |
| `GKA` | 13 |  |
| `GcK` | 11 |  |

⚠ 另外这些**不**解析到 HGNC:4195，别混进来：`Glucokinase regulator` → `HGNC:4196`（GCKR/GKRP，不是 GCK）、全大写 `GCK` → **空（未解析）**、`HPK/GCK-like kinase` → `HGNC:6866`（MAP4K 家族）。


## 7. 数据完整性实测

| 检查项 | 数量 | 结论 |
| --- | ---: | --- |
| map 里引用了但 `patents` 中不存在的 `patent_id` | 26 | ⚠ 外键悬空 |
| map 里引用了但 `compounds` 中不存在的 `compound_id` | 0 | ✅ |
| `patents` 里一个化合物都没有的专利 | 0 | 见下方说明 |

**「一个化合物都没有的专利 = 0」是个重要事实**：说明 bulk 导出**只收录了含化学结构的专利**，不是 SureChEMBL 系统里的全部文档。这解释了为什么本地实测的专利数（约 4,491 万）远小于官方宣传的 1.166 亿——两者口径不同。**做规模陈述以本地实测为准。**


## 8. 上手：可直接运行的查询

环境：micromamba `GKA_in_Brain`（`duckdb` + `pyarrow`）。**所有查询都直接打 parquet，不需要先导入数据库。**


### 8.1 基本设置

```python
import duckdb
con = duckdb.connect()
con.execute("SET threads=8; SET memory_limit='12GB';")
con.execute("SET enable_progress_bar=false;")   # 不关会把 stdout 刷爆

D = "/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17"
```


### 8.2 找出在权利要求里主张 glucokinase 的专利

注意三点：按 `resolved_form` 锚定、按 `field_id=2` 限定权利要求、按 `family_id` 去重。
```python
con.sql(f'''
WITH e AS (                                   -- 1. 按归一 ID 锚定，不用字符串匹配
  SELECT id FROM '{D}/biomedical_entities.parquet'
  WHERE resolved_form = 'HGNC:4195'
),
p AS (                                        -- 2. 只要权利要求里提到的
  SELECT DISTINCT l.patent_id
  FROM '{D}/biomedical_locations.parquet' l
  JOIN e ON e.id = l.entity_id
  WHERE l.field_id = 2
)
SELECT COUNT(*) AS n_docs,                    -- 3. 同族去重后才是「发明数」
       COUNT(DISTINCT pt.family_id) AS n_families
FROM p JOIN '{D}/patents.parquet' pt ON pt.id = p.patent_id
''').show()
```


### 8.3 取某篇专利在权利要求里保护的化合物结构

```python
con.sql(f'''
SELECT c.id, c.smiles, c.inchi_key, c.mol_weight
FROM '{D}/patent_compound_map.parquet' m
JOIN '{D}/compounds.parquet' c ON c.id = m.compound_id
WHERE m.patent_id = 12345678
  AND m.field_id = 2                          -- 只要权利要求
''').show()
```


### 8.4 反查：某个结构出现在哪些专利里

先用 InChIKey 定位 compound_id，再查关联表——**别拿 SMILES 做字符串匹配**，同一结构有多种 SMILES 写法。
```python
con.sql(f'''
WITH c AS (
  SELECT id FROM '{D}/compounds.parquet'
  WHERE inchi_key = 'KJSGTWFWVTYPFZ-AWEZNQCLSA-N'   -- MK-0941
)
SELECT pt.patent_number, pt.country, pt.publication_date, m.field_id
FROM '{D}/patent_compound_map.parquet' m
JOIN c ON c.id = m.compound_id
JOIN '{D}/patents.parquet' pt ON pt.id = m.patent_id
ORDER BY pt.publication_date
''').show()
```


### 8.5 按分类号过滤（LIST 类型，不能用 LIKE）

```python
con.sql(f'''
SELECT COUNT(*) FROM '{D}/patents.parquet'
WHERE list_contains(cpc, 'A61P3/10')          -- 抗糖尿病用途
''').show()
```


### 8.6 相似性检索用 FPSim2，不要自己算指纹

```python
# 官方已提供 3,099 万化合物的预计算指纹（1.26 GB）
from FPSim2 import FPSim2Engine
fpe = FPSim2Engine("fpsim2_fingerprints.h5")
results = fpe.similarity("CCS(=O)(=O)c1ccc(...)cn1", 0.7, n_workers=4)
# 返回 (compound_id, similarity)，再拿 compound_id 回 patent_compound_map
```

> FPSim2 需另装（`pip install FPSim2`），本项目环境暂未安装。


## 9. 用这份数据前必须知道的坑

按踩坑代价从大到小排。

1. **`field_id` 不过滤，结论就是错的**
   说明书关联数是权利要求的 6.5 倍。「专利里出现过某化合物」和「这篇专利要保护它」完全是两回事——说明书含背景技术，会大段引用他人化合物。**判断专利主张范围只能用 `field_id = 2`。**

2. **不按 `family_id` 去重，数字会虚高；去重时不排除哨兵值，又会算少**
   同一项发明在多国多阶段公开，每次都是一篇独立文档，「有多少个 GKA 发明」必须 `COUNT(DISTINCT family_id)`。但 `family_id = -1` 是「未分配同族」的**哨兵值**（7 万余篇），直接 `COUNT(DISTINCT)` 会把它们错当成同一个发明。正确写法：`COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)`。

3. **按 `original_text` 做文本匹配一定会漏**
   专利全文来自扫描件 OCR，同一概念有几十种字面写法。必须按 `resolved_form` 锚定。但反过来，缩写形（如 `GK`）会引入大量歧义，**锚定后仍需逐个 surface form 评估**。

4. **这里没有活性数据**
   SureChEMBL 只说「化合物出现在专利里」，不说它有没有效、多强。活性要回 ChEMBL 或读专利原文。**不能把「出现在 GCK 专利里」当成「是 GKA」。**

5. **结构是全自动抽取的，没有人工审编**
   含化学结构图识别（OSR）与 OCR，必然有错误结构。`mol_weight = 0` 与 `> 2000` 的都存在。拿来做候选池前要自己过 RDKit 合法性检查并设分子量窗口。

6. **马库什结构会主导计数**
   组合库专利可以枚举出上万个化合物，任何「按化合物计数」的统计都会被这类专利带偏。要么按专利归一，要么单列。

7. **`resolved_form` 不是单一命名空间**
   `HGNC:4195` 带前缀，`Q14397` 是裸 UniProt accession。解析前先判断命名空间。

8. **`cpc` / `ipc` / `assignee` 是 LIST 类型**
   `LIKE` 匹配不上，要用 `list_contains()` / `UNNEST()`。申请人名称未做机构消歧，同一公司多种写法。

9. **`WO` 不是国家**
   是 PCT 国际申请，与各国文档并存，不能和 `US`/`CN` 并列做「国别分布」解读。

10. **规模数字以本地实测为准**
   官方宣传 1.166 亿专利 / 4,770 万化合物，本地 bulk 快照实测是 4,491 万 / 3,099 万——bulk 只含有化学标注的专利，口径不同。

11. **快照两周一覆盖**
   `latest/` 会被覆盖，必须固定到日期目录，否则结论无法复现。


## 10. 与 ChEMBL 侧的衔接

本项目的 ChEMBL 侧锚点是 **UniProt `P35557`**，SureChEMBL 侧是 **`HGNC:4195`**，两者是同一个基因（GCK）的不同命名空间。**跨库对齐走基因层，不要字符串比对蛋白名。**

化合物层面的对齐用 **InChIKey**——两边都有，且都是标准 InChI 生成的。注意 InChIKey 前 14 位相同只代表骨架相同，**立体异构体会被合并**，而手性通常决定 GKA 的活性。

|  | ChEMBL | SureChEMBL |
| --- | --- | --- |
| 靶点锚点 | `P35557`（UniProt） | `HGNC:4195` |
| 化合物对齐键 | `standard_inchi_key` | `inchi_key` |
| 能回答 | 这个化合物**有多强** | 这个化合物**在谁的专利里** |
| 数据质量 | 人工审编 | 全自动抽取，需自己过滤 |


---

> 本报告由 `surechembl_profile.py` 自动生成，耗时 203.3 秒（deep）。
