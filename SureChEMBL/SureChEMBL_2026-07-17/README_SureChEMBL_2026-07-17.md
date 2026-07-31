# SureChEMBL 2.0 — 2026-07-17 全量快照

> 本目录记录 SureChEMBL 专利化学数据的下载、校验与结构信息，供 GKA_in_Brain 项目使用。

---

## 1. 下载信息

| 项目 | 内容 |
|---|---|
| 系统版本 | **SureChEMBL 2.0**（2025-05-07 发布的重构版） |
| 数据发布日期 | **2026-07-17**（FTP 上线 2026-07-29） |
| **本地下载日期** | **2026-07-31** |
| 下载源 | `https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/2026-07-17/` |
| 下载方式 | `curl -C -`（断点续传），脚本见本目录 `.download.sh` |
| 本地路径 | `/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17/` |
| 许可 | **CC BY 4.0**（见 `LICENCE`） |

**⚠ 固定日期目录，不要用 `latest/`。** SureChEMBL **每两周发一版全量快照**（不是增量），
`latest/` 两周就会被覆盖。本项目锚定 `2026-07-17` 这一版，任何结论都对应它。

### 目录内容

```
SureChEMBL_2026-07-17/
├── patents.parquet                  5.51 GB   专利元数据
├── patent_compound_map.parquet      4.63 GB   专利 ↔ 化合物关联（带出现位置）
├── compounds.parquet                3.91 GB   化合物结构
├── biomedical_locations.parquet     1.60 GB   生物医学实体 ↔ 专利/部分 的位置与次数
├── fpsim2_fingerprints.h5           1.26 GB   FPSim2 指纹库，可直接做相似性检索
├── biomedical_entities.parquet     31.66 MB   生物医学实体词表
├── biomedical_types.parquet          3.1 KB   实体类型字典（4 行）
├── fields.parquet                    1.7 KB   专利部分字典（6 行）
├── LICENCE                          18.2 KB   CC BY 4.0
├── .server_manifest.tsv                       下载前抓的服务端 Content-Length / ETag
├── .download.sh                               断点续传下载脚本（重跑安全）
└── README_SureChEMBL_2026-07-17.md            本文件
```

合计 **16.94 GB**（18,190,499,097 bytes）。

---

## 2. 校验

EBI **没有为 bulk_data 提供官方 checksums 文件**（与 ChEMBL 的 `checksums.txt` 不同）。
因此本地采用两级校验：

1. **大小比对**：下载前先 `curl -I` 抓服务端 `Content-Length` 与 `ETag`，
   存入 `.server_manifest.tsv`；下载后逐文件比对字节数。
2. **本地 SHA256**：下载完成后计算并记录在下表，供日后确认文件未被改动。

服务端元数据（下载前抓取，2026-07-31）：

| 文件 | Content-Length | Last-Modified (UTC) |
|---|---:|---|
| `patents.parquet` | 5,913,041,063 | Fri, 17 Jul 2026 20:52:19 |
| `patent_compound_map.parquet` | 4,976,853,813 | Fri, 17 Jul 2026 20:37:32 |
| `compounds.parquet` | 4,199,702,118 | Fri, 17 Jul 2026 16:42:58 |
| `biomedical_locations.parquet` | 1,714,784,809 | Sat, 18 Jul 2026 07:24:33 |
| `fpsim2_fingerprints.h5` | 1,352,894,071 | Fri, 17 Jul 2026 16:03:42 |
| `biomedical_entities.parquet` | 33,201,669 | Fri, 17 Jul 2026 15:57:18 |
| `biomedical_types.parquet` | 3,179 | Fri, 17 Jul 2026 15:57:28 |
| `fields.parquet` | 1,718 | Fri, 17 Jul 2026 15:57:53 |
| `LICENCE` | 18,657 | Wed, 29 Jul 2026 16:25:24 |

本地 SHA256：**（下载完成后填写）**

---

## 3. 数据规模（SureChEMBL 2.0 官方公布）

| 项目 | 数量 |
|---|---|
| 专利文档 | ~116.6 M |
| 独立化合物 | ~47.7 M |

覆盖 **5 家专利局**：

| 专利局 | 专利数 | 独立化合物 |
|---|---|---|
| CNIPA（中国，2.0 新增） | 51.1 M | 3.8 M |
| USPTO（美国） | 21.2 M | **19.9 M** |
| WIPO / EPO / JPO | 其余 | — |

注意 CNIPA 专利数最多但抽出的化合物最少，USPTO 反之——**专利数不代表化学信息量**。

---

## 4. 表结构

六张数据表 + 两张字典，schema 镜像 SureChEMBL 内部关系库，**不做冗余展开**
（这是 2.0 相对旧版 MAP files 的主要改进，文件体积因此小很多）。

### 4.1 `compounds`

化合物结构：`id`、`smiles`、`inchi`、`inchi_key`、`mol_weight`。

### 4.2 `patents`

专利元数据：`id`、`patent_number`、`country`、`publication_date`、`family_id`、
`cpc`、`ipcr`、`ipc`、`ecla`、`assignee`、`title`。

### 4.3 `patent_compound_map`

专利 ↔ 化合物：`patent_id`、`compound_id`、**`field_id`**。

### 4.4 `fields`（6 行，已实测）

| id | field_name | 含义 |
|---:|---|---|
| 1 | `desc` | 说明书 |
| 2 | `clms` | **权利要求** |
| 3 | `abst` | 摘要 |
| 4 | `ttl` | 标题 |
| 5 | `image` | 图片 |
| 6 | `molattachment` | MOL 附件 |

### 4.5 `biomedical_entities`（1,059,724 行，已实测）

`id`、`type_id`、`original_text`、`corrected_text`、`resolved_form`。

### 4.6 `biomedical_types`（4 行，已实测）

| id | type_name | 来源 | 本快照实际条目数 |
|---:|---|---|---:|
| 1 | `GeneOrProtein` | UniProt、HGNC、Entrez Gene | 913,773 |
| 2 | `Disease` | MeSH、Disease Ontology、Wikipedia | 137,749 |
| 3 | `Mechanism` | 化学物质的作用/机制 | 8,202 |
| 4 | `Physquant` | 各类物理量 | **0** |

### 4.7 `biomedical_locations`

实体 ↔ 专利 ↔ 部分 的映射，带出现次数。

### 4.8 `fpsim2_fingerprints.h5`

FPSim2 格式指纹库，可直接对 4,770 万化合物做相似性检索，不必自己算指纹。

---

## 5. 已验证的陷阱（做检索前必读）

### 5.1 ⚠ `field_id` 决定语义：出现在权利要求里 ≠ 出现在说明书里

`patent_compound_map` 的 `field_id` 指出化合物出现在专利的哪个部分。
**出现在 `clms`（权利要求，`field_id = 2`）的化合物是这篇专利要保护的；
出现在 `desc`（说明书）里的可能只是背景技术引用的他人化合物。**
不区分就会把「专利里提到过」当成「这篇专利的化合物」。

这与 Step1_03 的教训同构：**位置/归属不等于目的**，参见 CLAUDE.md 里
「挂在靶点下的 assay 未必真在测该靶点」。

### 5.2 ⚠ 「GCK」歧义在这里同样存在，而且更隐蔽

`biomedical_entities` 里实测：

| entity | resolved_form | 说明 |
|---|---|---|
| `glucokinase` | `HGNC:4195` | ✅ 就是我们的靶点 GCK |
| `Glucokinase regulator` | `HGNC:4196` | GCKR，即 GKRP，**不是 GCK** |
| `GKRP` | `Q14397` | 同上，UniProt 号 |
| `GCK` | **（空）** | **未解析！** 无法判断是 glucokinase 还是 MAP4K2 |
| `GCK-1` | `H2L099` | 解析到了别的蛋白 |
| `HPK/GCK-like kinase` | `HGNC:6866` | MAP4K 家族，**不是 glucokinase** |
| `glucokinase regulatory protein` | （空） | 未解析 |

**缩写 `GCK` 在 SureChEMBL 里是未解析状态**，比 ChEMBL 那次（误挂 52 个 MAP4K2 assay）
更难发现——ChEMBL 至少给了个错的靶点可以查，这里直接是空。
**锚定必须用 `glucokinase` / `HGNC:4195`，不能用缩写 `GCK`。**

### 5.3 ⚠ `resolved_form` 不是单一命名空间

同一列里混着 `HGNC:4195`（HGNC，带前缀）、`Q14397`(UniProt，无前缀)、
`H2L099`（UniProt，无前缀）、空串。
**解析前要先判断命名空间**，不能直接当 UniProt accession 用。
我们在 ChEMBL 侧的锚点是 UniProt **P35557**，与 HGNC:4195 是同一基因，
跨库对齐时要走基因层而不是字符串比对。

### 5.4 没有 `glucokinase activator` 这样的实体

`biomedical_entities` 里搜不到「glucokinase activator」——
`Mechanism` 类型只有 8,202 条，粒度不足以直接给出「GKA」这个概念。
**要找 GKA 专利，得靠「化合物结构相似性」或「glucokinase 实体 + 化学结构」组合，
不能指望一个现成的机制标签。**

---

## 6. 使用方式

环境：micromamba `GKA_in_Brain`（见项目 CLAUDE.md），已装 `pyarrow` 25.0.0 与 `duckdb` 1.5.5。

```python
# 小表直接读
import pyarrow.parquet as pq
t = pq.read_table("fields.parquet")

# 大表不要 read_table 全量加载（patents.parquet 5.5 GB）——用 duckdb 直接查
import duckdb
duckdb.sql("""
    SELECT country, COUNT(*) FROM 'patents.parquet' GROUP BY country ORDER BY 2 DESC
""").show()
```

```bash
# 重新下载 / 断点续传（重跑安全，已完整的文件自动跳过）
./.download.sh
```

---

## 7. 引用

SureChEMBL 数据按 **CC BY 4.0** 使用，需署名。原始论文：

> Papadatos G, Davies M, Dedman N, et al. SureChEMBL: a large-scale, chemically
> annotated patent document database. *Nucleic Acids Research*. 2016;44(D1):D1220-8.

---

## 8. 已废弃的下载方式（别再照旧教程找）

| 方式 | 状态 |
|---|---|
| MAP files（季度 TSV） | **已废弃**，最后一版停在 2023 |
| Compound data dumps（季度 SDF/TXT） | **已废弃** |
| Data Client（私有 FTP 日更） | **已废弃**，最后更新 2024-06 |

2.0 起统一为本目录这套两周一次的 Parquet 全量快照。
