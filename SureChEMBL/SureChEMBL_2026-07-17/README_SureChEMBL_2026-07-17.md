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
├── SHA256SUMS.txt                             本地 SHA256 清单
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

**下载完整性：✅ 9/9 文件字节数与服务端完全一致**（2026-07-31 核对）。

本地 SHA256（`sha256sum` 生成，同目录 `SHA256SUMS.txt` 可直接 `sha256sum -c`）：

| 文件 | SHA256 |
|---|---|
| `patents.parquet` | `7d1842bf6f0a086e9efe79b4060e14c6e26f74a367dca44e3bf08418eedf0054` |
| `patent_compound_map.parquet` | `2a67f7a843737b57bb61f82a878ffaf357ac78996569cfe04aeeba6ebbc19807` |
| `compounds.parquet` | `b4156878e42749fa5376a11df4e6861785a638183db820cbbb3c59ce8b29c5e1` |
| `biomedical_locations.parquet` | `c6256a2f545696fb00a97c75d5f91a4dc8f1fbd03891958a8921cca897191542` |
| `fpsim2_fingerprints.h5` | `5f63cbcace9afbb9e9b654822ad76cb978d10a01fbb6e5a8db661b5da41b6d60` |
| `biomedical_entities.parquet` | `b8d9a0e3757f5a9b52c8753037979f650a8eecaa81d0bbfa9880004040a99504` |
| `biomedical_types.parquet` | `b2095c9fd00165866e0b0aaf1286de3de53713d20e62bb6be29f234bc1405861` |
| `fields.parquet` | `4ca5b062eafd1ee41614530c145026189512242251c25e17e6f185c99300d610` |
| `LICENCE` | `9ba9550ad48438d0836ddab3da480b3b69ffa0aac7b7878b5a0039e7ab429411` |

### 下载踩过的坑

- **EBI 的连接会中途掉**（实测 `OpenSSL SSL_read: No route to host, errno 113`），
  必须断点续传 + 多轮补传，`.download.sh` 已处理。
- **本机 curl 是 7.68.0，不支持 `--retry-all-errors`**（7.71 才加入）。
  用了会直接报 `option is unknown` 且每轮空转，一个字节都下不来。
- **别用后台任务的退出码判断下载成功**——管道末端 `tail`/`echo` 会把退出码吃掉，
  报 0 但实际失败。**只认字节数核对**。

---

## 3. 数据规模

### 3.1 本地实测（以此为准）

| 表 | 行数 |
|---|---:|
| `patent_compound_map` | **1,537,106,020** |
| `biomedical_locations` | 453,327,904 |
| `patents` | 44,912,542 |
| `compounds` | 30,990,818 |
| `biomedical_entities` | 1,059,724 |
| `fields` | 6 |
| `biomedical_types` | 4 |

### 3.2 ⚠ 与官方宣传数字对不上

官方博客/新闻稿称 SureChEMBL 2.0 有 **~116.6 M 专利文档、~47.7 M 独立化合物**，
但**本地这份 bulk 快照实测只有 44.9 M 专利、31.0 M 化合物**。

差异很可能是「系统里索引过的文档总数」与「bulk 导出里带化学标注的文档数」口径不同，
官方文档没有说明。**做任何规模陈述都以上面 3.1 的实测行数为准，不要引用宣传数字。**

官方公布的专利局分布（备查，未在本地核实）：CNIPA 51.1 M 专利 / 3.8 M 化合物，
USPTO 21.2 M 专利 / 19.9 M 化合物。**专利数不代表化学信息量**，CNIPA 专利最多但
抽出的化合物最少。

---

## 4. 表结构

六张数据表 + 两张字典，schema 镜像 SureChEMBL 内部关系库，**不做冗余展开**
（这是 2.0 相对旧版 MAP files 的主要改进，文件体积因此小很多）。

以下 schema 均为**本地实测**（`pq.ParquetFile(...).schema_arrow`），不是抄文档。

### 4.1 `compounds`（30,990,818 行，2,385 个 row group）

```
id:int64 | smiles:large_string | inchi:large_string
inchi_key:large_string | mol_weight:double
```

### 4.2 `patents`（44,912,542 行，1,000 个 row group）

```
id:int64 | patent_number:large_string | country:large_string
publication_date:date32[day] | family_id:int64
cpc:list<string> | ipcr:list<string> | ipc:list<string> | ecla:list<string>
assignee:list<string> | title:string
```

**分类号与申请人是 list 类型**，不是字符串——过滤要用 `list_contains` 之类的数组函数，
`LIKE` 匹配不上。`publication_date` 是原生 date32，可直接比较。

### 4.3 `patent_compound_map`（**1,537,106,020 行**，5,006 个 row group）

```
patent_id:int64 | compound_id:int64 | field_id:int64
```

**15.4 亿行，这是全库最大的表。** 绝不能整表读进内存，必须用 duckdb 做谓词下推，
先按 `compound_id` 或 `patent_id` 过滤再取。

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

### 4.7 `biomedical_locations`（453,327,904 行，433 个 row group）

```
entity_id:int64 | patent_id:int64 | field_id:int64 | count:int64
```

实体 ↔ 专利 ↔ 部分 的映射，`count` 是该实体在该部分出现的次数。同样是亿级表，
查询前必须先按 `entity_id` 收窄。

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

实测 glucokinase（全部 33 个 HGNC:4195 实体）的分布，差距是 11 倍：

| field | 专利数 | 提及次数 |
|---|---:|---:|
| `desc` 说明书 | **28,272** | 123,965 |
| `clms` 权利要求 | **2,597** | 10,442 |
| `abst` 摘要 | 1,462 | 3,189 |
| `ttl` 标题 | 1,085 | 1,925 |

**2.8 万篇专利提到 glucokinase，但只有 2,597 篇在权利要求里提。**
不区分 `field_id` 就会把一万倍于真实规模的噪声当成 GKA 专利。

### 5.2 ⚠ 锚定要用 `resolved_form`，不是文本匹配 `original_text`

**`resolved_form = 'HGNC:4195'` 对应 33 个不同的 `original_text`**，
按字符串搜 `glucokinase` 只能捞到其中一部分。33 个里包含：

- **OCR 破碎形**：`glucokmase`(80 篇)、`giucokinase`(51)、`glu- cokinase`(37)、
  `gluco- kinase`(20)、`Glucoki nase`、`gl uc ok i na se`、`Glucokinas e` …
  专利全文是扫描件 OCR 出来的，这类变体**只能靠 `resolved_form` 兜住**。
- **ChEMBL 那边的名字**：`Hexokinase 4`(103 篇)、`Hexokinase-4`(25 篇)——
  与 CLAUDE.md 里「GCK 在 ChEMBL 的 `pref_name` 是 Hexokinase-4」对上了。
- **缩写**：`GK`(1,343 篇)、`GKAs`(131)、`GKA`(13)、`GlcK`(48)、`glk`(7)、
  `gki`(7)、`gukA`(6)、`GlkA`(8)、`Hk4`(2)、`GluK`(1)、`GcK`(11)。

### 5.3 ⚠ 缩写形是假阳性来源，尤其 `GK`

| surface form | 专利数 | 风险 |
|---|---:|---|
| `glucokinase`（全称） | 29,191 | 可靠 |
| **`GK`** | **1,343** | **⚠ 最大风险**。`GK` 在糖尿病文献里更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型），也可能是 glycerol kinase。**恰好和我们的领域重叠，必须逐篇核** |
| `4` | 4 | **标注错误**：单个数字 `4` 被解析成了 glucokinase（entity id 528883） |
| `GlkA` / `glk` / `gukA` | 6–8 | 细菌 glucokinase 基因名，不是人 GCK |

另外这些**不**解析到 HGNC:4195，别混进来：

| entity | resolved_form | 是什么 |
|---|---|---|
| `Glucokinase regulator` | `HGNC:4196` | GCKR / GKRP，**不是 GCK** |
| `GKRP` | `Q14397` | 同上 |
| `GCK`（全大写，id 21805） | **（空）** | **未解析**——无法判断是 glucokinase 还是 MAP4K2 |
| `HPK/GCK-like kinase` | `HGNC:6866` | MAP4K 家族，**不是 glucokinase** |
| `GCK-1` | `H2L099` | 别的蛋白 |

注意大小写：全大写 `GCK` 未解析，但 `GcK`(id 788019) 解析到了 HGNC:4195。
**别拿 `GCK` 当锚点**——这正是 CLAUDE.md 里 ChEMBL 误挂 52 个 MAP4K2 assay 的同一个坑，
而且这里更隐蔽：ChEMBL 至少给了个错靶点能查，这里直接是空值。

**实测的取舍**：在 `clms` 层，只用全称 `glucokinase` 得 2,553 篇，
用全部 33 个 surface form 得 2,597 篇——**多出来的 32 个形只贡献 +44 篇（+1.7%），
却带进 `GK` 那 1,343 篇的歧义**。做权利要求级检索时建议只用全称 + OCR 变体，
把纯缩写单独拉一张表人工核。

### 5.4 ⚠ `resolved_form` 不是单一命名空间

同一列里混着 `HGNC:4195`（HGNC，带前缀）、`Q14397`(UniProt，无前缀)、
`H2L099`（UniProt，无前缀）、空串。
**解析前要先判断命名空间**，不能直接当 UniProt accession 用。
我们在 ChEMBL 侧的锚点是 UniProt **P35557**，与 HGNC:4195 是同一基因，
跨库对齐时要走基因层而不是字符串比对。

### 5.5 没有 `glucokinase activator` 这样的实体

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
