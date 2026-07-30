# ChEMBL 37 — SQLite 本地数据库

> 本目录记录 ChEMBL 37 SQLite 版数据库的下载、校验与结构信息，供 GKA_in_Brain 项目使用。

---

## 1. 下载信息

| 项目 | 内容 |
|---|---|
| 数据库版本 | **ChEMBL 37** |
| 官方数据准备日期 | 2026-05-01 |
| FTP 发布日期 | 2026-05-29 |
| **本地下载日期** | **2026-07-30** |
| 下载源 | `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz` |
| 下载方式 | `wget -c`（断点续传） |
| 本地路径 | `/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/` |
| 许可 | CC BY-SA 3.0（见 `LICENSE`、`REQUIRED.ATTRIBUTION`） |

### 目录内容

```
ChEMBL_37/
├── chembl_37_sqlite.tar.gz          5.76 GB   官方压缩包（已校验）
├── chembl_37/
│   └── chembl_37_sqlite/
│       ├── chembl_37.db            30.48 GB   ← SQLite 数据库本体
│       └── INSTALL_sqlite                     官方安装说明
├── chembl_37_release_notes.txt                官方 release notes
├── chembl_37_schema.pdf                       schema ER 图（PDF）
├── schema_documentation.txt                   schema 字段级文档（纯文本）
├── checksums.txt                              官方 SHA256 清单
├── LICENSE / README / REQUIRED.ATTRIBUTION    官方许可与引用要求
└── README_ChEMBL_37.md                        本文件
```

---

## 2. 官方 Checksum 与校验结果

官方 `checksums.txt` 中 SQLite 包的条目：

```
33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281  chembl_37_sqlite.tar.gz
```

本地校验（2026-07-30）：

```bash
$ echo "33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281  chembl_37_sqlite.tar.gz" | sha256sum -c -
chembl_37_sqlite.tar.gz: OK
```

**状态：✅ 校验通过**，文件大小 5,764,252,857 bytes，与官方一致。

同一 release 其他格式的官方 SHA256（备查，本地未下载）：

| 文件 | SHA256 |
|---|---|
| `chembl_37_mysql.tar.gz` | `e6eb871ee46121404d0caa79ebc00d451f5bd9b4b14688e223659cca216ec467` |
| `chembl_37_postgresql.tar.gz` | `d5fb08a6c9bf6cf7319dbdf74bdb5a7d4ce767aace670250ba1435cab4d58c58` |
| `chembl_37.sdf.gz` | `2cdf4d33c5f426130c610b8627c3702804e6f10d149d76b8a1bea674fa605a77` |
| `chembl_37_chemreps.txt.gz` | `ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7` |
| `chembl_37.fps.gz` | `c33bfac42abfea96840279ec35eeb3364c1063aa3e35c02369618beaf6389529` |
| `chembl_37.fa.gz` | `8f59596c4ee8f6cc7abcc59a4ea6f785ce428945de322fe9be6fb50808a7a9ae` |
| `chembl_37_bio.fa.gz` | `d258764963dd8edd1e76753c7e7fc3a228326dcb80e470e7073306292d9af6e7` |
| `chembl_37_hmmr.fa.gz` | `96e518cf3bb668c59a1c606a03d08034a8df5f736cf08292501df7f230b7c364` |

---

## 3. SQLite 数据库

| 项目 | 内容 |
|---|---|
| 文件 | `chembl_37/chembl_37_sqlite/chembl_37.db` |
| 大小 | 30,480,314,368 bytes（≈ 30.5 GB，解压后） |
| 官方构建用 SQLite 版本 | 3.26.0 |
| 表数量 | 72（含 `sqlite_stat1`） |

### 使用方式

```bash
# 命令行（需系统装有 sqlite3；本机当前未安装）
sqlite3 /ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db
```

```python
# Python（本机 python3 内置 sqlite3 3.31.1，已验证可读）
import sqlite3
DB = "/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db"
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)   # 只读打开，避免误写
cur = con.cursor()
cur.execute("SELECT chembl_id, pref_name FROM molecule_dictionary LIMIT 5")
print(cur.fetchall())
```

### 本地实测行数（与 release notes 完全一致，可作完整性旁证）

| 表 | 行数 |
|---|---|
| `activities` | 24,527,044 |
| `compound_records` | 3,824,604 |
| `molecule_dictionary` | 2,921,148 |
| `assays` | 1,970,438 |
| `docs` | 101,100 |
| `target_dictionary` | 18,552 |

---

## 4. Release Notes 摘要（ChEMBL 37）

完整原文见 `chembl_37_release_notes.txt`。

### 4.1 规模

数据准备日期 2026-05-01：

- 2,921,148 compounds（其中 2,897,819 有 mol file）
- 3,824,604 compound records（非唯一化合物）
- 24,527,044 activities
- 1,970,438 assays
- 18,552 targets
- 101,100 documents

主要 bioassay 数据来源（按 assay 数）：Scientific Literature（1,835,048）、Chemical Probe data from Literature（55,081）、EUbOPEN Chemogenomic Library Literature（23,222）、SureChEMBL Patent Bioactivity（16,624）、SGC Frankfurt Donated Chemical Probes（15,037）、BindingDB Patent Bioactivity（13,835）、PubChem BioAssays（2,999，但贡献 7,434,992 条 activity）。

### 4.2 本版更新的数据源

BindingDB Patent Bioactivity Data、British National Formulary (BNF)、Clinical Candidate Compounds、EUbOPEN Chemogenomic Library、INN、SGC Frankfurt Donated Chemical Probes、Scientific Literature、SureChEMBL Patent Bioactivity Data、USAN、WHO ATC。

新增数据集：12 个 DCP probe 数据集（`CHEMBL6193868`–`CHEMBL6193879`，如 THNAN69、JYQ-173、BAY-439、SM311、JNJ-79883960、RO4938581、ACBI3 等）。

### 4.3 数据变更与修订

- **专利数据去重**：BindingDB Patent（`src_id=37`）与 SureChEMBL Patent（`src_id=38`）识别出 13 组源内/源间重复文献，数据点数相同时保留最早条目，否则保留数据点更多的条目。
- **Chemical probe 标记**：`MOLECULE_DICTIONARY.CHEMICAL_PROBE` 依据 chemicalprobes.org 与 probes-drug.org（2026-04-02 抓取）更新。
- **新增 organism targets**：遗留的以生物体为对象的 assay（多为抗菌 assay）重新映射到具体 organism target，`TARGET_DICTIONARY` 新增约 400 条；`ASSAY_ORGANISM`/`ASSAY_TAX_ID`、`ORGANISM_CLASS` 按 NCBI taxonomy 对齐。
- **Assay 菌株注释**：约 7,300 个 assay 从 `DESCRIPTION` 中补全了 `ASSAY_STRAIN`。
- **组织首选名**：`TARGET_DICTIONARY` 与 `TISSUE_DICTIONARY` 的 tissue `PREF_NAME` 以 UBERON 名称统一。
- **靶向蛋白降解（TPD）数据**：新流程识别出约 29,000 条 TPD bioactivity，用新字段 `ACTIVITIES.MODALITY = "Targeted protein degradation"` 标注；约 1,200 个遗留 TPD assay（约 3,000 条 activity）的 target 映射被修正；PROTEIN-PROTEIN INTERACTION 类型 target 的 `PREF_NAME` 按 UniProt 推荐名标准化（格式：降解效应蛋白/致病蛋白）。该部分仍在迭代中。
- **结构整理**：合并约 400 组 InChI 未识别、但经 SMILES tautomer hash 发现的互变异构体重复；整理 >3,000 条 BindingDB 遗留结构。
- **药物与临床候选**：ChEMBL 36 因 ATC/BNF 源调整导致 approved drug 数量大幅变化；本版经人工审阅将其中 169 个化合物恢复为 `MAX_PHASE = 4`（判定标准：至少一个司法辖区有监管批准证据，如 FDA/EMA/MHRA/PMDA；仅文献证据不足以支撑）。该审阅仍在持续，版本间 ATC `MAX_PHASE` 差异仍会存在。
- **`STANDARD_TEXT_VALUE`**：`ACTIVITIES`、`ACTIVITY_PROPERTIES`、`ASSAY_PARAMETERS` 中的高频 text_value 现按字典法自动标准化。

---

## 5. Schema

参考文件：

- `chembl_37_schema.pdf` — 官方 ER 图
- `schema_documentation.txt` — 逐表逐字段说明（含主键/外键/类型/注释）
- 在线：https://chembl.gitbook.io/chembl-interface-documentation/db-schema-description

### 5.1 相对 ChEMBL 36 的 schema 变更

**新增字段**

- `ACTIVITIES.MODALITY` — 标注化合物设计模态（如 "Targeted Protein Degradation"）。与 `ACTION_TYPE` 不同，`MODALITY` 记录的化合物对靶点可能是活性也可能是无活性。本版初始仅承载 TPD 注释结果。

**已废弃的表**

- `CURATION_LOOKUP`
- `PROTEIN_CLASS_SYNONYMS`

**已废弃的字段**

- `CELL_DICTIONARY.CL_LINCS_ID`
- `ASSAYS.CURATED_BY`

### 5.2 GKA 项目常用核心表与连接路径

```
compound_structures ──molregno── molecule_dictionary ──molregno── compound_records
                                        │                              │ record_id
                                        │ molregno                     │
                                        └──────── activities ──────────┘
                                                   │ assay_id      │ doc_id
                                                 assays            docs
                                                   │ tid
                                          target_dictionary ──tid── target_components
                                                                          │ component_id
                                                                  component_sequences
```

典型查询骨架（按靶点取带 pChEMBL 的活性数据）：

```sql
SELECT md.chembl_id, cs.canonical_smiles, act.standard_type,
       act.standard_relation, act.standard_value, act.standard_units,
       act.pchembl_value, td.pref_name
FROM activities act
JOIN assays a               ON act.assay_id = a.assay_id
JOIN target_dictionary td   ON a.tid        = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON md.molregno  = cs.molregno
WHERE td.chembl_id = 'CHEMBL3820'          -- 人 Hexokinase-4 / glucokinase (GCK)
  AND act.pchembl_value IS NOT NULL
  AND act.standard_relation = '='
  AND a.confidence_score >= 8
  AND act.data_validity_comment IS NULL
  AND act.potential_duplicate = 0;
```

### 5.3 全部表清单（72）

`ACTION_TYPE`, `ACTIVITIES`, `ACTIVITY_PROPERTIES`, `ACTIVITY_SMID`, `ACTIVITY_STDS_LOOKUP`, `ACTIVITY_SUPP`, `ACTIVITY_SUPP_MAP`, `ASSAY_CLASS_MAP`, `ASSAY_CLASSIFICATION`, `ASSAY_PARAMETERS`, `ASSAY_TYPE`, `ASSAYS`, `ATC_CLASSIFICATION`, `BINDING_SITES`, `BIO_COMPONENT_SEQUENCES`, `BIOASSAY_ONTOLOGY`, `BIOTHERAPEUTIC_COMPONENTS`, `BIOTHERAPEUTICS`, `CELL_DICTIONARY`, `CHEMBL_ID_LOOKUP`, `CHEMBL_RELEASE`, `COMPONENT_CLASS`, `COMPONENT_DOMAINS`, `COMPONENT_GO`, `COMPONENT_SEQUENCES`, `COMPONENT_SYNONYMS`, `COMPOUND_PROPERTIES`, `COMPOUND_RECORDS`, `COMPOUND_STRUCTURAL_ALERTS`, `COMPOUND_STRUCTURES`, `CONFIDENCE_SCORE_LOOKUP`, `DATA_VALIDITY_LOOKUP`, `DEFINED_DAILY_DOSE`, `DOCS`, `DOMAINS`, `DRUG_INDICATION`, `DRUG_MECHANISM`, `DRUG_WARNING`, `FORMULATIONS`, `GO_CLASSIFICATION`, `INDICATION_REFS`, `LIGAND_EFF`, `MECHANISM_REFS`, `METABOLISM`, `METABOLISM_REFS`, `MOLECULE_ATC_CLASSIFICATION`, `MOLECULE_DICTIONARY`, `MOLECULE_HIERARCHY`, `MOLECULE_SYNONYMS`, `ORGANISM_CLASS`, `PATENT_USE_CODES`, `PESTICIDE_CLASS_MAPPING`, `PESTICIDE_CLASSIFICATION`, `PREDICTED_BINDING_DOMAINS`, `PRODUCT_PATENTS`, `PRODUCTS`, `PROTEIN_CLASSIFICATION`, `RELATIONSHIP_TYPE`, `SITE_COMPONENTS`, `SOURCE`, `SQLITE_STAT1`, `STRUCTURAL_ALERT_SETS`, `STRUCTURAL_ALERTS`, `TARGET_COMPONENTS`, `TARGET_DICTIONARY`, `TARGET_RELATIONS`, `TARGET_TYPE`, `TISSUE_DICTIONARY`, `USAN_STEMS`, `VARIANT_SEQUENCES`, `VERSION`, `WARNING_REFS`

---

## 6. 引用要求

ChEMBL 数据遵循 `LICENSE` 中的协议（CC BY-SA）。按 `REQUIRED.ATTRIBUTION` 要求，使用 ChEMBL 数据发表时应引用：

> Mendez D, Gaulton A, Bento AP, Chambers J, De Veij M, Félix E, Magariños MP, Mosquera JF, Mutowo P, Nowotka M, Gordillo-Marañón M, Hunter F, Junco L, Mugumbate G, Rodriguez-Lopez M, Atkinson F, Bosc N, Radoux CJ, Segura-Cabrera A, Hersey A, Leach AR. **ChEMBL: towards direct deposition of bioassay data.** *Nucleic Acids Res.* 2019;47(D1):D930–D940. DOI: [10.1093/nar/gky1075](https://doi.org/10.1093/nar/gky1075)

若将 ChEMBL 并入其他工作，应保留 ChEMBL ID，并明确标注所用的 release 号（本处为 **ChEMBL 37**）。

联系方式：chembl-help@ebi.ac.uk

---

*文档生成日期：2026-07-30*
