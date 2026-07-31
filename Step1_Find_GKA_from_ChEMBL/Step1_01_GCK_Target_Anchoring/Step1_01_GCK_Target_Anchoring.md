# Step1_01 GCK 靶点锚定结果

- 锚点 UniProt accession：**P35557**
- ChEMBL 版本：**CHEMBL_37（2026-05-01）**
- 数据库文件：`/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db`
- 运行时间：2026-07-30 21:23:17
- 命中 target 数：**2**

## 蛋白组件

- `component_id` = **2138**（PROTEIN）
- 描述：Hexokinase-4
- 物种：Homo sapiens（tax_id 9606）
- 序列来源：SWISS-PROT 2026_01

> 注：这里的版本取自 `component_sequences.db_version`（逐条记录）。ChEMBL 的 `version` 表另有一条全局声明，两者在 ChEMBL 37 中并不一致（全局声明为 Swiss-Prot 2025_03，而绝大多数组件记录标注为 2026_01）。以逐条记录的值为准。

## 命中的 ChEMBL target

| target_chembl_id | tid | pref_name | organism | target_type | homologue | assays | assays(conf≥8) | activities | 有 pChEMBL |
| --- | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `CHEMBL3885579` | 117123 | Glucokinase/Glucokinase regulatory protein | Homo sapiens | PROTEIN-PROTEIN INTERACTION | 0 | 1 | 0 | 40 | 37 |
| `CHEMBL3820` | 20095 | Hexokinase-4 | Homo sapiens | SINGLE PROTEIN | 0 | 227 | 227 | 3,222 | 1,389 |

> `homologue = 0` 表示该组件是这个 target 的直接组成部分，而非同源映射。

## 说明

同一个蛋白组件会挂在多个 target 上：既有把它单独作为作用对象的 `SINGLE PROTEIN` 靶点，也有把它与其他蛋白的相互作用作为作用对象的 `PROTEIN-PROTEIN INTERACTION` 靶点。下游取活性数据时要**显式决定保留哪些**，不能默认只有一个。

