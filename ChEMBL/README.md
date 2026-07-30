# chembl_profile.py

调查一个 ChEMBL SQLite 数据库，生成 Markdown 概览报告。报告面向**有生物/生信背景、但没有药化背景**的读者：既有统计数字，也解释 assay、IC50、pChEMBL、confidence_score 这些概念。

只依赖 Python 3 标准库，以只读方式打开数据库。

## 运行

```bash
cd /ShangGaoAIProjects/GKA_in_Brain/ChEMBL

# 完整报告（约 70 秒）
python3 chembl_profile.py ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db \
        -o chembl_37_profile_report.md
```

## 参数

| 参数 | 说明 |
|---|---|
| `db_path` | 必填，ChEMBL SQLite 文件路径 |
| `-o, --out` | 输出路径（默认：数据库同目录下 `chembl_profile_report.md`） |
| `--quick` | 快速模式，约 18 秒。行数取自统计表为估计值，跳过 activities 全表统计 |
| `--deep` | 加做 activities×assays×targets 大连接，多出 3 张靶点排行对比表（推荐） |
| `--top N` | 各排行榜取前 N 项，默认 20 |
| `--stdout` | 报告同时打印到标准输出 |

## 说明

- 适用于 ChEMBL 33+ 各版本：脚本先探测表/列是否存在，缺失的章节自动跳过。
- 当前报告：`chembl_37_profile_report.md`。
- 数据库本体（~30 GB）与压缩包不入 git，下载与校验方式见 `ChEMBL_37/README_ChEMBL_37.md`。
