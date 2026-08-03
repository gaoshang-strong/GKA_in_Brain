# B3DB 快照说明

## ⚠ 这不是本步对照集的数据来源

本步的对照集来自 **ChEMBL 37 assay `CHEMBL1798466`**（Fridén et al. 2009 的
`K(p,uu,brain)` 实测值），与本目录下的 B3DB 文件无关。

B3DB 在这里只有一个用途：**记录 `logBB` 与 `Kp,uu` 两个判据在同一批化合物上的差异**，
产物是 `Step3_00_B3DB_Comparison.csv`。

**不对 B3DB 的标签作对错判定。** 两者测的不是同一个量
（`logBB` 是脑/血浆总浓度比，`Kp,uu` 是游离浓度比），
分组不同是预期之内的；判据取决于问题，差异本身留给下游判断。

## 来源

| | |
|---|---|
| 数据集 | B3DB (Blood-Brain Barrier Database) |
| 仓库 | https://github.com/theochem/B3DB |
| 论文 | Meng F, Xi Y, Huang J, Ayers PW. *A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors.* **Sci Data 2021;8:289.** DOI: 10.1038/s41597-021-01069-5 |
| 许可 | 见上游仓库 |
| 下载日期 | **2026-08-02** |

## 文件与校验

```
B3DB_classification.tsv    2,635,683 B   7,807 化合物，BBB+/BBB- 分类标签
B3DB_regression.tsv          378,224 B   1,058 化合物，数值 logBB
```

```
sha256  47a160aea1551423ead2aab70b301fb1646aea03cf151fd614eb37c78fa55b82  B3DB_classification.tsv
sha256  1be4e33ab1fa1d99897541b6e0a9a00cd6061cf970a36dc838e640407735eba1  B3DB_regression.tsv
```

`Step3_00_Build_BBB_Control_Set.py` 里写死了 classification 的 SHA256（`B3DB_SHA256`），
运行时校验，**对不上会打印警告但不中止**——因为快照变更不影响对照集本身，
只影响差异比对那一张表。

## 重新下载

```bash
curl -sL -o B3DB_classification.tsv \
  https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_classification.tsv
curl -sL -o B3DB_regression.tsv \
  https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_regression.tsv
sha256sum B3DB_classification.tsv B3DB_regression.tsv
```

⚠ **上游会不定期追加新数据**（仓库自述 "occasionally uploaded with new experimental data"），
`main` 分支不是固定版本。SHA256 对不上说明快照已变更，
比对结果可能随之改变——与 SureChEMBL「每两周覆盖 `latest/`」是同一类问题，
**任何结论都对应 2026-08-02 这一版**。

## 本地实测的结构（与上游自述一致）

| | |
|---|---|
| 分类表规模 | 7,807 行 |
| 其中有数值 `logBB` | **1,058（13.6%）** |
| 其余 6,749 | 仅有继承自上游文献的分类标签，无数值可核对 |
| 标签分布 | BBB+ 4,956 / BBB− 2,851 |
| `group` 列 | A 1,058（有数值）/ B 3,621 / C 3,077 / D 51 |

`reference` 列里的 `R1`–`R50` 是 50 个原始数据集的编号，
对照表见上游 `raw_data/raw_data_summary.tsv`。
