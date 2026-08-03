# Fridén 2009 Kp,uu,brain 数据快照

## ⚠ 这是本步对照集的**唯一**数据来源

与同目录下的 B3DB 不同：B3DB 只用于记录判据差异，**这份才是对照集的来源**。

## 原始文献

Fridén M, Winiwarter S, Jerndal G, Bengtsson O, Wan H, Bredberg U, et al.
*Structure-brain exposure relationships in rat and human using a novel data set of
unbound drug concentrations in brain interstitial and cerebrospinal fluids.*
**J Med Chem. 2009 Oct 22;52(20):6233-43.**
PMID: **19764786** · DOI: **10.1021/jm901036q**

⚠ **原文非开放获取**（ACS 期刊，PubMed Central 中只有引用它的文章，没有它本身），
补充材料无法直接下载。**ChEMBL 是这批数据的公开获取途径。**

## 获取路径

```
ChEMBL 37  →  assay CHEMBL1798466  (doc CHEMBL1795238)
```

| | |
|---|---|
| 实验体系 | Sprague-Dawley 大鼠（NCBI tax 10116） |
| 给药 | 2/3 化合物盒式给药，4 h 恒速静脉输注，1 (ml/kg)/hr ≈ 2 (µmol/kg)/hr |
| 测量量 | `K(p,uu,brain)` = 游离脑间质浓度 / 游离血浆浓度，**无量纲** |
| `standard_type` | 全部单一为 `K(p,uu,brain)`，无大小写变体 |
| 记录数 | **42 行 / 42 个不重复分子**，无重复 |
| 数值缺失 | 1 个（Cefotaxime，`standard_value` 为空） |
| `data_validity_comment` | 全空（无质量标记） |
| `potential_duplicate` | 全部为 0 |

## 快照文件

```
Friden2009_CHEMBL1798466_raw.tsv     35,992 B     42 行 × 36 列
sha256  8926fad7dceb502e1c2d63b7d86a078c5887de06c4cbe1e78ec1cd078e6fc464
```

**未经任何处理的原始导出**——分组、打标、自检都在
`Step3_00_Build_BBB_Control_Set.py` 里做，快照保持 ChEMBL 原值。
含 assay/doc 全部著录项、活性原值与标准化值、`compound_key`、
SMILES / InChI / InChIKey。

导出日期：**2026-08-02**，ChEMBL **37**。

## 为什么要冻结这份快照

1. **原文非开放获取**，快照是这批数据在本项目内的可追溯凭证。
2. **不必依赖 29 GB 的 ChEMBL 库**即可复现 Step3_00 —— 脚本会在
   `--db` 路径不存在时自动回落到快照。
3. ChEMBL 后续版本可能修订这些记录，**任何结论都对应 ChEMBL 37 这一版**。

## 两条路径等价性

脚本支持两条来源，已验证**产物逐字节一致**：

```bash
# 走 ChEMBL 库（默认，库存在时）
python Step3_00_Build_BBB_Control_Set.py

# 强制走快照（库不存在时自动，也可显式指定）
python Step3_00_Build_BBB_Control_Set.py --use-raw
```

两处排序都加了 `molecule_chembl_id` 作次级键——
Moxalactam 与 Nelfinavir 的 Kp,uu 都是 0.019，
不指定次级键时 SQL 与 pandas 的并列顺序不同，会导致 `control_id` 错位。

运行时校验 SHA256，对不上打印警告但不中止。
