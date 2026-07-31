从 SureChEMBL 2026-07-17 快照里检索**所有与葡萄糖激酶（GCK）相关的专利**，
作为 Step2 的锚定表。

这一步是**召回优先**：先把网撒全，不判断方向（激活/抑制）、不判断是不是 GKA 专利。
**只用 bulk parquet 的实体标注按规则检索，不调 API、不用 LLM 解析文本。**

输入：`SureChEMBL/SureChEMBL_2026-07-17/`
     `biomedical_entities.parquet` / `biomedical_locations.parquet` /
     `patents.parquet` / `patent_compound_map.parquet` / `compounds.parquet` / `fields.parquet`

交叉引用：`Step1_06_GKA_Physicochemical_Properties.csv`（782 个 ChEMBL 候选，用来交叉与自检）

## 一、锚定规则

**核心：按 `resolved_form` 锚定，不做文本匹配。**
专利全文来自扫描件 OCR，`glucokinase` 有几十种破碎写法
（`glucokmase`、`giucokinase`、`gl uc ok i na se`…），文本匹配一定漏。

### 纳入

| 规则 | 说明 | 预期文档数 |
|---|---|---:|
| `resolved_form = 'HGNC:4195'` | 人 GCK 基因，覆盖 33 个 surface form | 29,332 |
| 白名单（未解析但明确是 GCK） | `glucokinase (hexokinase 4)` | 102 |
| 白名单 | `Glucokinase (hexokinase 4, maturity onset diabetes of the young 2)`（MODY2 就是 GCK） | 40 |

### 排除（都实测过，名字带 glucokinase / hexokinase 但不是人 GCK）

| 实体 | resolved_form | 文档数 | 为什么排除 |
|---|---|---:|---|
| `glucokinase regulatory protein` | （未解析） | 302 | **GKRP**，另一类机制 |
| `Glucokinase (hexokinase 4) regulator` | （未解析） | 243 | 同上 |
| `Glucokinase regulator` | `HGNC:4196` | 85 | 同上 |
| `GKRP` | `Q14397` | 162 | 同上 |
| `ADP-dependent glucokinase` 及变体 | 未解析 / `HGNC:25250` | 424 | **ADPGK，另一个酶** |
| `ATP-dependent glucokinase` 及变体 | （未解析） | 27 | 细菌/古菌的酶 |
| `polyphosphate glucokinase` | （未解析） | 126 | 细菌的酶 |
| `glucokinase 1` / `glucokinase-1` | `Q9GTW9` / 未解析 | 35 | 非人物种 |
| `glucokinase-associated dual specificity phosphatase` | `Q9JIM4` | 4 | 另一个蛋白 |
| `glucokinase activity, related sequence 1/2` | （未解析） | 2 | 假基因 |
| 己糖激酶家族其他成员 | `HGNC:4922/4923/4925`（HK1/2/3）、`HGNC:6315`（酮己糖激酶） | ~4,000 | 不是 GCK |
| **泛称 `hexokinase`** | （未解析） | **29,158** | 规模与锚点集相当，但绝大多数是 HK1/HK2（肿瘤代谢），**无法分辨**。单独统计，不进主表 |

**GKRP 相关的另出一张表**，与 ChEMBL 侧把 `CHEMBL3885579`(GCK–GKRP PPI) 单列一致——
GKRP 解离剂与直接激活 GCK 是两类机制，不能无脑合并。

### 不排除但打标

这些 surface form 会带进假阳性，**标记而不删除**，由下游决定：

| surface form | 文档数 | 风险 |
|---|---:|---|
| `GK` | 1,343 | ⚠ 糖尿病文献里 `GK` 更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型），与本领域高度重叠 |
| `4` | 4 | 单个数字被解析成基因，**标注错误** |
| `glk` / `GlkA` / `gukA` / `GlcK` | ~70 | 细菌 glucokinase 基因名 |
| `GKA` / `GKAs` | 144 | 其实是 "glucokinase activator" 的缩写，**方向信号**，正向标记 |

## 二、不做 field 过滤

**这一步不按 `field_id` 筛**。`clms`（权利要求）是精度轴、`desc`（说明书）是召回轴，
在锚定步骤两者都要。但**每篇专利的 4 个部分命中情况记成列**（`hit_ttl` / `hit_abst` /
`hit_desc` / `hit_clms` + 提及次数），下游精筛直接用。

## 三、去重

`family_id` 有两类无效值，必须先排除：`-1` 是「未分配同族」哨兵、`NULL` 是缺失。

```sql
COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)
```

主表**一行一篇专利文档**（保持一行一个实体），`family_id` 作列，同族统计写进报告。

## 四、与 Step1 交叉

每篇专利标注：含多少个 Step1_06 的 782 个候选化合物（全部部分 / 仅权利要求）。
这是连接 ChEMBL 侧与专利侧的唯一桥梁，用 **InChIKey** 对齐。

⚠ SureChEMBL 的 `inchi_key` **不唯一**（30,990,818 行 / 29,874,136 唯一），
join 会产生重复行，要先去重。

## 五、自检：阳性对照

已知 GKA 的专利必须能被捞到。对照集 **8 个**：

| 分子 | ChEMBL | InChIKey |
|---|---|---|
| Ro-28-1675 | `CHEMBL1096435` | 从 Step1_06 读 |
| Piraglitin | `CHEMBL1783734` | 同上 |
| Neriglitin | `CHEMBL2165615` | 同上 |
| PF-04991532 | `CHEMBL2165620` | 同上 |
| AZD-1656 | `CHEMBL3219124` | 同上 |
| MK-0941 | `CHEMBL3580737` | 同上 |
| **BMS-820132** | `CHEMBL5072532` | `OYUDYQMFVRHPIY-UHFFFAOYSA-N` |
| **多格列艾汀 Dorzagliatin** | `CHEMBL4297508` | `HMUMWSORCUWQJO-QAPCUYQASA-N` |

后两个是 Step1 之后补进来的。**多格列艾汀在 ChEMBL 里一条活性记录都没有**
（已在库上核实），所以它不在 Step1 的 782 个候选里——
一个已上市的 GKA 在 ChEMBL 查不到活性，这正是要做专利挖掘的理由。

## 六、必须写进报告的已知偏倚

标注覆盖按专利局严重不均（实测）：

| | 命中文档 | 全库该国专利 |
|---|---:|---:|
| US | 15,899 | 969 万 |
| EP | 7,307 | 527 万 |
| WO | 5,150 | 301 万 |
| **CN** | **849** | **2,388 万** |
| **JP** | **31** | **306 万** |

原因：**JPO 不提供全文**（只有著录项 + 英文标题摘要），**CNIPA 只有英文机翻全文**，
实体标注管道在这两家基本失效。**这批数据实质上是 US / EP / WO 的视图**，
不能据此谈「全球 GKA 专利版图」。

## 输出

- `Step2_01_GCK_Related_Patents.csv` —— 主产物，一行一篇专利文档
- `Step2_01_GCKR_Related_Patents.csv` —— GKRP 相关，单独一张
- `Step2_01_GCK_Related_Patent_Retrieval.md` —— 报告

代码：`Step2_01_GCK_Related_Patent_Retrieval.py`
