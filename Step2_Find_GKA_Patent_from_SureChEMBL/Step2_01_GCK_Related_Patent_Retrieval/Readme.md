从 SureChEMBL 2026-07-17 快照检索**所有与葡萄糖激酶（GCK）相关的专利**，作为 Step2 的锚定表。

**召回优先**：不判断方向（激活/抑制）、不判断是不是 GKA 专利、不做 `field_id` 过滤。
**纯规则检索**，只用 bulk parquet，不调 API、不用 LLM。

## 零、两个原则

### 0.1 SureChEMBL 侧必须独立于 ChEMBL

目标是**两个库各自独立找出 GKA，最后 Union**：

```
ChEMBL 侧    ──独立──→  分子集合 A（Step1，787 个）
SureChEMBL 侧 ──独立──→  专利集合 B → 权要化合物
                              └──→ 最后 Union / 比较
```

**Step1 的候选不进检索式**，只用于两件事：事后测召回、最后比较。
表里 `val_*` 前缀的列都是验证列，**没有参与筛选**。

理由：如果用 ChEMBL 的分子去捞专利，找到的必然是已知化学空间的邻域，
「SureChEMBL 独立贡献了多少」这个问题就变成循环论证了。

> 同理，**结构相似性检索（`fpsim2`）不属于本步骤**——它以 ChEMBL 为种子，
> 是扩展臂不是独立锚点，留到后面单独做。

### 0.2 单一锚点必然漏，所以用并集

上一版只用实体标注，实测漏掉了标题里明写 "GLUCOKINASE ACTIVATOR" 的专利。
原因见下面的坑 ①。**本版用两个独立锚点取并集，再用同族扩展补。**

## 一、检索策略

```
锚点 A  patents.title 正则           ← 独立于标注管道
锚点 B  biomedical_locations 实体     ← resolved_form = 'HGNC:4195' + 白名单
           ↓ 并集
扩展    family_id 展开               ← 同族成员全收，实测 +6,423 篇
           ↓
标注（不筛）  CPC / 命中来源 / field 分布 / 化合物数 / 风险标记
```

### 锚点 A：标题正则

`patents.title` 是 `patents.parquet` 的原生列，44,912,542 行全在本地，
**完全不经过生物医学标注管道**，所以标注缺失的专利照样能命中。

匹配词：`glucokinase`、`glucose kinase`、`hexokinase 4`、`hexokinase IV`（不分大小写）。

实测：标题含 `glucokinase` 的 **1,017 篇 / 266 同族**；
其中 **51 篇是实体法完全捞不到的**。
另外标题含 `glucokinase activator` 的 **631 篇 / 141 同族**——
这个是**方向信号**，标 `title_says_activator`。

### 锚点 B：实体标注

**按 `resolved_form` 锚定，不做文本匹配**——`glucokinase` 在库里有 33 种写法，
含 OCR 破碎形（`glucokmase`、`giucokinase`、`gl uc ok i na se`、`Glucokinas e`），
文本匹配会漏掉一大半。

| 纳入 | 实体数 | 文档数 |
|---|---:|---:|
| `resolved_form = 'HGNC:4195'` | 33 | 29,191 |
| 白名单 `glucokinase (hexokinase 4)`（未解析） | 1 | 102 |
| 白名单 `Glucokinase (hexokinase 4, MODY2)`（未解析，MODY2 就是 GCK） | 1 | 40 |

### 扩展：同族展开

同族成员共享同一份说明书和权利要求，但**标注管道是逐篇文档跑的**——
会出现美国那篇标注全、欧洲那篇标注空。命中任一成员即收全族。

实测：命中 7,618 个同族 → 全部成员 35,749 篇，
其中 **6,423 篇没被任何锚点直接命中**。这一步免费且逻辑无懈可击。

### 排除（都实测过，名字像但不是人 GCK）

| 实体 | resolved_form | 文档数 | 为什么排除 |
|---|---|---:|---|
| `glucokinase regulatory protein` / `Glucokinase (hexokinase 4) regulator` / `GLKRP` | 未解析 | 302 / 243 | **GKRP**，另一类机制 |
| `Glucokinase regulator` / `GKRP` | `HGNC:4196` / `Q14397` | 85 / 162 | 同上 |
| `ADP-dependent glucokinase` 及变体 | 未解析 / `HGNC:25250` | 424 | **ADPGK**，另一个酶 |
| `ATP-dependent glucokinase` 及变体 | 未解析 | 27 | 细菌/古菌的酶 |
| `polyphosphate glucokinase` | 未解析 | 126 | 细菌的酶 |
| `glucokinase 1` / `glucokinase-1` | `Q9GTW9` / 未解析 | 35 | 非人物种 |
| HK1 / HK2 / HK3 / 酮己糖激酶 | `HGNC:4922/4923/4925/6315` | ~4,000 | 不是 GCK |
| **泛称 `hexokinase`** | 未解析 | **29,158** | 规模与锚点集相当，绝大多数是 HK1/HK2（肿瘤代谢），**无法分辨** |

GKRP 相关**另出一张表**，与 ChEMBL 侧单列 `CHEMBL3885579`(GCK–GKRP PPI) 一致。

### 打标不排除

| surface form | 文档数 | 标记 |
|---|---:|---|
| `GK` | 1,343 | ⚠ 糖尿病文献里更常指 **Goto-Kakizaki 大鼠**（2 型糖尿病模型） |
| `4` | 4 | ⚠ 单个数字被解析成基因，**标注错误** |
| `glk` / `GlkA` / `gukA` / `GlcK` | ~70 | 细菌 glucokinase 基因名 |
| `GKA` / `GKAs` | 144 | ✅「glucokinase activator」缩写，**方向正向信号** |

## 二、⚠ 已验证的坑（做之前必读）

### ① 标注管道会整篇缺失，这是上一版最大的错误

**`biomedical_locations` 为空的专利，实体锚定永远找不到，无论标题写得多明白。**

实测四篇（都在权要里含已知 GKA 化合物）：

| 专利 | 标题 | 化合物关联 | 生物医学标注 |
|---|---|---:|---:|
| `EP-4725482-A1` | GLUCOKINASE ACTIVATOR FOR **COGNITIVE DISORDERS AND NEUROLOGICAL**… | 75 条 | **0 条** |
| `US-20260200881-A1` | SULFOXIDE AND SULFONE **GLUCOKINASE ACTIVATORS**… | 238 条 | **0 条** |
| `CN-118453592-A` | **Glucokinase activator** composition for treating diabetes | 151 条 | **0 条** |
| `US-12064416-B2` | Pharmaceutical combination containing **glucose kinase** activator | 158 条 | **0 条** |

化学抽取跑了（几十上百个化合物），生物医学标注一条没产出。
**第一篇正是本项目最该找到的专利——GKA 用于认知障碍与神经系统疾病。**

对策：锚点 A（标题）+ 同族扩展。产物里 `has_biomedical_annotation` 列直接暴露这个缺口。

### ② `field_id` 决定语义，差 6.5 倍

`desc`（说明书）含背景技术，会大段引用**他人的**化合物；
真正主张保护的只在 `clms`（权利要求）。全库 `desc` 12.18 亿关联 vs `clms` 1.87 亿。
实测 glucokinase 在 `desc` 命中 28,272 篇、`clms` 只有 2,597 篇，**差 11 倍**。

本步骤**不过滤**，但四个部分的命中次数分别记成列，下游精筛直接用。

### ③ `family_id` 有哨兵值 `-1`

`-1` 是「未分配同族」的占位（全库 71,862 篇），
直接 `COUNT(DISTINCT family_id)` 会把它们错当成同一个发明。

```sql
COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)
```

另：`family_id` 为空的 1,172,063 篇与 `publication_date` 为空的**完全是同一批**，
是整块元数据缺失，同族扩展对它们无效。

### ④ 「GCK」缩写歧义，与 ChEMBL 那次同源

全大写 `GCK` 在库里 `resolved_form` **是空的**（无法判断是 glucokinase 还是 MAP4K2），
`HPK/GCK-like kinase` → `HGNC:6866` 是 MAP4K 家族。**不能拿 `GCK` 当锚点。**
注意大小写：`GcK`（id 788019）反而正确解析到了 HGNC:4195。

### ⑤ CPC 只能当过滤器，不能当锚点

`A61P 3/10`（抗糖尿病）覆盖了 91% 的 GKA 专利（573/631），召回极好，
但**全库有 204,448 篇**，单独用等于没筛。本步骤只把它记成标注列。

### ⑥ 专利局覆盖严重不均——结论不能外推

| | 命中 | 全库该国 |
|---|---:|---:|
| US / EP / WO | 15,899 / 7,307 / 5,150 | 969 万 / 527 万 / 301 万 |
| **CN** | **849** | **2,388 万** |
| **JP** | **31** | **306 万** |

**JPO 不提供全文**（只有著录项+英文标题摘要），**CNIPA 只有英文机翻全文**。
这批数据实质上是 **US / EP / WO 的视图**。
不能说「中国/日本没有 GKA 专利」——**是看不见，不是没有**。

### ⑦ `inchi_key` 在 SureChEMBL 里不唯一

30,990,818 行 / 29,874,136 唯一，约 110 万重复。join 前必须 `DISTINCT`，
否则行数被放大。（本步骤只在验证列用到。）

### ⑧ 孤儿 patent_id

`biomedical_locations` 里引用了但 `patents` 表中不存在的 `patent_id`，
实测 96 个。join `patents` 后数量会掉，是正常的。

## 三、换靶点怎么复用

靶点相关的东西全收在脚本的 `TARGET` 配置块里，检索逻辑不用动。
换靶点时这几项**必须重新推**，不能照抄：

1. **`resolved_form` 用哪个命名空间**——全库 53% 的实体压根没解析，
   HGNC 前缀 16.7 万、裸 accession 19.9 万、MeSH 13.2 万。先查你的靶点在哪一类。
2. **缩写歧义图谱**——`GK` = Goto-Kakizaki 是 GCK+糖尿病特有的。
   注意实体表里 `WHITE`、`GLASS`、`CLAMP`、`BLIND` 都是真实基因名，
   **如果靶点基因名是常用英文词，实体锚定会淹在噪声里**。
3. **近邻蛋白黑名单**——名字像但不是的那些。
4. **阳性对照**——没有已知临床药就没法标定阈值，这是硬限制。

## 输出

- `Step2_01_GCK_Related_Patents.csv` —— 主产物，一行一篇专利文档
- `Step2_01_GCKR_Related_Patents.csv` —— GKRP 相关，单独一张
- `Step2_01_GCK_Related_Patent_Retrieval.md` —— 报告

代码：`Step2_01_GCK_Related_Patent_Retrieval.py`
