# Step2_01 GCK 相关专利检索

- 数据源：**SureChEMBL 2.0，`2026-07-17` 全量快照**
- 快照路径：`/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17`
- 运行时间：2026-07-31 20:02:46
- 靶点：**GCK / glucokinase**（`HGNC:4195`）
- 命中：**35,793** 篇专利文档 → 同族去重 **7,620** 个发明

> **召回优先**：不判断方向、不判断是不是 GKA 专利、**不做 `field_id` 过滤**。
> **独立于 ChEMBL**：Step1 候选不进检索式，`val_*` 列仅作验证。

## 一、两个独立锚点 + 同族扩展

| 命中来源 | 专利数 | 同族 | 说明 |
| --- | ---: | ---: | --- |
| `entity` | 28,340 | 7,407 | 只有实体标注命中（标题没写靶点名，常见于化合物专利） |
| `family` | 6,422 | 2,927 | 本身没命中，靠同族成员带进来 |
| `title+entity` | 975 | 260 | 两个锚点都命中，最可靠 |
| `title` | 56 | 37 | **只有标题命中——实体标注没抓到** |

锚点 A（标题正则）命中 1,031 篇，锚点 B（实体）命中 29,411 篇，并集 29,467 篇；同族展开后 **35,889** 篇（补 6,422 篇）。

### ⚠ 标注管道整篇缺失的证据

本表 **1,335** 篇专利的 `biomedical_locations` **一条记录都没有**（`has_biomedical_annotation = FALSE`），其中 **49** 篇是标题正则捞回来的——**只用实体锚定这些会被整个漏掉**。

标题明写 activator、却完全没有生物医学标注的（按公开日倒序）：

| 专利 | 国 | 公开日 | 标题 | 化合物数 |
| --- | --- | --- | --- | ---: |
| `US-20260200881-A1` | US | 2026-07-16 | SULFOXIDE AND SULFONE GLUCOKINASE ACTIVATORS AND METHODS OF US | 236 |
| `EP-4714442-A3` | EP | 2026-06-03 | PRODRUG OF PYRROLIDONE DERIVATIVES AS GLUCOKINASE ACTIVATOR | 1 |
| `EP-4725482-A1` | EP | 2026-04-15 | GLUCOKINASE ACTIVATOR FOR COGNITIVE DISORDERS AND NEURODEGENER | 72 |
| `EP-4714442-A2` | EP | 2026-03-25 | PRODRUG OF PYRROLIDONE DERIVATIVES AS GLUCOKINASE ACTIVATOR | 156 |
| `EP-4682144-A1` | EP | 2026-01-21 | SOLID FORM OF PYRROLIDONE DERIVATIVE AS GLUCOKINASE ACTIVATOR | 108 |
| `EP-3804714-B1` | EP | 2024-12-25 | PHARMACEUTICAL COMBINATION AND COMPOSITION, AND COMBINATION PR | 148 |
| `EP-3804716-B1` | EP | 2024-12-25 | PHARMACEUTICAL COMBINATION, COMPOSITION, AND COMBINATION PREPA | 146 |
| `EP-3804715-B1` | EP | 2024-12-25 | PHARMACEUTICAL COMBINATION, COMPOSITION AND COMPOUND PREPARATI | 153 |

这些专利的**化学侧是完整的**（几十上百个化合物），只是文本标注为空。

## 二、方向信号（不做判定，只标记）

| 信号 | 专利数 | 同族 | 说明 |
| --- | ---: | ---: | --- |
| 标题写着 activator | 638 | 143 | `title_says_activator`，**最强的方向信号** |
| 出现 `GKA` / `GKAs` 缩写 | 139 | — | 缩写本身就是 glucokinase activator |

> bulk 数据判不了方向（无全文，`Mechanism` 实体全是工业化学词）。这两列是**规则能拿到的全部方向信息**，真正的方向判定要读权利要求原文。

## 三、命中位置分布（不过滤，只记录）

- `ttl` 标题：1,066 篇（3.0%）
- `abst` 摘要：1,435 篇（4.0%）
- `desc` 说明书：28,287 篇（79.0%）
- **`clms` 权利要求**：2,587 篇（7.2%）

`clms` 是精度轴、`desc` 是召回轴，锚定步骤两个都留。全库这两者差 6.5 倍（12.18 亿 vs 1.87 亿关联）。

## 四、专利局分布 ⚠

| 专利局 | 命中 | 同族 | 全库该国 | 命中率 |
| --- | ---: | ---: | ---: | ---: |
| `US` | 16,806 | 5,911 | 9,691,977 | 0.2% |
| `EP` | 8,598 | 4,368 | 5,265,735 | 0.2% |
| `WO` | 6,341 | 5,328 | 3,008,081 | 0.2% |
| `CN` | 3,618 | 2,113 | 23,884,165 | 0.0% |
| `JP` | 430 | 363 | 3,062,582 | 0.0% |

**JPO 不提供全文**（只有著录项 + 英文标题摘要），**CNIPA 只有英文机翻全文**——标注管道在这两家基本失效。**这批数据实质上是 US / EP / WO 的视图**，不能说「中国/日本没有 GKA 专利」，是看不见不是没有。

## 五、风险标记

| surface form | 专利数 | 备注 |
| --- | ---: | --- |
| `glucokinase` | 29,096 |  |
| `GK` | 1,338 | ⚠ 糖尿病文献里更常指 Goto-Kakizaki 大鼠（2 型糖尿病模型），与本领域高度重叠 |
| `GKAs` | 131 | ✅「glucokinase activators」的缩写，**方向正向信号** |
| `Hexokinase 4` | 103 |  |
| `glucokinase (hexokinase 4)` | 102 |  |
| `glucokmase` | 80 |  |
| `giucokinase` | 51 |  |
| `GlcK` | 48 | 细菌 glucokinase 基因名 |
| `Glucokinase (hexokinase 4, maturity onset diabetes of the young 2)` | 40 |  |
| `glu- cokinase` | 37 |  |
| `Hexokinase-4` | 25 |  |
| `gluco- kinase` | 20 |  |
| `GcK` | 11 |  |
| `GKA` | 10 | ✅「glucokinase activator」的缩写，**方向正向信号** |
| `gluco kinase` | 10 |  |

- 带风险标记 **1,420** 篇
- **只靠风险形命中、标题也没写靶点名** 的 **10** 篇 ← 最可疑，下游优先人工核

另：`A61P 3/10`（抗糖尿病用途）命中 11,321 篇（31.6%）。CPC 召回好但全库 20 万篇，**只能当过滤器不能当锚点**，本表只作标注。

## 六、验证：与 ChEMBL 侧的召回对照

> 以下全部是**事后验证**，Step1 的分子**没有参与检索**。

| 口径 | 全库有多少 | 本表捞到 | 召回率 |
| --- | ---: | ---: | ---: |
| 含 Step1 候选化合物（任意部分） | 1,096 | 955 | 87.1% |
| **权利要求里含 Step1 候选** | 285 | 271 | **95.1%** |

第二行是更扎实的口径——权要里主张已知 GKA 化合物的专利，本表覆盖了多少。

阳性对照（12 个，来自 Step1 整合表）：

| 对照 | 结构在 SureChEMBL | 命中专利 |
| --- | :---: | ---: |
| AZD-1656 | ✅ | 274 |
| Cadisegliatin | ✅ | 267 |
| Dorzagliatin | ✅ | 159 |
| Piraglitin | ✅ | 158 |
| PF-04991532 | ✅ | 49 |
| Neriglitin | ✅ | 43 |
| Ro-28-1675 (参比化合物) | ✅ | 33 |
| MK-0941 (free base) | ✅ | 24 |
| LY-2608204 / Globalagliatin | ✅ | 11 |
| BMS-820132 | ✅ | 7 |
| MK-0941 (mesylate) | ✅ | 1 |
| Globalagliatin HCl | **❌ 无对应结构** | — |

**11/12 个对照捞到了专利。** 「无对应结构」的是盐型——SureChEMBL 里不单独注册盐，母体条目已命中，药物层没真丢。

## 七、GKRP 单独成表

**651** 篇（216 个同族），与主表重叠 **373** 篇（主表 `sibling_gckr` 列标出）。

GKRP / GCKR 解离剂与直接激活 GCK 是两类机制，与 ChEMBL 侧单列 `CHEMBL3885579`(GCK–GKRP PPI) 的处理一致。

## 八、这一步没做什么

| 没做 | 为什么 | 留给谁 |
| --- | --- | --- |
| 方向判定 | bulk 无全文；`Mechanism` 实体全是工业化学词、`resolved_form` 全空 | 后续读权利要求原文 |
| `field_id` 过滤 | 锚定要召回 | 下游用 `hit_*` 列 |
| 结构相似性检索 | **以 ChEMBL 为种子，不独立**，属扩展臂 | 单独一步 |
| 泛称 `hexokinase`（29,158 篇） | 无法分辨 HK1/2/3 与 GCK | 需结构或全文佐证 |

