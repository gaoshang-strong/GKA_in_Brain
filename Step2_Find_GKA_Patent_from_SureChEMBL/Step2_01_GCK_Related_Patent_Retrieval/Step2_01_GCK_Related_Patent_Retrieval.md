# Step2_01 GCK 相关专利检索

- 数据源：**SureChEMBL 2.0，`2026-07-17` 全量快照**
- 快照路径：`/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17`
- 运行时间：2026-07-31 18:23:00
- 交叉引用：`Step1_GKA_Candidates_with_Properties.csv`（Step1_06 的 782 个候选）
- 命中专利文档：**29,315** 篇 → 同族去重后 **7,608** 个发明

> **召回优先**：本步骤不判断方向（激活/抑制）、不判断是不是 GKA 专利、**不做 `field_id` 过滤**。纯规则检索，未调 API、未用 LLM。

## 一、锚定规则与实测规模

锚点实体 **35** 个。核心是 `resolved_form = 'HGNC:4195'`，**不做文本匹配**——专利全文来自 OCR，`glucokinase` 有几十种破碎写法。

| 纳入方式 | 实体数 | 说明 |
| --- | ---: | --- |
| `resolved` | 33 | `resolved_form = 'HGNC:4195'`，覆盖全部 OCR 变体与缩写 |
| `whitelist` | 2 | 未解析但明确是 GCK 的写法，手工白名单 |

白名单逐条（未解析，人工确认是 GCK）：

- `glucokinase (hexokinase 4)`
- `Glucokinase (hexokinase 4, maturity onset diabetes of the young 2)`　← MODY2 就是 GCK

### 明确排除的（名字带 glucokinase / hexokinase 但不是人 GCK）

| 类别 | 判据 | 涉及专利 | 排除理由 |
| --- | --- | ---: | --- |
| ADP-dependent glucokinase 及变体 | `HGNC:25250` | 996 | ADPGK，另一个酶 |
| ATP-dependent glucokinase 及变体 | `ATP-dependent glucokinase`, `ATP dependent glucokinase`… | 27 | 细菌/古菌的酶 |
| polyphosphate glucokinase | `polyphosphate glucokinase` | 126 | 细菌的酶 |
| glucokinase 1 / -1 | `Q9GTW9` | 35 | 非人物种 |
| glucokinase-associated dual specificity phosphatase | `Q9JIM4` | 4 | 另一个蛋白 |
| glucokinase 假基因 | `glucokinase activity, related sequence 1`, `glucokinase activity, related sequence 2`… | 2 | 假基因/类似序列 |
| 己糖激酶家族其他成员 | `HGNC:4922`, `HGNC:4923`, `HGNC:4925`, `HGNC:6315` | 6,276 | HK1 / HK2 / HK3 / 酮己糖激酶，不是 GCK |
| **泛称 `hexokinase`** | 未解析 | **29,158** | 规模与锚点集相当，但绝大多数是 HK1/HK2（肿瘤代谢），**无法分辨**，不进主表 |

## 二、规模与去重

**29,315 篇专利文档 → 7,608 个同族**（平均每族 3.85 篇），另有 38 篇没有有效 `family_id`（哨兵 `-1` 或 `NULL`）。

> 去重必须 `COUNT(DISTINCT family_id) FILTER (WHERE family_id > 0)`——`-1` 是「未分配同族」的哨兵，不排除会把它们错当成同一个发明。

## 三、⚠ 专利局覆盖严重不均（必读）

| 专利局 | 命中文档 | 同族 | 全库该国专利 | 命中率 |
| --- | ---: | ---: | ---: | ---: |
| `US` | 15,940 | 5,758 | 9,691,977 | 0.2% |
| `EP` | 7,326 | 4,149 | 5,265,735 | 0.1% |
| `WO` | 5,169 | 4,835 | 3,008,081 | 0.2% |
| `CN` | 849 | 596 | 23,884,165 | 0.0% |
| `JP` | 31 | 27 | 3,062,582 | 0.0% |

**JP 与 CN 的命中率低到不能用**。原因不在检索式，在数据源本身：

- **JPO 不提供全文**，SureChEMBL 只拿到著录项 + 英文标题摘要
- **CNIPA 只有英文机器翻译全文**，实体标注管道在机翻文本上效果差

> **这批数据实质上是 US / EP / WO 的视图。**不能据此谈「全球 GKA 专利版图」，也不能说「中国/日本没有 GKA 专利」——是看不见，不是没有。

## 四、命中位置分布

| 位置 | 命中专利数 | 占比 |
| --- | ---: | ---: |
| `ttl` 标题 | 1,066 | 3.6% |
| `abst` 摘要 | 1,435 | 4.9% |
| `desc` 说明书 | 28,287 | 96.5% |
| **`clms` 权利要求** | 2,587 | 8.8% |

本步骤**不按位置过滤**，四个位置的命中次数都写进主表的 `hit_ttl` / `hit_abst` / `hit_desc` / `hit_clms` 列，下游精筛直接用。`clms` 是精度轴、`desc` 是召回轴，别在锚定步骤就二选一。

## 五、surface form 分解与风险标记

| surface form | 命中专利数 | 备注 |
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
| `gluco-kinase` | 8 |  |
| `GlkA` | 8 | 细菌 glucokinase 基因名 |
| `glk` | 7 | 细菌 glucokinase 基因名 |
| `Glucoki- nase` | 7 |  |
| `gki` | 6 | 缩写，需核 |

- 带风险标记的专利 **1,420** 篇（`risk_flags` 列）
- **只靠风险形命中、没有任何可靠写法** 的 **10** 篇 ← **这批最可疑，下游应优先人工核或直接排除**
- 出现 `GKA` / `GKAs` 缩写的 **139** 篇（`has_activator_abbrev = TRUE`）——这两个缩写本身就是「glucokinase activator」，是**方向正向信号**，虽然 SureChEMBL 把它解析成了基因

## 六、与 Step1（ChEMBL 侧）的交叉

用 **InChIKey** 对齐（唯一的跨库桥梁）。417 个 Step1 候选在 SureChEMBL 里找到了对应结构。

| 项目 | 专利数 |
| --- | ---: |
| 含至少 1 个 Step1 候选化合物 | 778 |
| 权利要求里含 Step1 候选化合物 | 254 |

> ⚠ SureChEMBL 的 `inchi_key` 不唯一（30,990,818 行 / 29,874,136 唯一），join 前必须 `DISTINCT`，否则行数会被放大。本脚本已处理。

## 七、阳性对照自检

已知 GKA 的专利必须能被这套规则捞到。落不进说明规则有问题，不是数据有问题。

| 对照 | 结构在 SureChEMBL | 命中专利数 | 其中权要也提到 GCK |
| --- | :---: | ---: | ---: |
| Cadisegliatin | ✅ | 180 | 44 |
| AZD-1656 | ✅ | 171 | 24 |
| Dorzagliatin | ✅ | 114 | 54 |
| Piraglitin | ✅ | 100 | 20 |
| PF-04991532 | ✅ | 37 | 11 |
| Ro-28-1675 (参比化合物) | ✅ | 30 | 14 |
| Neriglitin | ✅ | 30 | 11 |
| MK-0941 (free base) | ✅ | 19 | 13 |
| LY-2608204 / Globalagliatin | ✅ | 9 | 6 |
| BMS-820132 | ✅ | 7 | 4 |
| MK-0941 (mesylate) | ✅ | 1 | 0 |
| Globalagliatin HCl | **❌ 无对应结构** | — | — |

**自检结论：11/12 个对照捞到了专利。**

其中 **1** 个对照的 InChIKey 在 SureChEMBL 里**没有对应结构**，从没进入匹配池——这不是检索式的问题：

- **Globalagliatin HCl**

**盐型通常不会在 SureChEMBL 里单独注册**（专利写的是游离碱结构）。这正是 CLAUDE.md 那条的实证：**跨库对齐结构要先归到母体再取 InChIKey**，拿盐型的 key 去对必然对不上。母体条目已经命中，所以药物层面没有真丢。

## 八、GKRP 单独成表

**651** 篇（216 个同族），与主表重叠 **350** 篇。

GKRP（`HGNC:4196` / `Q14397` / `glucokinase regulatory protein` 等）**不并入主表**——GKRP 解离剂与直接激活 GCK 是两类机制，与 ChEMBL 侧把 `CHEMBL3885579`(GCK–GKRP PPI) 单列的处理一致。重叠的那部分两张表都有，按需取用。

## 九、这一步没做什么

| 没做 | 为什么 | 留给谁 |
| --- | --- | --- |
| 方向判定（激活 / 抑制） | bulk 数据无全文，`Mechanism` 实体类型是工业化学词，`resolved_form` 全空，判不了 | 后续步骤（读权利要求原文） |
| `field_id` 过滤 | 锚定步骤要召回，`clms` 与 `desc` 都留着 | 下游精筛，用主表的 hit_* 列 |
| 结构相似性检索 | 本步骤只做实体锚定 | 后续用 `fpsim2_fingerprints.h5` |
| 泛称 `hexokinase` 的那 29,158 篇 | 无法分辨 HK1/2/3 与 GCK | 若要补召回，需结构或全文佐证 |

