# 数据缺陷记录：MAP4K2 被误映射到葡萄糖激酶

**结论先行**：Step1_02 取到的 228 个 assay 中，**52 个（23%）实际测量的是
MAP4K2 / Germinal Center Kinase，而非葡萄糖激酶**。它们在 ChEMBL 37 中被挂到了
`CHEMBL3820`（Hexokinase-4）名下。

这不是分类错误，而是**靶点与 assay 的匹配错误**——比方向判错更根本：
Step1_03 的规则与 LLM 都只回答「这个实验测的是激活还是抑制」，
没有任何环节质疑「这个实验测的到底是不是这个蛋白」。

---

## 1. 问题的来源：GCK 是个歧义缩写

| 缩写 | 含义 | 基因 | UniProt | 长度 | 蛋白类型 |
|---|---|---|---|---|---|
| GCK | **Glucokinase**（葡萄糖激酶） | `GCK` | P35557 | 465 aa | 己糖激酶家族，磷酸化**葡萄糖** |
| GCK | **Germinal Center Kinase** | `MAP4K2` | Q12851 | 820 aa | STE20 家族**蛋白激酶**，磷酸化**蛋白质** |

在激酶研究文献里，"GCK" 默认指 Germinal Center Kinase。激酶选择性谱
（kinase selectivity panel）中出现的 "GCK" 一律是 MAP4K2——
**葡萄糖激酶不是蛋白激酶，不会出现在激酶谱里**。

ChEMBL 中 MAP4K2 本身有独立靶点 `CHEMBL5330`（Q12851），
说明这批数据本应挂在那里。

---

## 2. 调查过程

### 2.1 触发

人工提出质疑：Germinal Center Kinase / MAP4K2 是否被错误映射到了葡萄糖激酶。

此前 Step1_03 的结果中已有可疑迹象，但当时未被识别为靶点问题：
「GCK 抑制」类有 56 个 assay 却只有 204 条活性（每个 assay 平均不到 4 条，
是典型的选择性谱单点数据形态），且大量来自 `DONATED_PROBES` /
`LIT_CHEM_PROBES` 这类化学探针数据源。

### 2.2 决定性证据：探针肽序列比对

有两个 assay 的描述中直接写出了质谱探针肽。ChEMBL 的
`component_sequences.sequence` 存有两个蛋白的完整序列，可直接比对：

```python
pep = "DIKGANLLLTLQGDVK"
# P35557 葡萄糖激酶 (465 aa):  未命中
# Q12851 MAP4K2     (820 aa):  命中，136-151 位
#   上下文 …LKGLHHLHSQGKIHRDIKGANLLLTLQGDVKLADFGVSGELTASVA…
```

`HRD…IKGAN…DFG` 是蛋白激酶催化环的标志性基序。第二个探针肽
`DTVTSELAAVKIVK` 结果相同：MAP4K2 命中，葡萄糖激酶未命中。

全库检索该肽段，只命中人和小鼠的 MAP4K2，无其他蛋白。

### 2.3 蛋白长度矛盾

`CHEMBL3829755` 的描述为 "Inhibition of human GCK (**2 to 812 residues**)"。

葡萄糖激酶全长仅 **465** aa，不可能有第 812 位残基；MAP4K2 全长 **820** aa，
"2 to 812" 恰好是其近全长构建体。

### 2.4 文献语境

52 个可疑 assay 来自 35 篇文献，逐一检查标题，**全部是蛋白激酶论文**：

> LRRK2 · DYRK · EGFR · STK4 · PLK1 · GSK-3β · RIP2 · BTK · CDK · CK2 ·
> TAK1 · RIPK · ROCK · MAP4K4 · JAK/HDAC · TNIK · CHK1 · IGF-1R · hepcidin

其中一篇标题直接写着 "Investigating small molecules to inhibit
**germinal center kinase**-like kinase (GLK/MAP4K3)"。

没有任何一篇是糖尿病 / 代谢 / 葡萄糖激酶激活剂方向的论文。

### 2.5 激酶谱特征

这 35 篇文献在 ChEMBL 中合计涉及 **1,221 个不同靶点、14,816 个 assay**。
共现频次最高的靶点全是蛋白激酶：

| 次数 | 靶点 |
|---:|---|
| 100 | CHK1 |
| 92 | ABL1 |
| 67 | EGFR |
| 62 | 胰岛素受体 |
| 59 | IGF-1R / p38 |
| 58 | Lck / Lyn / RSK3 |

这是典型的激酶选择性谱结构。

### 2.6 化合物零重叠

| 集合 | 数量 |
|---|---:|
| 可疑 assay 涉及的化合物 | 57 |
| 其中也出现在 MAP4K2 (`CHEMBL5330`) assay 中的 | 5 |
| 其中也出现在**真正**葡萄糖激酶 assay 中的 | **0** |

若这批数据真是葡萄糖激酶的，应与 GKA 文献有化合物重叠。实测为零。

---

## 3. 判别规则

五条证据指向同一条极干净的判别规则：

> **描述中使用缩写 "GCK" 而不出现 "glucokinase" 字样的 assay，
> 全部是 MAP4K2。**

52 个全部符合，无例外，也不存在两种写法同时出现的边界情况。
反之，写明 "glucokinase" 的 176 个 assay 无一落入激酶谱文献。

辅助佐证（在 `Step1_03_GCK_Assay_Classification.py` 中一并记录为证据）：

- 残基范围超过 465（葡萄糖激酶全长）
- 描述含 MAP4K2 特征探针肽
- 描述含 "kinase activity assay" 等蛋白激酶专属读数
  （葡萄糖激酶磷酸化的是葡萄糖，不是蛋白）

---

## 4. 影响范围

| 类别 | 修正前 | 修正后 | 说明 |
|---|---:|---:|---|
| GCK 激活 | 142 | **142** | 不受影响 |
| GCK 抑制 | 51 | **7** | 44 个是 MAP4K2 |
| GCK 结合 | 21 | 18 | 3 个是 MAP4K2 |
| GCK–GKRP 相互作用 | 3 | 3 | 不受影响 |
| 细胞或表型效应 | 1 | 1 | 不受影响 |
| 无法判断 | 10 | 5 | 5 个是 MAP4K2 |

**GKA 主线未受污染**：142 个激活 assay、2,861 条活性全部来自写明
"glucokinase" 的描述。被污染的 52 个 assay 仅带 72 条活性，
且高度集中在抑制类——这也解释了为什么「GCK 抑制」类的活性密度异常低。

---

## 5. 处理方式

**标记而不删除**，保留完整证据链以便复核。
`Step1_03_GCK_Assay_Classification.py` 中新增靶点身份校验，输出两列：

- `target_identity_suspect`：`TRUE` / `FALSE`
- `target_identity_evidence`：触发了哪些检查、命中的原文片段

被标记的记录：

- `review_required` 强制置为 `TRUE`
- 不进入第 2 阶段 LLM（其问题在靶点身份，不在方向判定）
- 单独导出到 `Step1_03_target_mismapped.csv`
- 下游筛选 GKA 时应以 `target_identity_suspect == FALSE` 为前置条件

---

## 6. 教训

1. **靶点身份校验必须是独立的一步。** 原有链路
   `accession → target → assay` 隐含假设「挂在该靶点下的 assay 都测的是该蛋白」，
   而这个假设在 ChEMBL 中并不总成立。
2. **基因符号歧义是系统性风险**，不限于 GCK。任何用缩写指代靶点的 assay 描述
   都值得怀疑，尤其当缩写同时是另一个蛋白的常用名时。
3. **数据形态本身是信号。** 「assay 很多但每个只有几条活性」是选择性谱的特征，
   与该靶点的专门研究文献形态截然不同，当时应当引起警觉。
4. **序列是最可靠的裁判。** ChEMBL 自带 `component_sequences.sequence`，
   描述中出现的肽段、残基范围都可以直接比对，不必依赖外部资源。
