# Step1_03 GCK assay 分类结果（第 1 阶段：规则）

- 输入：`Step1_02_GCK_Assay_Mapping.csv`
- 运行时间：2026-07-30 22:37:56
- assay 总数：**228**，活性合计 **3,262**
- 需要进入第 2 步（LLM）/ 第 4 步（人工）的：**19**

> 本阶段只用规则分类。规则命中不明确的一律标 `review_required = TRUE`，**不强行给标签**。每条判定都附证据原文，可逐条复核。

## ⚠ 靶点身份校验

**52 / 228** 个 assay 的靶点身份存疑——它们挂在葡萄糖激酶 `CHEMBL3820` 下，但实际测量的很可能是 **MAP4K2 / Germinal Center Kinase**（「GCK」是歧义缩写）。

完整调查过程与五条独立证据见 [`Step1_03_Target_Mismapping_MAP4K2.md`](Step1_03_Target_Mismapping_MAP4K2.md)，名单见 `Step1_03_target_mismapped.csv`。

**下游筛选 GKA 时必须加上 `target_identity_suspect == FALSE`。**

## 分类分布

「剔除后」列已排除靶点身份存疑的记录，是实际可用于 GKA 的数量。

| 类别 | assay 数 | 活性数 | 其中需复核 | 剔除误映射后 |
| --- | ---: | ---: | ---: | ---: |
| GCK 激活 | 142 | 2,861 | 6 | **142**（2,861 条活性） |
| GCK 抑制 | 56 | 204 | 53 | **9**（139 条活性） |
| GCK 结合 | 18 | 104 | 1 | **18**（104 条活性） |
| GCK–GKRP 相互作用 | 1 | 40 | 0 | **1**（40 条活性） |
| 无法判断 | 11 | 53 | 11 | **6**（46 条活性） |

## 置信度分布

| 置信度 | assay 数 | 含义 |
| --- | ---: | --- |
| high | 198 | 描述中有明确措辞，且与实测 standard_type 不冲突 |
| medium | 17 | 描述无措辞靠 standard_type 推断，或存在复合/多方向情形 |
| low | 13 | 信号缺失或相互矛盾 |

## 待第 2 步处理的 assay

下列记录规则无法明确判定，已导出到 `Step1_03_pending_llm.csv`。按 README，第 2 步交由 LLM 判断，且**必须输出证据句与置信度**；仍不能确定的进入人工审核。

| assay_chembl_id | 规则暂定 | 置信度 | 规则给出的理由 | 描述 |
| --- | --- | --- | --- | --- |
| `CHEMBL1054393` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Vmax:1） | Activity at human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) |
| `CHEMBL2168355` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Activity:1） | Activity at human recombinant glucokinase assessed as increase in Vmax at 1 uM |
| `CHEMBL2168357` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Activity:1） | Activity at human recombinant glucokinase assessed as decrease in Km at 1 uM |
| `CHEMBL2432300` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Vmax:18） | Activity of human glucokinase at 50 uM |
| `CHEMBL3117036` | GCK 抑制 | medium | 描述命中抑制类措辞：「Inhibition of fluorescein-labeled human GK interaction with biotin-lab…」；描述提及 GKRP 但未构成明确的相互作用表述，需人工确认；standard_type 与描述一致（IC50:30→抑制） | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP incubated for 20 mins pr |
| `CHEMBL3226993` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（EC50:27→激活） | Activity at glucokinase (unknown origin) |
| `CHEMBL3223676` | GCK 抑制 | low | 描述命中抑制类措辞：「Inhibition of recombinant human pancreatic glucokinase by G-6-P dehydr…」；⚠ 描述判定为「GCK 抑制」，但实测 standard_type 指向其他方向（EC50:17→激活） | Inhibition of recombinant human pancreatic glucokinase by G-6-P dehydrogenase assay |
| `CHEMBL3239304` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（S50:18→激活） | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal saturation concentration |
| `CHEMBL3239305` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Vmax:1） | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal activity |
| `CHEMBL3377919` | 无法判断 | low | 描述无方向措辞，standard_type 也全为方向中性（Km:24） | Activity of human recombinant Glucokinase measured over 5 mins by spectrophotometry |
| `CHEMBL3364949` | GCK 抑制 | medium | 描述命中抑制类措辞：「Inhibition of fluorescein-labeled human GK interaction with biotin-lab…」；描述提及 GKRP 但未构成明确的相互作用表述，需人工确认；standard_type 与描述一致（IC50:45→抑制） | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP compound incubated for 2 |
| `CHEMBL3369166` | GCK 抑制 | medium | 描述命中抑制类措辞：「…mouse hepatocytes to IC50 for inhibition of fluorescein-labeled human GK interaction with biotin-lab…」；描述提及 GKRP 但未构成明确的相互作用表述，需人工确认；检测到细胞/体内语境：「…atio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepa…」；同时含方向措辞与表型读数「…atio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepa…」，属复合实验，需人工确认归类 | Ratio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepatocytes to IC50 for inhibition of fl |
| `CHEMBL3751340` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（EC50:1→激活）；描述出现 GKRP，但属于「absence/presence of GKRP」实验条件表述，并非在测 GK–GKRP 相互作用，需人工确认 | Effect on human glucokinase activity after 60 mins by luciferase-based luminescence assay in absence of human  |
| `CHEMBL3888434` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（EC50:38→激活） | Coupled Enzymatic Assay: The assay is carried out according to the protocol outlined in Hariharan et al (1997) |
| `CHEMBL4022295` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（Ratio EC50:1→激活） | Potency index, ratio of (S)-N-(5-chlorothiazol-2-yl)-2-(4-(cyclopropylsulfonyl)-6-methyl-2-oxopyridin-1(2H)-yl |
| `CHEMBL5048925` | GCK 结合 | low | 描述命中结合类措辞：「Displacement of fluorescent labeled derivative from recombinant human hepat…」；⚠ 描述判定为「GCK 结合」，但实测 standard_type 指向其他方向（IC50:2→抑制） | Displacement of fluorescent labeled derivative from recombinant human hepatic glucokinase incubated for 30 min |
| `CHEMBL5261324` | GCK 激活 | medium | 描述无方向措辞，依据实测 standard_type 推断（FC:1→激活） | Agonist activity at glucokinase (unknown origin) assessed as fold increase at 10 uM by luciferase reporter gen |
| `CHEMBL5733424` | GCK 抑制 | medium | 描述命中抑制类措辞：「Inhibition Assay: The in vitro activity of the compounds described her…」；standard_type 同时指向多个方向（IC50:3→抑制、k_off:3→结合、kon:3→结合），已按描述判定 | Inhibition Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK, and other k |
| `CHEMBL5735123` | GCK 抑制 | medium | 描述命中抑制类措辞：「…compounds described herein in inhibiting TAK1, HCK and other kinases were obtained using an Invitrog…」；standard_type 同时指向多个方向（IC50:3→抑制、k_off:3→结合、kon:3→结合），已按描述判定 | In Vitro Activity Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK and o |

## 全部分类明细

| assay_chembl_id | 类别 | 置信度 | 复核 | 活性数 | 描述 |
| --- | --- | --- | --- | ---: | --- |
| `CHEMBL1825590` | GCK 激活 | high | FALSE | 70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5  |
| `CHEMBL1825591` | GCK 激活 | high | FALSE | 70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5  |
| `CHEMBL1825592` | GCK 激活 | high | FALSE | 70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 20 |
| `CHEMBL1825593` | GCK 激活 | high | FALSE | 70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5  |
| `CHEMBL2353001` | GCK 激活 | high | FALSE | 51 | Activation of glucokinase (unknown origin) |
| `CHEMBL3768920` | GCK 激活 | high | FALSE | 50 | Activation of recombinant human glucokinase assessed as NADPH formation using glucose as substrate incubated f |
| `CHEMBL2353000` | GCK 激活 | high | FALSE | 48 | Activation of glucokinase (unknown origin) relative to control |
| `CHEMBL4193347` | GCK 激活 | high | FALSE | 47 | Activation of recombinant human liver glucokinase 2 assessed as reduction in NADH production in presence of 5  |
| `CHEMBL4193348` | GCK 激活 | high | FALSE | 47 | Activation of recombinant human liver glucokinase 2 assessed as reduction in NADH production in presence of 5  |
| `CHEMBL3131381` | GCK 激活 | high | FALSE | 44 | Activation of glucokinase (unknown origin) using glucose as substrate at 10 uM relative to control |
| `CHEMBL3768919` | GCK 激活 | high | FALSE | 43 | Activation of recombinant human glucokinase assessed as NADPH formation using glucose as substrate incubated f |
| `CHEMBL3889107` | GCK 激活 | high | FALSE | 42 | Glucokinase-Activating Assay: To test the exemplified compounds, the following assay was employed. Recombinant |
| `CHEMBL3226995` | GCK 激活 | high | FALSE | 38 | Activation of glucokinase (unknown origin) |
| `CHEMBL3888434` | GCK 激活 | medium | TRUE | 38 | Coupled Enzymatic Assay: The assay is carried out according to the protocol outlined in Hariharan et al (1997) |
| `CHEMBL1032034` | GCK 激活 | high | FALSE | 34 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosp |
| `CHEMBL1032035` | GCK 激活 | high | FALSE | 34 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosp |
| `CHEMBL1042963` | GCK 激活 | high | FALSE | 34 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) at 10 uM by G6PDH-coupled spectrometry |
| `CHEMBL4614223` | GCK 激活 | high | FALSE | 34 | Activation of glucokinase (unknown origin) |
| `CHEMBL959301` | GCK 激活 | high | FALSE | 33 | Activation of human glucokinase by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric ass |
| `CHEMBL959302` | GCK 激活 | high | FALSE | 33 | Activation of human glucokinase by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric ass |
| `CHEMBL959303` | GCK 激活 | high | FALSE | 33 | Activation of human glucokinase assessed as maximal response by glucose-6-phosphate dehydrogenase coupled cont |
| `CHEMBL959304` | GCK 激活 | high | FALSE | 33 | Activation of human glucokinase assessed as maximal response by glucose-6-phosphate dehydrogenase coupled cont |
| `CHEMBL1788085` | GCK 激活 | high | FALSE | 33 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometry |
| `CHEMBL1788086` | GCK 激活 | high | FALSE | 33 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometry relative to control |
| `CHEMBL2167270` | GCK 激活 | high | FALSE | 32 | Activation of human recombinant glucokinase assessed as concentration required for 1.5 fold increase in enzyma |
| `CHEMBL3223666` | GCK 激活 | high | FALSE | 32 | Activation of recombinant human glucokinase assessed as formation of glucose-6-phosphate by G6PDH/NADP coupled |
| `CHEMBL1053631` | GCK 激活 | high | FALSE | 31 | Activation of His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) assessed |
| `CHEMBL1260806` | GCK 激活 | high | FALSE | 31 | Activation of glucokinase assessed as concentration required to 50% increase in enzyme activity |
| `CHEMBL1056161` | GCK 激活 | high | FALSE | 30 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE |
| `CHEMBL1056163` | GCK 激活 | high | FALSE | 30 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE |
| `CHEMBL1056162` | GCK 激活 | high | FALSE | 29 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE |
| `CHEMBL1056164` | GCK 激活 | high | FALSE | 29 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE |
| `CHEMBL1167927` | GCK 激活 | high | FALSE | 29 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH by spectrometry |
| `CHEMBL2025594` | GCK 激活 | high | FALSE | 29 | Activation of human recombinant glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH assesse |
| `CHEMBL2410527` | GCK 激活 | high | FALSE | 29 | Activation of human recombinant glucokinase by matrix assay in presence of glucose |
| `CHEMBL4272644` | GCK 激活 | high | FALSE | 29 | Activation of human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate at 10 uM preincub |
| `CHEMBL4272645` | GCK 激活 | high | FALSE | 29 | Activation of human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate at 10 uM preincub |
| `CHEMBL2318649` | GCK 激活 | high | FALSE | 28 | Induction of human pancreatic glucokinase activity at 5 mM glucose concentration |
| `CHEMBL2410525` | GCK 激活 | high | FALSE | 28 | Activation of human recombinant glucokinase assessed as ratio of the enzyme velocity at maximum compound conce |
| `CHEMBL2410526` | GCK 激活 | high | FALSE | 28 | Activation of human recombinant glucokinase assessed as ratio of the glucokinase Km at maximum compound concen |
| `CHEMBL1167928` | GCK 激活 | high | FALSE | 27 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH by spectrometry rel |
| `CHEMBL2025595` | GCK 激活 | high | FALSE | 27 | Activation of human recombinant glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH assesse |
| `CHEMBL3226993` | GCK 激活 | medium | TRUE | 27 | Activity at glucokinase (unknown origin) |
| `CHEMBL4427549` | GCK 激活 | high | FALSE | 27 | Activation of human recombinant glucokinase using 5 mM glucose monitored over 5 mins in presence of NAD+ by gl |
| `CHEMBL3384773` | GCK 激活 | high | FALSE | 26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) using 5 mM glucose by spectrophotometry in p |
| `CHEMBL3384774` | GCK 激活 | high | FALSE | 26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) using 5 mM glucose by spectrophotometry in p |
| `CHEMBL3384775` | GCK 激活 | high | FALSE | 26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) assessed as Vmax by spectrophotometry in pre |
| `CHEMBL4022283` | GCK 激活 | high | FALSE | 26 | Activation of full length human C-terminal FLAG-tagged glucokinase (12 to 465 residues) expressed in Escherich |
| `CHEMBL4022294` | GCK 激活 | high | FALSE | 26 | Activation of full length human C-terminal FLAG-tagged glucokinase (12 to 465 residues) expressed in Escherich |
| `CHEMBL4427552` | GCK 激活 | high | FALSE | 26 | Activation of human recombinant glucokinase assessed as substrate Vmax ratio monitored over 5 mins in presence |
| `CHEMBL5048861` | GCK 激活 | high | FALSE | 26 | Activation of recombinant full length human liver glucokinase expressed in Escherichia coli incubated for 10 m |
| `CHEMBL5048862` | GCK 激活 | high | FALSE | 26 | Activation of recombinant full length human liver glucokinase expressed in Escherichia coli assessed as maxima |
| `CHEMBL957641` | GCK 激活 | high | FALSE | 25 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosp |
| `CHEMBL957643` | GCK 激活 | high | FALSE | 25 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosp |
| `CHEMBL1058313` | GCK 激活 | high | FALSE | 25 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucos |
| `CHEMBL1058314` | GCK 激活 | high | FALSE | 25 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucos |
| `CHEMBL3377918` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant Glucokinase measured over 5 mins by G6-PD coupled assay in presence of 5 mM gl |
| `CHEMBL4427348` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant glucokinase using 5 mM glucose as substrate in presence of NAD+ by glucose 6-p |
| `CHEMBL4427349` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant glucokinase using 5 mM glucose as substrate in presence of NAD+ and 4% HSA by  |
| `CHEMBL4427350` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant glucokinase assessed as enzyme affinity for glucose using 5 mM glucose as subs |
| `CHEMBL4427351` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant glucokinase assessed as Vmax ratio using 0.16 to 80 mM of glucose as substrate |
| `CHEMBL4427551` | GCK 激活 | high | FALSE | 25 | Activation of human recombinant glucokinase assessed as substrate affinity monitored over 5 mins in presence o |
| `CHEMBL957642` | GCK 激活 | high | FALSE | 24 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli assessed as maxima |
| `CHEMBL957569` | GCK 激活 | high | FALSE | 24 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli assessed as maxima |
| `CHEMBL3377920` | GCK 激活 | high | FALSE | 24 | Activity of activated human recombinant Glucokinase assessed as maximum glucose phosphorylation measured over  |
| `CHEMBL3578684` | GCK 激活 | high | FALSE | 24 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate by G6PDH coupled assay |
| `CHEMBL3578770` | GCK 激活 | high | FALSE | 24 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate assessed as enzyme hal |
| `CHEMBL840244` | GCK 激活 | high | FALSE | 23 | Potency in Glucokinase activation assay |
| `CHEMBL1781588` | GCK 激活 | high | FALSE | 23 | Activation of His-tagged recombinant glucokinase expressed in Escherichia coli using [14C]-glucose substrate b |
| `CHEMBL3225185` | GCK 激活 | high | FALSE | 23 | Activation of recombinant human His6-tagged glucokinase by G6PDH-coupled spectrophotometry |
| `CHEMBL3225186` | GCK 激活 | high | FALSE | 22 | Activation of recombinant human His6-tagged glucokinase by G6PDH-coupled spectrophotometry relative to control |
| `CHEMBL2346010` | GCK 激活 | high | FALSE | 21 | Activation of glucokinase (unknown origin) using glucose as substrate assessed as stimulation concentration re |
| `CHEMBL3243402` | GCK 激活 | high | FALSE | 20 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose by spectrophotometry |
| `CHEMBL4427550` | GCK 激活 | high | FALSE | 20 | Activation of human recombinant glucokinase using 5 mM glucose monitored over 5 mins in presence of NAD+ and g |
| `CHEMBL2215199` | GCK 激活 | high | FALSE | 19 | Activation of human glucokinase |
| `CHEMBL2215687` | GCK 激活 | high | FALSE | 19 | Ratio of recombinant human glucokinase Km at maximum activator concentration to recombinant human glucokinase  |
| `CHEMBL2216127` | GCK 激活 | high | FALSE | 19 | Ratio of recombinant human glucokinase velocity at maximum activator concentration to recombinant human glucok |
| `CHEMBL2216128` | GCK 激活 | high | FALSE | 19 | Activation of recombinant human glucokinase assessed measuring rate of glucose 6-phosphate formation using G6P |
| `CHEMBL2432299` | GCK 激活 | high | FALSE | 18 | Activation of human glucokinase assessed as glucose concentration at enzyme's half-maximal phosphorylation rat |
| `CHEMBL2432301` | GCK 激活 | high | FALSE | 18 | Allosteric activation of human glucokinase using glucose as substrate measured every 10 secs for 5 mins |
| `CHEMBL3239304` | GCK 激活 | medium | TRUE | 18 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal saturation concentration |
| `CHEMBL3583541` | GCK 激活 | high | FALSE | 17 | Activation of recombinant human glucokinase assessed as formation of glucose-6-phosphate by G6PDH/NADP coupled |
| `CHEMBL5261510` | GCK 激活 | high | FALSE | 17 | Activation of recombinant human glucokinase assessed as increase in maximal activation by measuring accumulati |
| `CHEMBL1065927` | GCK 激活 | high | FALSE | 16 | Activation of human liver glucokinase expressed in CHO cells at 2.5 mM glucose concentration by glucose-6-phos |
| `CHEMBL1038417` | GCK 激活 | high | FALSE | 16 | Activation of human liver glucokinase expressed in CHO cells at 10 mM glucose concentration by glucose-6-phosp |
| `CHEMBL1038416` | GCK 激活 | high | FALSE | 14 | Activation of human liver glucokinase expressed in CHO cells at 2.5 mM glucose concentration by glucose-6-phos |
| `CHEMBL1038418` | GCK 激活 | high | FALSE | 14 | Activation of human liver glucokinase expressed in CHO cells at 10 mM glucose concentration by glucose-6-phosp |
| `CHEMBL4014754` | GCK 激活 | high | FALSE | 14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL4014755` | GCK 激活 | high | FALSE | 14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL4014756` | GCK 激活 | high | FALSE | 14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL4014757` | GCK 激活 | high | FALSE | 14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL990114` | GCK 激活 | high | FALSE | 13 | Activation of GST fused human liver glucokinase assessed as generation of NADPH by G6PDH coupled assay |
| `CHEMBL2168593` | GCK 激活 | high | FALSE | 13 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometric analysis |
| `CHEMBL990115` | GCK 激活 | high | FALSE | 12 | Activation of GST fused human liver glucokinase assessed as generation of NADPH by G6PDH coupled assay relativ |
| `CHEMBL990116` | GCK 激活 | high | FALSE | 12 | Activation of GST fused human liver glucokinase assessed as concentration required to double glucokinase activ |
| `CHEMBL1035543` | GCK 激活 | high | FALSE | 11 | Activation of glucokinase |
| `CHEMBL4624346` | GCK 激活 | high | FALSE | 9 | Activation of human glucokinase |
| `CHEMBL4014791` | GCK 激活 | high | FALSE | 7 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL4014792` | GCK 激活 | high | FALSE | 7 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in pr |
| `CHEMBL4614239` | GCK 激活 | high | FALSE | 7 | Activation of glucokinase (unknown origin) in presence of 4% human serum albumin |
| `CHEMBL865109` | GCK 激活 | high | FALSE | 6 | Activation of glucokinase |
| `CHEMBL1058318` | GCK 激活 | high | FALSE | 6 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucos |
| `CHEMBL1042962` | GCK 激活 | high | FALSE | 6 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) by G6PDH-coupled spectrometry |
| `CHEMBL3095515` | GCK 激活 | high | FALSE | 6 | Activation of recombinant human glucokinase by G6PDH/NADP coupled assay |
| `CHEMBL1825594` | GCK 激活 | high | FALSE | 5 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5  |
| `CHEMBL1825595` | GCK 激活 | high | FALSE | 5 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 20 |
| `CHEMBL3095513` | GCK 激活 | high | FALSE | 5 | Activation of recombinant human glucokinase assessed as increase in enzyme Vmax for glucose by G6PDH/NADP coup |
| `CHEMBL3095514` | GCK 激活 | high | FALSE | 5 | Activation of recombinant human glucokinase assessed as decrease in enzyme Km for glucose by G6PDH/NADP couple |
| `CHEMBL3224178` | GCK 激活 | high | FALSE | 5 | Activation of glucokinase (unknown origin) |
| `CHEMBL3224190` | GCK 激活 | high | FALSE | 5 | Activation of glucokinase (unknown origin) |
| `CHEMBL4427376` | GCK 激活 | high | FALSE | 5 | Ratio of EC50 for activation of human recombinant glucokinase in presence of 4% HSA to EC50 for activation of  |
| `CHEMBL4428666` | GCK 激活 | high | FALSE | 5 | Activation of human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21 |
| `CHEMBL4193353` | GCK 激活 | high | FALSE | 4 | Activation of recombinant human liver glucokinase 2 assessed as assessed as reduction in Km for glucose in pre |
| `CHEMBL2412836` | GCK 激活 | high | FALSE | 3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 |
| `CHEMBL2412837` | GCK 激活 | high | FALSE | 3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 |
| `CHEMBL2412840` | GCK 激活 | high | FALSE | 3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 |
| `CHEMBL2412842` | GCK 激活 | high | FALSE | 3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 |
| `CHEMBL4349275` | GCK 激活 | high | FALSE | 3 | Activation of glucokinase (unknown origin) at 10 uM relative to control |
| `CHEMBL4427586` | GCK 激活 | high | FALSE | 3 | Ratio of EC50 for activation of human recombinant glucokinase in presence of 4% human serum albumin to EC50 fo |
| `CHEMBL3131525` | GCK 激活 | high | FALSE | 2 | Activation of human glucokinase |
| `CHEMBL3131380` | GCK 激活 | high | FALSE | 2 | Activation of glucokinase (unknown origin) using glucose as substrate |
| `CHEMBL5048855` | GCK 激活 | high | FALSE | 2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli incubated for 10 m |
| `CHEMBL5048926` | GCK 激活 | high | FALSE | 2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli assessed as maxima |
| `CHEMBL5048927` | GCK 激活 | high | FALSE | 2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli incubated for 10 m |
| `CHEMBL5048928` | GCK 激活 | high | FALSE | 2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli assessed as maxima |
| `CHEMBL5048931` | GCK 激活 | high | FALSE | 2 | Activation of human glucokinase assessed as decrease in glucose Km |
| `CHEMBL5048932` | GCK 激活 | high | FALSE | 2 | Activation of human glucokinase assessed as increase rate of glucose turnover |
| `CHEMBL1054394` | GCK 激活 | high | FALSE | 1 | Activation of His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) assessed |
| `CHEMBL2025596` | GCK 激活 | high | FALSE | 1 | Ratio of EC50 for GK50 to compound for activation of human recombinant glucokinase expressed in Escherichia co |
| `CHEMBL2167261` | GCK 激活 | high | FALSE | 1 | Activation of human recombinant glucokinase |
| `CHEMBL2318453` | GCK 激活 | high | FALSE | 1 | Induction of human hepatic glucokinase activity at 5 mM glucose concentration |
| `CHEMBL2318648` | GCK 激活 | high | FALSE | 1 | Induction of human pancreatic glucokinase activity assessed as maximal reaction rate at 5 mM glucose concentra |
| `CHEMBL2318650` | GCK 激活 | high | FALSE | 1 | Induction of human pancreatic glucokinase activity assessed as half maximal saturation concentration of glucos |
| `CHEMBL2346006` | GCK 激活 | high | FALSE | 1 | Activation of human recombinant glucokinase using glucose as substrate |
| `CHEMBL2412841` | GCK 激活 | high | FALSE | 1 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 |
| `CHEMBL3377935` | GCK 激活 | high | FALSE | 1 | Activation of human recombinant Glucokinase measured over 5 mins by G6-PD coupled assay in presence of 5 mM gl |
| `CHEMBL3578772` | GCK 激活 | high | FALSE | 1 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate assessed as enzyme Vma |
| `CHEMBL3751340` | GCK 激活 | medium | TRUE | 1 | Effect on human glucokinase activity after 60 mins by luciferase-based luminescence assay in absence of human  |
| `CHEMBL4022295` | GCK 激活 | medium | TRUE | 1 | Potency index, ratio of (S)-N-(5-chlorothiazol-2-yl)-2-(4-(cyclopropylsulfonyl)-6-methyl-2-oxopyridin-1(2H)-yl |
| `CHEMBL4349277` | GCK 激活 | high | FALSE | 1 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate dehyd |
| `CHEMBL5261324` | GCK 激活 | medium | TRUE | 1 | Agonist activity at glucokinase (unknown origin) assessed as fold increase at 10 uM by luciferase reporter gen |
| `CHEMBL5261513` | GCK 激活 | high | FALSE | 1 | Activation of recombinant human glucokinase assessed as increase in maximal activation by measuring accumulati |
| `CHEMBL3364949` | GCK 抑制 | medium | TRUE | 45 | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP compound incubated for 2 |
| `CHEMBL3117036` | GCK 抑制 | medium | TRUE | 30 | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP incubated for 20 mins pr |
| `CHEMBL4813785` | GCK 抑制 | high | FALSE | 19 | Inhibition of human hexokinase-4 at 20 uM relative to control |
| `CHEMBL3223676` | GCK 抑制 | low | TRUE | 17 | Inhibition of recombinant human pancreatic glucokinase by G-6-P dehydrogenase assay |
| `CHEMBL5733424` | GCK 抑制 | medium | TRUE | 9 | Inhibition Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK, and other k |
| `CHEMBL5735123` | GCK 抑制 | medium | TRUE | 9 | In Vitro Activity Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK and o |
| `CHEMBL4332276` | GCK 抑制 | high | FALSE | 7 | Inhibition of recombinant human N-terminal GST-fused glucokinase expressed in Escherichia coli using glucose-6 |
| `CHEMBL4804384` | GCK 抑制 | high | TRUE | 6 | Inhibition of GCK in human NCI-H929 cells by mass spectroscopic analysis |
| `CHEMBL4804566` | GCK 抑制 | high | TRUE | 5 | Inhibition of GCK in human NCI-H929 cells at 10 uM by mass spectroscopic analysis relative to control |
| `CHEMBL1942767` | GCK 抑制 | high | TRUE | 4 | Inhibition of human GCK in HL-60 cells lysate assessed as reduction of labeling of acyl-phosphate ATP probe at |
| `CHEMBL3630507` | GCK 抑制 | high | TRUE | 3 | Inhibition of GCK (unknown origin) at 10 uM after 120 mins P33 radiolabeled kinase activity assay |
| `CHEMBL4118352` | GCK 抑制 | high | TRUE | 3 | Inhibition of GCK in human SKCO1 cells at 1 uM after 4 hrs using biotin labeled DIKGANLLLTLQGDVK probe by mass |
| `CHEMBL3369166` | GCK 抑制 | medium | TRUE | 2 | Ratio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepatocytes to IC50 for inhibition of fl |
| `CHEMBL4327702` | GCK 抑制 | high | TRUE | 2 | Inhibition of human GCK assessed as residual activity at 1 uM using MBP as substrate by [gamma-33P]-ATP assay  |
| `CHEMBL4674372` | GCK 抑制 | high | TRUE | 2 | Inhibition of human GCK using MBP as substrate assessed as residual activity at 1 uM by [gamma-33P]-ATP assay  |
| `CHEMBL3750455` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL3807041` | GCK 抑制 | high | FALSE | 1 | Inhibition of HK4 (unknown origin) |
| `CHEMBL3829755` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK (2 to 812 residues) assessed as remaining enzyme activity at 50 uM after 30 mins by 33 |
| `CHEMBL4034328` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK Lysine 1 labelling site (unknown origin) at 10 uM |
| `CHEMBL4034329` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK Lysine 2 labelling site (unknown origin) at 10 uM |
| `CHEMBL4045684` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 1 uM |
| `CHEMBL4045685` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 1 uM |
| `CHEMBL4045926` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 0.1 uM |
| `CHEMBL4045927` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 0.1 uM |
| `CHEMBL4050296` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK assessed as residual activity at 1 uM in presence of 33P-ATP by filter-binding assay r |
| `CHEMBL4057298` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4057494` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 0.49 uM relative to control |
| `CHEMBL4120341` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4120568` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 0.19 uM relative to control |
| `CHEMBL4120792` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 0.164 uM relative to control |
| `CHEMBL4328066` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK assessed as residual activity at 100 uM using MBP as substrate by [gamma-33P]-ATP assa |
| `CHEMBL4328430` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4352383` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK at 10 uM using MBP as substrate in presence of [gamma-33P]-ATP |
| `CHEMBL4379338` | GCK 抑制 | high | TRUE | 1 | Inhibition of recombinant human GCK (1 to 473 residues) using myelin basic protein as substrate after 40 mins  |
| `CHEMBL4674071` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4674112` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK using MBP as substrate at 0.5 uM by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4681971` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK human PBMC lysates at 1 uM by ActivX screen assay relative to control |
| `CHEMBL4720720` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4721072` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK assessed as residual activity using MBP as substrate at 1 uM by [gamma-33P]-ATP assay  |
| `CHEMBL4880638` | GCK 抑制 | medium | TRUE | 1 | GCK (h) Millipore kinase activity assay |
| `CHEMBL4884756` | GCK 抑制 | medium | TRUE | 1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885047` | GCK 抑制 | medium | TRUE | 1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885337` | GCK 抑制 | medium | TRUE | 1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4887878` | GCK 抑制 | medium | TRUE | 1 | GCK(h) Millipore kinase panel |
| `CHEMBL5038899` | GCK 抑制 | high | TRUE | 1 | Inhibition of human recombinant GCK assessed as reduction in substrate phosphorylation at 50 to 500 nM using A |
| `CHEMBL5058199` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL5059136` | GCK 抑制 | medium | TRUE | 1 | GCK(h) Kinase panel |
| `CHEMBL5156973` | GCK 抑制 | high | TRUE | 1 | Inhibition of human recombinant GCK (1 to 473 residues) using myelin basic protein as substrate in presence of |
| `CHEMBL5464035` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464534` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464912` | GCK 抑制 | high | TRUE | 1 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5657730` | GCK 抑制 | high | TRUE | 1 | Inhibition of human GCK assessed as remaining activity at 10 uM in presence of ATP by radiometric kinase assay |
| `CHEMBL5679451` | GCK 抑制 | high | TRUE | 1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5679684` | GCK 抑制 | high | TRUE | 1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5681798` | GCK 抑制 | high | TRUE | 1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5682058` | GCK 抑制 | high | TRUE | 1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL3582523` | GCK 结合 | high | FALSE | 17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as disso |
| `CHEMBL3583538` | GCK 结合 | high | FALSE | 17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as on ra |
| `CHEMBL3583539` | GCK 结合 | high | FALSE | 17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as off r |
| `CHEMBL3583540` | GCK 结合 | high | FALSE | 17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as disso |
| `CHEMBL4428660` | GCK 结合 | high | FALSE | 6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escheric |
| `CHEMBL4428665` | GCK 结合 | high | FALSE | 6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escheric |
| `CHEMBL4428667` | GCK 结合 | high | FALSE | 6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escheric |
| `CHEMBL4428668` | GCK 结合 | high | FALSE | 6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escheric |
| `CHEMBL3751341` | GCK 结合 | high | FALSE | 2 | Binding affinity to human glucokinase by surface plasmon resonance analysis |
| `CHEMBL5048925` | GCK 结合 | low | TRUE | 2 | Displacement of fluorescent labeled derivative from recombinant human hepatic glucokinase incubated for 30 min |
| `CHEMBL2412838` | GCK 结合 | high | FALSE | 1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 at 1 |
| `CHEMBL2412839` | GCK 结合 | high | FALSE | 1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 at 1 |
| `CHEMBL2412843` | GCK 结合 | high | FALSE | 1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 by i |
| `CHEMBL3119606` | GCK 结合 | high | FALSE | 1 | Binding affinity to biotinylated GK (unknown origin) by SPR analysis |
| `CHEMBL4428661` | GCK 结合 | high | FALSE | 1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia col |
| `CHEMBL4428662` | GCK 结合 | high | FALSE | 1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia col |
| `CHEMBL4428663` | GCK 结合 | high | FALSE | 1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia col |
| `CHEMBL4428664` | GCK 结合 | high | FALSE | 1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia col |
| `CHEMBL3579280` | GCK–GKRP 相互作用 | high | FALSE | 40 | Inhibition of human biotin-labeled GKRP and fluorescein-labeled human GK interaction preincubated for 20 mins  |
| `CHEMBL3377919` | 无法判断 | low | TRUE | 24 | Activity of human recombinant Glucokinase measured over 5 mins by spectrophotometry |
| `CHEMBL2432300` | 无法判断 | low | TRUE | 18 | Activity of human glucokinase at 50 uM |
| `CHEMBL5464286` | 无法判断 | low | TRUE | 2 | % Activity remaining of GCK in the Dundee kinase panel at 1.0 µM |
| `CHEMBL6194772` | 无法判断 | low | TRUE | 2 | Effect of GCK(h) at compound concentration of 1.0 uM using the Cerep Kinase panel |
| `CHEMBL1054393` | 无法判断 | low | TRUE | 1 | Activity at human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) |
| `CHEMBL2168355` | 无法判断 | low | TRUE | 1 | Activity at human recombinant glucokinase assessed as increase in Vmax at 1 uM |
| `CHEMBL2168357` | 无法判断 | low | TRUE | 1 | Activity at human recombinant glucokinase assessed as decrease in Km at 1 uM |
| `CHEMBL3239305` | 无法判断 | low | TRUE | 1 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal activity |
| `CHEMBL4881060` | 无法判断 | low | TRUE | 1 | GCK(h) Eurofins Kinase panel |
| `CHEMBL4884319` | 无法判断 | low | TRUE | 1 | GCK(h) Eurofins kinase panel |
| `CHEMBL5724151` | 无法判断 | low | TRUE | 1 | activity of GCK(h) at 1.0 µM in the Eurofins Kinase panel |

