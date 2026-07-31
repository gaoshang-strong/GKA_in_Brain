# Step1_02 GCK assay 清单

- ChEMBL 版本：**CHEMBL_37（2026-05-01）**
- 数据库文件：`/ShangGaoAIProjects/GKA_in_Brain/ChEMBL/ChEMBL_37/chembl_37/chembl_37_sqlite/chembl_37.db`
- 运行时间：2026-07-30 21:55:11
- 输入 tid：20095, 117123
- assay 总数：**228**，其上活性数据点合计 **3,262**

> 本步骤只清点、不判定方向。下面的分布统计与完整清单用于人工判断哪些 assay 真正用于测量 GKA。

## 按靶点

| target | pref_name | target_type | assay 数 | 活性数 |
| --- | --- | --- | ---: | ---: |
| `CHEMBL3820` | Hexokinase-4 | SINGLE PROTEIN | 227 | 3,222 |
| `CHEMBL3885579` | Glucokinase/Glucokinase regulatory protein | PROTEIN-PROTEIN INTERACTION | 1 | 40 |

## 实验类型 assay_type

> 注意：绝大多数 GCK assay 被归为 B(Binding)，**该字段在这里几乎没有区分力**，不能用它判断实验测的是激活还是抑制。

| 取值 | assay 数 | 活性数 |
| --- | ---: | ---: |
| B | 226 | 3,237 |
| F | 1 | 6 |
| A | 1 | 19 |

## 实验格式 bao_format

> BAO 本体对实验形式的描述，比 assay_type 细一些。

| 取值 | assay 数 | 活性数 |
| --- | ---: | ---: |
| single protein format | 145 | 1,879 |
| assay format | 60 | 1,166 |
| cell-based format | 13 | 84 |
| tissue-based format | 8 | 92 |
| subcellular format | 1 | 1 |
| protein format | 1 | 40 |

## 物种 assay_organism

| 取值 | assay 数 | 活性数 |
| --- | ---: | ---: |
| Homo sapiens | 225 | 3,197 |
| (空) | 3 | 65 |

## 靶点可信度 confidence_score

> 官方含义：9 = 直接指认到单一蛋白，8 = 同源单一蛋白，5 = 可能对应多个蛋白。它描述的是靶点归因的粒度与方式，**不是数据质量分**；每行对应的官方原文见 CSV 的 `confidence_desc` 列。

| 取值 | assay 数 | 活性数 |
| --- | ---: | ---: |
| 9 | 222 | 3,077 |
| 8 | 5 | 145 |
| 5 | 1 | 40 |

## 数据来源 src_short_name

| 取值 | assay 数 | 活性数 |
| --- | ---: | ---: |
| LITERATURE | 203 | 3,132 |
| DONATED_PROBES | 14 | 16 |
| LIT_CHEM_PROBES | 5 | 5 |
| BINDINGDB | 4 | 98 |
| PATENT | 2 | 11 |

## 全部 assay 上实测的 standard_type

> 这是**事实统计**，不是方向判定。但它是判断某个 assay 是否在测 GKA 的最强线索：`EC50`/`Emax`/`%max`/`S0.5`/`Vmax` 指向激活与酶动力学表征，`IC50`/`Inhibition` 指向抑制，`Kd` 指向结合。

| standard_type | 活性数 |
| --- | ---: |
| `EC50` | 1,378 |
| `Ratio` | 451 |
| `FC` | 330 |
| `Activity` | 264 |
| `Emax` | 230 |
| `IC50` | 139 |
| `%max` | 109 |
| `Inhibition` | 59 |
| `S0.5` | 50 |
| `Vmax` | 48 |
| `K` | 48 |
| `Kd` | 29 |
| `Km` | 26 |
| `max activation` | 22 |
| `S50` | 18 |
| `T1/2` | 17 |
| `Ratio EC50` | 10 |
| `TIME` | 7 |
| `% inhibition` | 6 |
| `k_off` | 6 |
| `kon` | 6 |
| `% residual kinase activity` | 2 |
| `Kcat` | 2 |
| `% Activity remaining` | 2 |
| `Effect` | 2 |
| `activity` | 1 |

## 突变体 / 实验参数 / 实验分类的覆盖情况

| 字段 | 有值的 assay 数 | 说明 |
| --- | ---: | --- |
| `variant_id` | 0 / 228 | 非空表示该实验用的是突变体蛋白（如激酶耐药突变）；GCK 的天然激活突变若被单独建为 variant 会出现在这里 |
| `assay_parameters` | 5 / 228 | 实验参数，已聚合为 JSON 列表 |
| `assay_classifications` | 0 / 228 | 实验分类，已聚合为 JSON 列表；该表主要覆盖体内/治疗领域实验 |

实际出现的参数类型：`EFO_ID`（17）、`EFO_TERM`（17）

> 注意：这里的参数**不是**葡萄糖浓度、孵育时间这类实验条件，而是 EFO/MONDO 疾病本体标注。GCK 这批 assay 的实验条件（葡萄糖浓度等）只写在 `assay_description` 自由文本里，没有被结构化到 `assay_parameters`。要按葡萄糖浓度分层，必须自己解析描述文本。

## 文献年份分布

| 年份 | assay 数 |
| ---: | ---: |
| 2005 | 1 |
| 2006 | 1 |
| 2008 | 3 |
| 2009 | 31 |
| 2010 | 3 |
| 2011 | 12 |
| 2012 | 16 |
| 2013 | 27 |
| 2014 | 22 |
| 2015 | 7 |
| 2016 | 27 |
| 2017 | 28 |
| 2018 | 9 |
| 2019 | 8 |
| 2020 | 6 |
| 2021 | 9 |
| 2022 | 12 |
| 2023 | 5 |
| 2025 | 1 |

## 完整 assay 清单

按活性数据量降序。完整字段见同目录 CSV。

| assay_chembl_id | type | conf | 活性数 | 化合物数 | standard_types | 描述 |
| --- | --- | ---: | ---: | ---: | --- | --- |
| `CHEMBL1825590` | B | 9 | 70 | 70 | Ratio:70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5 mM glucose measured as ratio of enzyme activation in treated to untreated control by G6PDH coupled continuous spectrophotometric assay |
| `CHEMBL1825591` | B | 9 | 70 | 70 | Ratio:70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5 mM glucose measured as ratio of enzyme activation in compound treated to Ro-28-1675 by G6PDH coupled continuous spectrophotometric assay |
| `CHEMBL1825592` | B | 9 | 70 | 70 | Ratio:70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 20 mM glucose measured as ratio of enzyme activation in treated to untreated control by G6PDH coupled continuous spectrophotometric assay |
| `CHEMBL1825593` | B | 9 | 70 | 70 | Ratio:70 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5 mM glucose measured as ratio of enzyme activation in compound treated to Ro-28-1675 by G6PDH coupled continuous spectrophotometric assay |
| `CHEMBL2353001` | B | 9 | 51 | 51 | EC50:51 | Activation of glucokinase (unknown origin) |
| `CHEMBL3768920` | B | 9 | 50 | 44 | Emax:50 | Activation of recombinant human glucokinase assessed as NADPH formation using glucose as substrate incubated for 30 mins in presence of NADP+ and glucose 6-phosphate dehydrogenase relative to control |
| `CHEMBL2353000` | B | 9 | 48 | 48 | %max:48 | Activation of glucokinase (unknown origin) relative to control |
| `CHEMBL4193347` | B | 9 | 47 | 47 | EC50:47 | Activation of recombinant human liver glucokinase 2 assessed as reduction in NADH production in presence of 5 mM glucose by G6PDH coupled assay |
| `CHEMBL4193348` | B | 9 | 47 | 47 | FC:47 | Activation of recombinant human liver glucokinase 2 assessed as reduction in NADH production in presence of 5 mM glucose by G6PDH coupled assay relative to control |
| `CHEMBL3364949` | B | 9 | 45 | 45 | IC50:45 | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP compound incubated for 20 mins prior to addition of fluorescein-labeled GK measured after 2 to 4 hrs by AlphaScreen assay |
| `CHEMBL3131381` | B | 9 | 44 | 44 | FC:44 | Activation of glucokinase (unknown origin) using glucose as substrate at 10 uM relative to control |
| `CHEMBL3768919` | B | 9 | 43 | 37 | EC50:43 | Activation of recombinant human glucokinase assessed as NADPH formation using glucose as substrate incubated for 30 mins in presence of NADP+ and glucose 6-phosphate dehydrogenase |
| `CHEMBL3889107` | B | 8 | 42 | 34 | EC50:42 | Glucokinase-Activating Assay: To test the exemplified compounds, the following assay was employed. Recombinant human liver glucokinase was expressed as a FLAG fusion protein in E. coli, and purified on ANTIFLAG M2 AFFINITY GEL (Sigma). The assay was carried out at 30° C. in a 96-well plate. In the plate was distributed 69 ul each of assay buffer (25 mM Hepes Buffer: pH=7.2, 2 mM MgCl2, 1 mM ATP, 0.5 mM TNAD, 1 mM dithiothreitol), to which was added 1 ul of a DMSO solution of the compound or DMSO as control. Then, 20 ul of pre-ice-cooled enzyme mixture (FLAG-GK, 20 U/ml G6PDH) was distributed thereto, to which was added 10 ul of 25 mM glucose as substrate to initiate the reaction (final glucose concentration=2.5 mM). After starting the reaction, the absorbance at 405 nm was measured every 30 seconds for 10 minutes to evaluate the compound based on the initial increase for 5 minutes. FLAG-GK was added so that the increase of absorbance after 5 minutes fell between 0.05 to 0.1 in the presence of 1% DMSO. The OD values of the respective compounds were measured in the respective concentrations, wherein the OD value of DMSO as control is regarded as 100%. From the OD values at the respective concentrations, Emax (%, 2.5 mM Glu) and EC50 (nM, 2.5 mM Glu) were calculated and used as indicators of the GK activation capability of the compounds. According to the above assay, the GK activation capability of the exemplified compounds of the present invention was determined. The following table shows the results. Where different enantiomers of the same compound were tested, two numbers are provided. Where 3 numbers are provided, for example, Ex. #30, data for the racemic form is also shown. |
| `CHEMBL3579280` | B | 5 | 40 | 40 | IC50:40 | Inhibition of human biotin-labeled GKRP and fluorescein-labeled human GK interaction preincubated for 20 mins prior to fluorescein-labeled human GK addition measured after 2 to 4 hrs by AlphaScreen assay |
| `CHEMBL3226995` | B | 9 | 38 | 36 | EC50:38 | Activation of glucokinase (unknown origin) |
| `CHEMBL3888434` | B | 8 | 38 | 34 | EC50:38 | Coupled Enzymatic Assay: The assay is carried out according to the protocol outlined in Hariharan et al (1997), Diabetes 46: 11-16. Briefly, the test compounds are incubated in a reaction mix containing 25 mM HEPES (pH 7.2), 10 mM MgCl2, 100 mM KCl, 5 mM ATP, 2 mM DTT, 0.5 mM NAD, 1 U/mL Leuconostoc mesenteroides G6PD, 0.3 U/mL of purified human recombinant GK, and different concentrations of glucose. Enzymatic activity is calculated from the initial reaction velocity, measured from the change in NADH absorbance as a function of time. |
| `CHEMBL1032034` | B | 9 | 34 | 32 | EC50:34 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose |
| `CHEMBL1032035` | B | 9 | 34 | 32 | EC50:34 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose |
| `CHEMBL1042963` | B | 9 | 34 | 34 | FC:34 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) at 10 uM by G6PDH-coupled spectrometry relative to untreated control |
| `CHEMBL4614223` | B | 9 | 34 | 34 | EC50:34 | Activation of glucokinase (unknown origin) |
| `CHEMBL959301` | B | 9 | 33 | 33 | EC50:33 | Activation of human glucokinase by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose |
| `CHEMBL959302` | B | 9 | 33 | 33 | EC50:33 | Activation of human glucokinase by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose |
| `CHEMBL959303` | B | 9 | 33 | 33 | FC:33 | Activation of human glucokinase assessed as maximal response by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose relative to control |
| `CHEMBL959304` | B | 9 | 33 | 33 | FC:33 | Activation of human glucokinase assessed as maximal response by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose relative to control |
| `CHEMBL1788085` | B | 9 | 33 | 33 | EC50:33 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometry |
| `CHEMBL1788086` | B | 9 | 33 | 33 | Emax:33 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometry relative to control |
| `CHEMBL2167270` | B | 9 | 32 | 32 | Activity:32 | Activation of human recombinant glucokinase assessed as concentration required for 1.5 fold increase in enzymatic activity |
| `CHEMBL3223666` | B | 9 | 32 | 32 | EC50:32 | Activation of recombinant human glucokinase assessed as formation of glucose-6-phosphate by G6PDH/NADP coupled assay |
| `CHEMBL1053631` | B | 9 | 31 | 31 | Activity:31 | Activation of His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) assessed as drug level required for half-maximal activation |
| `CHEMBL1260806` | B | 8 | 31 | 31 | Activity:31 | Activation of glucokinase assessed as concentration required to 50% increase in enzyme activity |
| `CHEMBL1056161` | B | 9 | 30 | 30 | EC50:30 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose |
| `CHEMBL1056163` | B | 9 | 30 | 30 | EC50:30 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose |
| `CHEMBL3117036` | B | 9 | 30 | 30 | IC50:30 | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP incubated for 20 mins prior to addition of fluorescein-labeled GK measured after 2 to 4 hrs by Alpha Screen assay |
| `CHEMBL1056162` | B | 9 | 29 | 29 | Emax:29 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide |
| `CHEMBL1056164` | B | 9 | 29 | 29 | Emax:29 | Activation of N-terminal His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide |
| `CHEMBL1167927` | B | 9 | 29 | 29 | EC50:29 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH by spectrometry |
| `CHEMBL2025594` | B | 9 | 29 | 29 | EC50:29 | Activation of human recombinant glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH assessed as glucose 6-phosphate formation by spectrometric analysis |
| `CHEMBL2410527` | B | 9 | 29 | 29 | EC50:29 | Activation of human recombinant glucokinase by matrix assay in presence of glucose |
| `CHEMBL4272644` | B | 9 | 29 | 29 | FC:29 | Activation of human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate at 10 uM preincubated for 10 mins followed by 5 mM glucose addition measured for 5 mins at 30 secs time interval by G6PDH/NADP coupled assay relative to control |
| `CHEMBL4272645` | B | 9 | 29 | 29 | FC:29 | Activation of human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate at 10 uM preincubated for 10 mins followed by 5 mM glucose addition measured for 5 mins at 30 secs time interval by G6PDH/NADP coupled assay relative to RO0281675 |
| `CHEMBL2318649` | B | 9 | 28 | 28 | EC50:28 | Induction of human pancreatic glucokinase activity at 5 mM glucose concentration |
| `CHEMBL2410525` | B | 9 | 28 | 28 | Ratio:28 | Activation of human recombinant glucokinase assessed as ratio of the enzyme velocity at maximum compound concentration to the enzyme Km in the absence of compound |
| `CHEMBL2410526` | B | 9 | 28 | 28 | Ratio:28 | Activation of human recombinant glucokinase assessed as ratio of the glucokinase Km at maximum compound concentration to the enzyme Km in the absence of compound |
| `CHEMBL1167928` | B | 9 | 27 | 27 | FC:27 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH by spectrometry relative to control |
| `CHEMBL2025595` | B | 9 | 27 | 27 | FC:27 | Activation of human recombinant glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH assessed as glucose 6-phosphate formation by spectrometric analysis relative to control |
| `CHEMBL3226993` | B | 9 | 27 | 27 | EC50:27 | Activity at glucokinase (unknown origin) |
| `CHEMBL4427549` | B | 9 | 27 | 27 | EC50:27 | Activation of human recombinant glucokinase using 5 mM glucose monitored over 5 mins in presence of NAD+ by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL3384773` | B | 9 | 26 | 26 | EC50:26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) using 5 mM glucose by spectrophotometry in presence of NAD+ and glucose 6-phosphate dehydrogenase |
| `CHEMBL3384774` | B | 9 | 26 | 26 | EC50:26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) using 5 mM glucose by spectrophotometry in presence of NAD+ and glucose 6-phosphate dehydrogenase in presence of 4% HSA |
| `CHEMBL3384775` | B | 9 | 26 | 26 | Vmax:26 | Activation of purified human glucokinase isoform 3 (13 to 466 aa) assessed as Vmax by spectrophotometry in presence of varying concentration of glucose relative to control |
| `CHEMBL4022283` | B | 9 | 26 | 26 | EC50:26 | Activation of full length human C-terminal FLAG-tagged glucokinase (12 to 465 residues) expressed in Escherichia coli DH10b using glucose as substrate after 60 mins in presence of ATP by kinase-glo luminescence assay |
| `CHEMBL4022294` | B | 9 | 26 | 26 | %max:26 | Activation of full length human C-terminal FLAG-tagged glucokinase (12 to 465 residues) expressed in Escherichia coli DH10b using glucose as substrate after 60 mins in presence of ATP by kinase-glo luminescence assay relative to control |
| `CHEMBL4427552` | B | 9 | 26 | 26 | Ratio:26 | Activation of human recombinant glucokinase assessed as substrate Vmax ratio monitored over 5 mins in presence of varying concentration of glucose and NAD+ by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL5048861` | B | 9 | 26 | 26 | EC50:26 | Activation of recombinant full length human liver glucokinase expressed in Escherichia coli incubated for 10 mins by plate reader analysis |
| `CHEMBL5048862` | B | 9 | 26 | 26 | Emax:26 | Activation of recombinant full length human liver glucokinase expressed in Escherichia coli assessed as maximal activation incubated for 10 mins by plate reader analysis relative to control |
| `CHEMBL957641` | B | 9 | 25 | 25 | EC50:25 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 2.5 mM glucose |
| `CHEMBL957643` | B | 9 | 25 | 25 | EC50:25 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay in presence of 10 mM glucose |
| `CHEMBL1058313` | B | 9 | 25 | 24 | EC50:25 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucose-6-phosphate dehydrogenase by spectrophotometry |
| `CHEMBL1058314` | B | 9 | 25 | 24 | Emax:25 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucose-6-phosphate dehydrogenase activity by spectrophotometry relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide |
| `CHEMBL3377918` | B | 9 | 25 | 25 | EC50:25 | Activation of human recombinant Glucokinase measured over 5 mins by G6-PD coupled assay in presence of 5 mM glucose |
| `CHEMBL4427348` | B | 9 | 25 | 25 | EC50:25 | Activation of human recombinant glucokinase using 5 mM glucose as substrate in presence of NAD+ by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL4427349` | B | 9 | 25 | 25 | EC50:25 | Activation of human recombinant glucokinase using 5 mM glucose as substrate in presence of NAD+ and 4% HSA by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL4427350` | B | 9 | 25 | 25 | S0.5:25 | Activation of human recombinant glucokinase assessed as enzyme affinity for glucose using 5 mM glucose as substrate in presence of NAD+ by by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL4427351` | B | 9 | 25 | 25 | Ratio:25 | Activation of human recombinant glucokinase assessed as Vmax ratio using 0.16 to 80 mM of glucose as substrate in presence of NAD+ by glucose 6-phosphate dehydrogenase coupled assay relative to control |
| `CHEMBL4427551` | B | 9 | 25 | 25 | S0.5:25 | Activation of human recombinant glucokinase assessed as substrate affinity monitored over 5 mins in presence of varying concentration of glucose and NAD+ by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL957642` | B | 9 | 24 | 24 | Activity:24 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli assessed as maximal activating response by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide in presence of 2.5 mM glucose |
| `CHEMBL957569` | B | 9 | 24 | 24 | Activity:24 | Activation of flag-tagged human recombinant liver glucokinase expressed in Escherichia coli assessed as maximal activating response by glucose-6-phosphate dehydrogenase coupled continuous spectrophotometric assay relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide in presence of 10 mM glucose |
| `CHEMBL3377919` | B | 9 | 24 | 24 | Km:24 | Activity of human recombinant Glucokinase measured over 5 mins by spectrophotometry |
| `CHEMBL3377920` | B | 9 | 24 | 24 | Ratio:24 | Activity of activated human recombinant Glucokinase assessed as maximum glucose phosphorylation measured over 5 mins by spectrophotometry relative to unactivated enzyme |
| `CHEMBL3578684` | B | 9 | 24 | 24 | EC50:24 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate by G6PDH coupled assay |
| `CHEMBL3578770` | B | 9 | 24 | 24 | Activity:24 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate assessed as enzyme half maximal saturation concentration at 1 uM by G6PDH coupled assay |
| `CHEMBL840244` | B | 9 | 23 | 23 | EC50:23 | Potency in Glucokinase activation assay |
| `CHEMBL1781588` | B | 8 | 23 | 23 | EC50:23 | Activation of His-tagged recombinant glucokinase expressed in Escherichia coli using [14C]-glucose substrate by spectrophotometrically |
| `CHEMBL3225185` | B | 9 | 23 | 23 | EC50:23 | Activation of recombinant human His6-tagged glucokinase by G6PDH-coupled spectrophotometry |
| `CHEMBL3225186` | B | 9 | 22 | 22 | max activation:22 | Activation of recombinant human His6-tagged glucokinase by G6PDH-coupled spectrophotometry relative to control |
| `CHEMBL2346010` | B | 9 | 21 | 21 | Activity:21 | Activation of glucokinase (unknown origin) using glucose as substrate assessed as stimulation concentration required to achieve 150% activity relative to control |
| `CHEMBL3243402` | B | 9 | 20 | 20 | EC50:20 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose by spectrophotometry |
| `CHEMBL4427550` | B | 9 | 20 | 20 | EC50:20 | Activation of human recombinant glucokinase using 5 mM glucose monitored over 5 mins in presence of NAD+ and glucose 6-phosphate 4% human serum albumin by glucose 6-phosphate dehydrogenase coupled assay |
| `CHEMBL2215199` | B | 9 | 19 | 19 | EC50:19 | Activation of human glucokinase |
| `CHEMBL2215687` | B | 9 | 19 | 19 | Ratio:19 | Ratio of recombinant human glucokinase Km at maximum activator concentration to recombinant human glucokinase Km in absence of activator |
| `CHEMBL2216127` | B | 9 | 19 | 19 | Ratio:19 | Ratio of recombinant human glucokinase velocity at maximum activator concentration to recombinant human glucokinase velocity in absence of activator |
| `CHEMBL2216128` | B | 9 | 19 | 19 | EC50:19 | Activation of recombinant human glucokinase assessed measuring rate of glucose 6-phosphate formation using G6PDH/NADP coupling |
| `CHEMBL4813785` | A | 9 | 19 | 19 | Inhibition:19 | Inhibition of human hexokinase-4 at 20 uM relative to control |
| `CHEMBL2432299` | B | 9 | 18 | 18 | Activity:18 | Activation of human glucokinase assessed as glucose concentration at enzyme's half-maximal phosphorylation rate at 50 uM |
| `CHEMBL2432300` | B | 9 | 18 | 18 | Vmax:18 | Activity of human glucokinase at 50 uM |
| `CHEMBL2432301` | B | 9 | 18 | 18 | EC50:18 | Allosteric activation of human glucokinase using glucose as substrate measured every 10 secs for 5 mins |
| `CHEMBL3239304` | B | 9 | 18 | 18 | S50:18 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal saturation concentration |
| `CHEMBL3223676` | B | 9 | 17 | 17 | EC50:17 | Inhibition of recombinant human pancreatic glucokinase by G-6-P dehydrogenase assay |
| `CHEMBL3582523` | B | 9 | 17 | 16 | T1/2:17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as dissociation half life by SPR method |
| `CHEMBL3583538` | B | 9 | 17 | 16 | K:17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as on rate constant for activator-enzyme binding kinetics by SPR method |
| `CHEMBL3583539` | B | 9 | 17 | 16 | K:17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as off rate constant for activator-enzyme binding kinetics by SPR method |
| `CHEMBL3583540` | B | 9 | 17 | 16 | Kd:17 | Binding affinity to biotinylated human recombinant glucokinase expressed in Escherichia coli assessed as dissociation rate constant by SPR method |
| `CHEMBL3583541` | B | 9 | 17 | 16 | EC50:17 | Activation of recombinant human glucokinase assessed as formation of glucose-6-phosphate by G6PDH/NADP coupled assay |
| `CHEMBL5261510` | B | 9 | 17 | 17 | EC50:17 | Activation of recombinant human glucokinase assessed as increase in maximal activation by measuring accumulation of NADH in presence of S0.5 of glucose and NAD by plate reader analysis |
| `CHEMBL1065927` | B | 9 | 16 | 16 | EC50:16 | Activation of human liver glucokinase expressed in CHO cells at 2.5 mM glucose concentration by glucose-6-phosphate coupled continuous spectrophotometric assay |
| `CHEMBL1038417` | B | 9 | 16 | 16 | EC50:16 | Activation of human liver glucokinase expressed in CHO cells at 10 mM glucose concentration by glucose-6-phosphate coupled continuous spectrophotometric assay |
| `CHEMBL1038416` | B | 9 | 14 | 14 | Emax:14 | Activation of human liver glucokinase expressed in CHO cells at 2.5 mM glucose concentration by glucose-6-phosphate coupled continuous spectrophotometric assay relative to 2-amino-N-(4-methyl-1,3-thiazol-2-yl)-5-[(4-methyl-4H-1,2,4-triazol-3-yl)sulfanyl]benzamide |
| `CHEMBL1038418` | B | 9 | 14 | 14 | Emax:14 | Activation of human liver glucokinase expressed in CHO cells at 10 mM glucose concentration by glucose-6-phosphate coupled continuous spectrophotometric assay relative to 2-amino-N-(4-methyl-1,3-thiazol-2-yl)-5-[(4-methyl-4H-1,2,4-triazol-3-yl)sulfanyl]benzamide |
| `CHEMBL4014754` | B | 9 | 14 | 9 | EC50:14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 2.5 mM glucose by G6PDH coupled spectrophotometric assay |
| `CHEMBL4014755` | B | 9 | 14 | 9 | %max:14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 2.5 mM glucose by G6PDH coupled spectrophotometric assay relative to control |
| `CHEMBL4014756` | B | 9 | 14 | 9 | EC50:14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 10 mM glucose by G6PDH coupled spectrophotometric assay |
| `CHEMBL4014757` | B | 9 | 14 | 9 | %max:14 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 10 mM glucose by G6PDH coupled spectrophotometric assay relative to control |
| `CHEMBL990114` | B | 9 | 13 | 13 | EC50:13 | Activation of GST fused human liver glucokinase assessed as generation of NADPH by G6PDH coupled assay |
| `CHEMBL2168593` | B | 9 | 13 | 13 | EC50:13 | Activation of human recombinant glucokinase using 6.5 mM glucose by spectrophotometric analysis |
| `CHEMBL990115` | B | 9 | 12 | 12 | FC:12 | Activation of GST fused human liver glucokinase assessed as generation of NADPH by G6PDH coupled assay relative to control |
| `CHEMBL990116` | B | 9 | 12 | 12 | Activity:12 | Activation of GST fused human liver glucokinase assessed as concentration required to double glucokinase activity by G6PDH coupled assay |
| `CHEMBL1035543` | B | 8 | 11 | 11 | EC50:11 | Activation of glucokinase |
| `CHEMBL4624346` | B | 9 | 9 | 9 | EC50:9 | Activation of human glucokinase |
| `CHEMBL5733424` | B | 9 | 9 | 3 | IC50:3; k_off:3; kon:3 | Inhibition Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK, and other kinases were obtained using an Invitrogen Select Screening assay as known in the art. |
| `CHEMBL5735123` | B | 9 | 9 | 3 | IC50:3; k_off:3; kon:3 | In Vitro Activity Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK and other kinases were obtained using an Invitrogen Select Screening assay as known in the art. |
| `CHEMBL4014791` | B | 9 | 7 | 7 | EC50:7 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 10 mM glucose by G6PDH coupled spectrophotometric assay |
| `CHEMBL4014792` | B | 9 | 7 | 7 | %max:7 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate in presence of 10 mM glucose by G6PDH coupled spectrophotometric assay relative to control |
| `CHEMBL4332276` | B | 9 | 7 | 7 | IC50:7 | Inhibition of recombinant human N-terminal GST-fused glucokinase expressed in Escherichia coli using glucose-6-phosphatedehydrogenase as substrate by spectrophotometric assay |
| `CHEMBL4614239` | B | 9 | 7 | 7 | EC50:7 | Activation of glucokinase (unknown origin) in presence of 4% human serum albumin |
| `CHEMBL865109` | F | 9 | 6 | 6 | EC50:6 | Activation of glucokinase |
| `CHEMBL1058318` | B | 9 | 6 | 6 | Emax:6 | Activation of flag-tagged recombinant human liver glucokinase expressed in Escherichia coli assessed as glucose-6-phosphate dehydrogenase activity at 30 uM spectrophotometry relative to 2-amino-5-(4-methyl-4H-1,2,4-triazol-3-ylthio)-N-(4-methylthiazol-2-yl)benzamide |
| `CHEMBL1042962` | B | 9 | 6 | 6 | EC50:6 | Activation of human glucokinase expressed in Escherichia coli BL21(DE3) by G6PDH-coupled spectrometry |
| `CHEMBL3095515` | B | 9 | 6 | 6 | EC50:6 | Activation of recombinant human glucokinase by G6PDH/NADP coupled assay |
| `CHEMBL4428660` | B | 9 | 6 | 6 | Kd:6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) by stopped-flow fluorometric method |
| `CHEMBL4428665` | B | 9 | 6 | 6 | K:6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as compound association rate constant by stopped-flow fluorometric method |
| `CHEMBL4428667` | B | 9 | 6 | 6 | K:6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as compound dissociation rate constant by stopped-flow fluorometric method |
| `CHEMBL4428668` | B | 9 | 6 | 6 | TIME:6 | Displacement of TAFMT from human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as compound residence time by stopped-flow fluorometric method |
| `CHEMBL4804384` | B | 9 | 6 | 6 | % inhibition:6 | Inhibition of GCK in human NCI-H929 cells by mass spectroscopic analysis |
| `CHEMBL1825594` | B | 9 | 5 | 5 | Activity:5 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 5 mM glucose by G6PDH coupled continuous spectrophotometric assay relative to Ro-28-1675 |
| `CHEMBL1825595` | B | 9 | 5 | 5 | Activity:5 | Activation of human recombinant glucokinase expressed in Escherichia coli coexpressing G6PDH at 10 uM using 20 mM glucose by G6PDH coupled continuous spectrophotometric assay relative to Ro-28-1675 |
| `CHEMBL3095513` | B | 9 | 5 | 5 | FC:5 | Activation of recombinant human glucokinase assessed as increase in enzyme Vmax for glucose by G6PDH/NADP coupled assay relative to control |
| `CHEMBL3095514` | B | 9 | 5 | 5 | FC:5 | Activation of recombinant human glucokinase assessed as decrease in enzyme Km for glucose by G6PDH/NADP coupled assay relative to control |
| `CHEMBL3224178` | B | 9 | 5 | 5 | EC50:5 | Activation of glucokinase (unknown origin) |
| `CHEMBL3224190` | B | 9 | 5 | 5 | EC50:5 | Activation of glucokinase (unknown origin) |
| `CHEMBL4427376` | B | 9 | 5 | 5 | Ratio EC50:5 | Ratio of EC50 for activation of human recombinant glucokinase in presence of 4% HSA to EC50 for activation of human recombinant glucokinase in absence of 4% HSA |
| `CHEMBL4428666` | B | 9 | 5 | 5 | Activity:5 | Activation of human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as conversion of glucose to glucose-6-phosphate by G6PDH/NADP-coupled enzyme assay |
| `CHEMBL4804566` | B | 9 | 5 | 5 | Inhibition:5 | Inhibition of GCK in human NCI-H929 cells at 10 uM by mass spectroscopic analysis relative to control |
| `CHEMBL1942767` | B | 9 | 4 | 4 | Inhibition:4 | Inhibition of human GCK in HL-60 cells lysate assessed as reduction of labeling of acyl-phosphate ATP probe at 100 nM |
| `CHEMBL4193353` | B | 9 | 4 | 4 | EC50:4 | Activation of recombinant human liver glucokinase 2 assessed as assessed as reduction in Km for glucose in presence of 5 mM glucose by G6PDH coupled assay |
| `CHEMBL2412836` | B | 9 | 3 | 3 | Activity:3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 assessed as increase in glucose Vmax by enzyme-kinetic study in presence of 200 mM glucose |
| `CHEMBL2412837` | B | 9 | 3 | 3 | Activity:3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 assessed as increase in glucose kcat by enzyme-kinetic study in presence of 200 mM glucose |
| `CHEMBL2412840` | B | 9 | 3 | 3 | Activity:3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 assessed as reduction in transition time at 20 uM by spectrophotometry in absence of glucose |
| `CHEMBL2412842` | B | 9 | 3 | 3 | Activity:3 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 assessed as reduction in lag of transient-state enzyme progress curve at 25 to 100 uM by spectrophotometry in absence of glucose |
| `CHEMBL3630507` | B | 9 | 3 | 3 | Inhibition:3 | Inhibition of GCK (unknown origin) at 10 uM after 120 mins P33 radiolabeled kinase activity assay |
| `CHEMBL4118352` | B | 9 | 3 | 3 | Inhibition:3 | Inhibition of GCK in human SKCO1 cells at 1 uM after 4 hrs using biotin labeled DIKGANLLLTLQGDVK probe by mass-spectrometric analysis relative to control |
| `CHEMBL4349275` | B | 9 | 3 | 3 | FC:3 | Activation of glucokinase (unknown origin) at 10 uM relative to control |
| `CHEMBL4427586` | B | 9 | 3 | 3 | Ratio EC50:3 | Ratio of EC50 for activation of human recombinant glucokinase in presence of 4% human serum albumin to EC50 for activation of human recombinant glucokinase in absence of 4% human serum |
| `CHEMBL3131525` | B | 9 | 2 | 2 | EC50:2 | Activation of human glucokinase |
| `CHEMBL3131380` | B | 9 | 2 | 2 | EC50:2 | Activation of glucokinase (unknown origin) using glucose as substrate |
| `CHEMBL3369166` | B | 9 | 2 | 2 | Ratio:2 | Ratio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepatocytes to IC50 for inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP |
| `CHEMBL3751341` | B | 9 | 2 | 2 | Kd:2 | Binding affinity to human glucokinase by surface plasmon resonance analysis |
| `CHEMBL4327702` | B | 9 | 2 | 2 | Activity:2 | Inhibition of human GCK assessed as residual activity at 1 uM using MBP as substrate by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4674372` | B | 9 | 2 | 2 | Activity:2 | Inhibition of human GCK using MBP as substrate assessed as residual activity at 1 uM by [gamma-33P]-ATP assay relative to control |
| `CHEMBL5048855` | B | 9 | 2 | 2 | EC50:2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli incubated for 10 mins in presence of 5 mM glucose by plate reader analysis |
| `CHEMBL5048925` | B | 9 | 2 | 2 | IC50:2 | Displacement of fluorescent labeled derivative from recombinant human hepatic glucokinase incubated for 30 mins in presence of 12 mM glucose by fluorescent polarization assay |
| `CHEMBL5048926` | B | 9 | 2 | 2 | Emax:2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli assessed as maximal activation incubated for 10 mins in presence of 5 mM glucose by plate reader analysis relative to control |
| `CHEMBL5048927` | B | 9 | 2 | 2 | EC50:2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli incubated for 10 mins in presence of 2 mM glucose by plate reader analysis |
| `CHEMBL5048928` | B | 9 | 2 | 2 | Emax:2 | Activation of recombinant human full length liver glucokinase expressed in Escherichia coli assessed as maximal activation incubated for 10 mins in presence of 2 mM glucose by plate reader analysis relative to control |
| `CHEMBL5048931` | B | 9 | 2 | 2 | Km:2 | Activation of human glucokinase assessed as decrease in glucose Km |
| `CHEMBL5048932` | B | 9 | 2 | 2 | Kcat:2 | Activation of human glucokinase assessed as increase rate of glucose turnover |
| `CHEMBL5464286` | B | 9 | 2 | 2 | % Activity remaining:2 | % Activity remaining of GCK in the Dundee kinase panel at 1.0 µM |
| `CHEMBL6194772` | B | 9 | 2 | 1 | Effect:2 | Effect of GCK(h) at compound concentration of 1.0 uM using the Cerep Kinase panel |
| `CHEMBL1054393` | B | 9 | 1 | 1 | Vmax:1 | Activity at human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) |
| `CHEMBL1054394` | B | 9 | 1 | 1 | Activity:1 | Activation of His-tagged human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) assessed as drug level required for half-maximal activation for sigmoid activation curve |
| `CHEMBL2025596` | B | 9 | 1 | 1 | Ratio EC50:1 | Ratio of EC50 for GK50 to compound for activation of human recombinant glucokinase expressed in Escherichia coli BL21(DE3) coexpressing G6PDH assessed as glucose 6-phosphate formation by spectrometric analysis |
| `CHEMBL2167261` | B | 9 | 1 | 1 | Activity:1 | Activation of human recombinant glucokinase |
| `CHEMBL2168355` | B | 9 | 1 | 1 | Activity:1 | Activity at human recombinant glucokinase assessed as increase in Vmax at 1 uM |
| `CHEMBL2168357` | B | 9 | 1 | 1 | Activity:1 | Activity at human recombinant glucokinase assessed as decrease in Km at 1 uM |
| `CHEMBL2318453` | B | 9 | 1 | 1 | EC50:1 | Induction of human hepatic glucokinase activity at 5 mM glucose concentration |
| `CHEMBL2318648` | B | 9 | 1 | 1 | Vmax:1 | Induction of human pancreatic glucokinase activity assessed as maximal reaction rate at 5 mM glucose concentration |
| `CHEMBL2318650` | B | 9 | 1 | 1 | Activity:1 | Induction of human pancreatic glucokinase activity assessed as half maximal saturation concentration of glucose at 5 mM glucose concentration |
| `CHEMBL2346006` | B | 9 | 1 | 1 | Activity:1 | Activation of human recombinant glucokinase using glucose as substrate |
| `CHEMBL2412838` | B | 9 | 1 | 1 | Kd:1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 at 11 uM by isothermal titration calorimetry in presence of 200 mM glucose |
| `CHEMBL2412839` | B | 9 | 1 | 1 | Kd:1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 at 100 uM by isothermal titration calorimetry in presence of 200 mM glucose |
| `CHEMBL2412841` | B | 9 | 1 | 1 | Activity:1 | Allosteric activation of recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 assessed as reduction in lag of transient-state enzyme progress curve by spectrophotometry in absence of glucose |
| `CHEMBL2412843` | B | 9 | 1 | 1 | Kd:1 | Binding affinity to recombinant wild-type human pancreatic glucokinase expressed in Escherichia coli K-12 by isothermal titration calorimetry in presence of 200 mM glucose |
| `CHEMBL3119606` | B | 9 | 1 | 1 | Activity:1 | Binding affinity to biotinylated GK (unknown origin) by SPR analysis |
| `CHEMBL3239305` | B | 9 | 1 | 1 | Vmax:1 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal activity |
| `CHEMBL3377935` | B | 9 | 1 | 1 | EC50:1 | Activation of human recombinant Glucokinase measured over 5 mins by G6-PD coupled assay in presence of 5 mM glucose and 4% HSA |
| `CHEMBL3578772` | B | 9 | 1 | 1 | Vmax:1 | Activation of recombinant human pancreatic glucokinase using 10 mM glucose as substrate assessed as enzyme Vmax by G6PDH coupled assay |
| `CHEMBL3750455` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL3751340` | B | 9 | 1 | 1 | EC50:1 | Effect on human glucokinase activity after 60 mins by luciferase-based luminescence assay in absence of human GKRP |
| `CHEMBL3807041` | B | 9 | 1 | 1 | IC50:1 | Inhibition of HK4 (unknown origin) |
| `CHEMBL3829755` | B | 9 | 1 | 1 | Activity:1 | Inhibition of human GCK (2 to 812 residues) assessed as remaining enzyme activity at 50 uM after 30 mins by 33P-ATP filter-binding assay |
| `CHEMBL4022295` | B | 9 | 1 | 1 | Ratio EC50:1 | Potency index, ratio of (S)-N-(5-chlorothiazol-2-yl)-2-(4-(cyclopropylsulfonyl)-6-methyl-2-oxopyridin-1(2H)-yl)-3-(tetrahydro-2H-pyran-4-yl)propanamide EC50 to compound EC50 for full length human C-terminal FLAG-tagged glucokinase (12 to 465 residues) |
| `CHEMBL4034328` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK Lysine 1 labelling site (unknown origin) at 10 uM |
| `CHEMBL4034329` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK Lysine 2 labelling site (unknown origin) at 10 uM |
| `CHEMBL4045684` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 1 uM |
| `CHEMBL4045685` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 1 uM |
| `CHEMBL4045926` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 0.1 uM |
| `CHEMBL4045927` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 0.1 uM |
| `CHEMBL4050296` | B | 9 | 1 | 1 | Activity:1 | Inhibition of human GCK assessed as residual activity at 1 uM in presence of 33P-ATP by filter-binding assay relative to control |
| `CHEMBL4057298` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4057494` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 0.49 uM relative to control |
| `CHEMBL4120341` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4120568` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 0.19 uM relative to control |
| `CHEMBL4120792` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 0.164 uM relative to control |
| `CHEMBL4328066` | B | 9 | 1 | 1 | Activity:1 | Inhibition of human GCK assessed as residual activity at 100 uM using MBP as substrate by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4328430` | B | 9 | 1 | 1 | IC50:1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4349277` | B | 9 | 1 | 1 | FC:1 | Activation of recombinant human glucokinase assessed as conversion of D-glucose to D-glucose-6-phosphate dehydrate at 10 uM by G6PDH/NADP coupled assay relative to control |
| `CHEMBL4352383` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of human GCK at 10 uM using MBP as substrate in presence of [gamma-33P]-ATP |
| `CHEMBL4379338` | B | 9 | 1 | 1 | IC50:1 | Inhibition of recombinant human GCK (1 to 473 residues) using myelin basic protein as substrate after 40 mins in presence of [gamma-33ATP] by radiometric scintillation counting analysis |
| `CHEMBL4428661` | B | 9 | 1 | 1 | Kd:1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) by stopped-flow fluorometric method |
| `CHEMBL4428662` | B | 9 | 1 | 1 | K:1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as association rate constant by stopped-flow fluorometric method |
| `CHEMBL4428663` | B | 9 | 1 | 1 | K:1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as dissociation rate constant by stopped-flow fluorometric method |
| `CHEMBL4428664` | B | 9 | 1 | 1 | TIME:1 | Binding affinity to human pancreatic N-terminal His6-tagged glucokinase isoform 1 expressed in Escherichia coli BL21(DE3) assessed as residence time in presence of glucose by stopped-flow fluorometric method |
| `CHEMBL4674071` | B | 9 | 1 | 1 | IC50:1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4674112` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of human GCK using MBP as substrate at 0.5 uM by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4681971` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK human PBMC lysates at 1 uM by ActivX screen assay relative to control |
| `CHEMBL4720720` | B | 9 | 1 | 1 | IC50:1 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4721072` | B | 9 | 1 | 1 | Activity:1 | Inhibition of human GCK assessed as residual activity using MBP as substrate at 1 uM by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4880638` | B | 9 | 1 | 1 | Inhibition:1 | GCK (h) Millipore kinase activity assay |
| `CHEMBL4881060` | B | 9 | 1 | 1 | % residual kinase activity:1 | GCK(h) Eurofins Kinase panel |
| `CHEMBL4884319` | B | 9 | 1 | 1 | % residual kinase activity:1 | GCK(h) Eurofins kinase panel |
| `CHEMBL4884756` | B | 9 | 1 | 1 | IC50:1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885047` | B | 9 | 1 | 1 | IC50:1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885337` | B | 9 | 1 | 1 | IC50:1 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4887878` | B | 9 | 1 | 1 | Inhibition:1 | GCK(h) Millipore kinase panel |
| `CHEMBL5038899` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of human recombinant GCK assessed as reduction in substrate phosphorylation at 50 to 500 nM using ATP and Ulight-CKKSRGDYMTMQIG (IRS-1) incubated for 60 min by LANCE detection method |
| `CHEMBL5058199` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL5059136` | B | 9 | 1 | 1 | Inhibition:1 | GCK(h) Kinase panel |
| `CHEMBL5156973` | B | 9 | 1 | 1 | IC50:1 | Inhibition of human recombinant GCK (1 to 473 residues) using myelin basic protein as substrate in presence of ATP incubated for 40 mins by radiometric based scintillation method |
| `CHEMBL5261324` | B | 9 | 1 | 1 | FC:1 | Agonist activity at glucokinase (unknown origin) assessed as fold increase at 10 uM by luciferase reporter gene assay |
| `CHEMBL5261513` | B | 9 | 1 | 1 | Activity:1 | Activation of recombinant human glucokinase assessed as increase in maximal activation by measuring accumulation of NADH at 60 uM in presence of S0.5 of glucose and NAD by plate reader analysis |
| `CHEMBL5464035` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464534` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464912` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5657730` | B | 9 | 1 | 1 | Activity:1 | Inhibition of human GCK assessed as remaining activity at 10 uM in presence of ATP by radiometric kinase assay relative to control |
| `CHEMBL5679451` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infected Sf21 cells at 0.3 uM by filter binding assay relative to control |
| `CHEMBL5679684` | B | 9 | 1 | 1 | Inhibition:1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infected Sf21 cells at 3 uM by filter binding assay relative to control |
| `CHEMBL5681798` | B | 9 | 1 | 1 | Activity:1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infected Sf21 cells assessed as residual activity at 0.1 uM relative to control |
| `CHEMBL5682058` | B | 9 | 1 | 1 | Activity:1 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infected Sf21 cells assessed as residual activity at 1 uM relative to control |
| `CHEMBL5724151` | B | 9 | 1 | 1 | activity:1 | activity of GCK(h) at 1.0 µM in the Eurofins Kinase panel |

