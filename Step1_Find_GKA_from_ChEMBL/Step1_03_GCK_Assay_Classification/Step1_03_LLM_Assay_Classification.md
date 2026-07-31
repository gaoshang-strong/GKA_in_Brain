# Step1_03 LLM 分类结果（第 2 阶段）

- 模型：`deepseek-v4-pro`
- 运行时间：2026-07-30 22:37:56
- 交给 LLM 的记录：**19** / 228
- API 调用 0 次（缓存命中 19 次），token 输入 11,927 / 输出 10,607

> LLM 只处理规则阶段判不明确的记录，不覆盖规则的 high 置信度结论。每条结论都要求给出**逐字引自原描述**的证据句；脚本会校验该证据是否真实存在，核验不通过的一律不采纳标签并转人工。

## 证据校验

- 通过：**19**
- 未通过（证据句无法在原文中找到，标签不予采纳）：**0**

## 最终分类分布

| 类别 | assay 数 | 其中仍需人工审核 |
| --- | ---: | ---: |
| GCK 激活 | 142 | 4 |
| GCK 抑制 | 54 | 50 |
| GCK 结合 | 18 | 0 |
| 无法判断 | 10 | 10 |
| GCK–GKRP 相互作用 | 3 | 0 |
| 细胞或表型效应 | 1 | 1 |

## 判定方式

| final_method | 数量 | 含义 |
| --- | ---: | --- |
| rule | 209 | 规则明确判定，未经 LLM |
| manual_review_required | 10 | 规则与 LLM 都未能确定，需人工审核 |
| llm | 9 | 规则不明确，采纳了 LLM 的判定（证据已核验） |

## LLM 逐条结果

| assay | 规则暂定 | LLM 判定 | 置信度 | 证据核验 | 证据句 | LLM 理由 |
| --- | --- | --- | --- | --- | --- | --- |
| `CHEMBL1054393` | 无法判断 | 无法判断 | low | TRUE | Activity at human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) | 描述仅提及活性测定且实测指标为 Vmax，未指明化合物对 GCK 的作用方向（激活或抑制），无法明确分类。 |
| `CHEMBL2168355` | 无法判断 | GCK 激活 | high | TRUE | increase in Vmax | Vmax最大反应速率的增加直接反映化合物增强酶活性，属于GCK激活。 |
| `CHEMBL2168357` | 无法判断 | 无法判断 | low | TRUE | Activity at human recombinant glucokinase assessed as decrease in Km at 1 uM | 描述仅提及降低 Km，未指明该变化对应激活还是抑制，方向不明确。 |
| `CHEMBL2432300` | 无法判断 | 无法判断 | low | TRUE | Activity of human glucokinase at 50 uM | 描述仅提及在50 uM某种条件下的人葡萄糖激酶活性，未说明该条件是指化合物还是底物，也未提及测量的是激活、抑制或结合，因此无法归类。 |
| `CHEMBL3117036` | GCK 抑制 | GCK–GKRP 相互作用 | high | TRUE | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP | 实验直接测量荧光标记GK与生物素标记GKRP之间的相互作用被化合物抑制的情况，属于GCK–GKRP相互作用检测。 |
| `CHEMBL3226993` | GCK 激活 | 无法判断 | high | TRUE | Activity at glucokinase (unknown origin) | 实验描述过于简略，仅提及对葡萄糖激酶的活性，未指明激活或抑制，且assay_type标记为B但实测standard_type均为EC50，存在不一致，无法可靠判断测量内容。 |
| `CHEMBL3223676` | GCK 抑制 | GCK 抑制 | high | TRUE | Inhibition of recombinant human pancreatic glucokinase by G-6-P dehydrogenase assay | 实验描述直接使用“Inhibition of ... glucokinase”，明确指向抑制活性检测。尽管实测参数为EC50，但抑制实验中也常使用EC50表示半数有效抑制浓度，描述 |
| `CHEMBL3239304` | GCK 激活 | 无法判断 | low | TRUE | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal satu | 描述仅提及测量葡萄糖半最大饱和浓度，未提及化合物处理或作用方向，无法确定实验目的。 |
| `CHEMBL3239305` | 无法判断 | 无法判断 | high | TRUE | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal acti | 实验描述仅评估了重组人胰腺葡萄糖激酶对底物葡萄糖的活性（半最大活性浓度），未涉及任何测试化合物或调控方向，因此无法归类到给定类别中。 |
| `CHEMBL3377919` | 无法判断 | 无法判断 | high | TRUE | Activity of human recombinant Glucokinase measured over 5 mins by spectrophotometry | 描述仅提及重组GCK活性测量，未提及任何化合物处理，standard_type为Km，属于酶学动力学参数，无法判断是否评估化合物对GCK的激活或抑制等效应 |
| `CHEMBL3364949` | GCK 抑制 | GCK–GKRP 相互作用 | high | TRUE | Inhibition of fluorescein-labeled human GK interaction with biotin-labeled human GKRP | 实验描述明确指出测量的是GCK与GKRP之间相互作用的抑制，属于蛋白-蛋白相互作用实验，而非GCK酶活性直接测定。 |
| `CHEMBL3369166` | GCK 抑制 | 无法判断 | low | TRUE | Ratio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepatocytes to IC50  | 实验同时涉及细胞表型（GK转位的EC50）和GCK–GKRP相互作用抑制（IC50），最终读取值为二者比值，无法归为单一类别。 |
| `CHEMBL3751340` | GCK 激活 | GCK 激活 | high | TRUE | Effect on human glucokinase activity after 60 mins by luciferase-based luminescence assay  | 实验明确测量对 GCK 酶活性的影响，standard_type 仅含 EC50，表明化合物增强酶活，且条件不含 GKRP，排除相互作用检测。 |
| `CHEMBL3888434` | GCK 激活 | GCK 激活 | medium | TRUE | Enzymatic activity is calculated from the initial reaction velocity, measured from the cha | 该实验测定重组人GCK的酶活性，标准化类型为EC50，表明测量的是化合物增强酶活性的能力，因此归类为GCK激活。 |
| `CHEMBL4022295` | GCK 激活 | 无法判断 | low | TRUE | Potency index, ratio of (S)-N-(5-chlorothiazol-2-yl)-2-(4-(cyclopropylsulfonyl)-6-methyl-2 | 描述仅提及EC50比值，未说明功能方向（激活或抑制），且assay_type标注为Binding但与EC50不符，信息不足以归类至具体测量类别。 |
| `CHEMBL5048925` | GCK 结合 | GCK 结合 | high | TRUE | Displacement of fluorescent labeled derivative from recombinant human hepatic glucokinase | 描述为基于荧光偏振的竞争结合实验，直接测量化合物从重组GK蛋白上置换荧光标记探针的能力，属于GCK结合测定。 |
| `CHEMBL5261324` | GCK 激活 | 细胞或表型效应 | medium | TRUE | by luciferase reporter gene assay | 该实验使用荧光素酶报告基因检测，读数反映的是细胞水平的响应，而非直接测量GCK酶活性，因此属于细胞或表型效应。 |
| `CHEMBL5733424` | GCK 抑制 | 无法判断 | high | TRUE | Inhibition Assay: The in vitro activity of the compounds described herein in inhibiting TA | 实验描述明确测量的是TAK1、HCK等激酶的抑制活性，未提及葡萄糖激酶（GCK），无法判断该实验是否实际测量了GCK的抑制或其他效应。 |
| `CHEMBL5735123` | GCK 抑制 | GCK 抑制 | medium | TRUE | the in vitro activity of the compounds described herein in inhibiting TAK1, HCK and other  | 实验描述明确提及测定化合物对TAK1、HCK及其他激酶的抑制活性，结合靶点为GCK，推断为GCK抑制实验。但未直接提及GCK，故置信度中等。 |

## 仍需人工审核

共 **65** 条。按 README 第 4 步，这些不允许由 LLM 强行分类。

| assay | 暂定类别 | 描述 |
| --- | --- | --- |
| `CHEMBL1054393` | 无法判断 | Activity at human recombinant liver glucokinase expressed in Escherichia coli BL21 (DE3) |
| `CHEMBL1942767` | GCK 抑制 | Inhibition of human GCK in HL-60 cells lysate assessed as reduction of labeling of acyl-phosphate ATP probe at |
| `CHEMBL2168357` | 无法判断 | Activity at human recombinant glucokinase assessed as decrease in Km at 1 uM |
| `CHEMBL2432300` | 无法判断 | Activity of human glucokinase at 50 uM |
| `CHEMBL3226993` | GCK 激活 | Activity at glucokinase (unknown origin) |
| `CHEMBL3239304` | GCK 激活 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal saturation concentration |
| `CHEMBL3239305` | 无法判断 | Activity of recombinant human pancreatic glucokinase assessed as glucose half-maximal activity |
| `CHEMBL3377919` | 无法判断 | Activity of human recombinant Glucokinase measured over 5 mins by spectrophotometry |
| `CHEMBL3369166` | GCK 抑制 | Ratio of EC50 for GK translocation from nucleus to cytoplasm of mouse hepatocytes to IC50 for inhibition of fl |
| `CHEMBL3630507` | GCK 抑制 | Inhibition of GCK (unknown origin) at 10 uM after 120 mins P33 radiolabeled kinase activity assay |
| `CHEMBL3750455` | GCK 抑制 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL3829755` | GCK 抑制 | Inhibition of human GCK (2 to 812 residues) assessed as remaining enzyme activity at 50 uM after 30 mins by 33 |
| `CHEMBL3888434` | GCK 激活 | Coupled Enzymatic Assay: The assay is carried out according to the protocol outlined in Hariharan et al (1997) |
| `CHEMBL4022295` | GCK 激活 | Potency index, ratio of (S)-N-(5-chlorothiazol-2-yl)-2-(4-(cyclopropylsulfonyl)-6-methyl-2-oxopyridin-1(2H)-yl |
| `CHEMBL4034328` | GCK 抑制 | Inhibition of GCK Lysine 1 labelling site (unknown origin) at 10 uM |
| `CHEMBL4034329` | GCK 抑制 | Inhibition of GCK Lysine 2 labelling site (unknown origin) at 10 uM |
| `CHEMBL4045684` | GCK 抑制 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 1 uM |
| `CHEMBL4045685` | GCK 抑制 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 1 uM |
| `CHEMBL4045926` | GCK 抑制 | Inhibition of GCK conserved Lys1 (DTVTSELAAVKIVK) in human PBMC at 0.1 uM |
| `CHEMBL4045927` | GCK 抑制 | Inhibition of GCK conserved Lys2 (DIKGANLLLTLQGDVK) in human PBMC at 0.1 uM |
| `CHEMBL4050296` | GCK 抑制 | Inhibition of human GCK assessed as residual activity at 1 uM in presence of 33P-ATP by filter-binding assay r |
| `CHEMBL4057298` | GCK 抑制 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4057494` | GCK 抑制 | Inhibition of GCK (unknown origin) at 0.49 uM relative to control |
| `CHEMBL4118352` | GCK 抑制 | Inhibition of GCK in human SKCO1 cells at 1 uM after 4 hrs using biotin labeled DIKGANLLLTLQGDVK probe by mass |
| `CHEMBL4120341` | GCK 抑制 | Inhibition of GCK (unknown origin) at 0.45 uM relative to control |
| `CHEMBL4120568` | GCK 抑制 | Inhibition of GCK (unknown origin) at 0.19 uM relative to control |
| `CHEMBL4120792` | GCK 抑制 | Inhibition of GCK (unknown origin) at 0.164 uM relative to control |
| `CHEMBL4327702` | GCK 抑制 | Inhibition of human GCK assessed as residual activity at 1 uM using MBP as substrate by [gamma-33P]-ATP assay  |
| `CHEMBL4328066` | GCK 抑制 | Inhibition of human GCK assessed as residual activity at 100 uM using MBP as substrate by [gamma-33P]-ATP assa |
| `CHEMBL4328430` | GCK 抑制 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4352383` | GCK 抑制 | Inhibition of human GCK at 10 uM using MBP as substrate in presence of [gamma-33P]-ATP |
| `CHEMBL4379338` | GCK 抑制 | Inhibition of recombinant human GCK (1 to 473 residues) using myelin basic protein as substrate after 40 mins  |
| `CHEMBL4674071` | GCK 抑制 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4674112` | GCK 抑制 | Inhibition of human GCK using MBP as substrate at 0.5 uM by [gamma-33P]-ATP assay relative to control |
| `CHEMBL4674372` | GCK 抑制 | Inhibition of human GCK using MBP as substrate assessed as residual activity at 1 uM by [gamma-33P]-ATP assay  |
| `CHEMBL4681971` | GCK 抑制 | Inhibition of GCK human PBMC lysates at 1 uM by ActivX screen assay relative to control |
| `CHEMBL4720720` | GCK 抑制 | Inhibition of human GCK using MBP as substrate by [gamma-33P]-ATP assay |
| `CHEMBL4721072` | GCK 抑制 | Inhibition of human GCK assessed as residual activity using MBP as substrate at 1 uM by [gamma-33P]-ATP assay  |
| `CHEMBL4804384` | GCK 抑制 | Inhibition of GCK in human NCI-H929 cells by mass spectroscopic analysis |
| `CHEMBL4804566` | GCK 抑制 | Inhibition of GCK in human NCI-H929 cells at 10 uM by mass spectroscopic analysis relative to control |
| `CHEMBL4880638` | GCK 抑制 | GCK (h) Millipore kinase activity assay |
| `CHEMBL4881060` | 无法判断 | GCK(h) Eurofins Kinase panel |
| `CHEMBL4884319` | 无法判断 | GCK(h) Eurofins kinase panel |
| `CHEMBL4884756` | GCK 抑制 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885047` | GCK 抑制 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4885337` | GCK 抑制 | GCK(M4K2LGY1) Takeda global kinase panel |
| `CHEMBL4887878` | GCK 抑制 | GCK(h) Millipore kinase panel |
| `CHEMBL5038899` | GCK 抑制 | Inhibition of human recombinant GCK assessed as reduction in substrate phosphorylation at 50 to 500 nM using A |
| `CHEMBL5058199` | GCK 抑制 | Inhibition of GCK (unknown origin) at 1 uM |
| `CHEMBL5059136` | GCK 抑制 | GCK(h) Kinase panel |
| `CHEMBL5156973` | GCK 抑制 | Inhibition of human recombinant GCK (1 to 473 residues) using myelin basic protein as substrate in presence of |
| `CHEMBL5261324` | 细胞或表型效应 | Agonist activity at glucokinase (unknown origin) assessed as fold increase at 10 uM by luciferase reporter gen |
| `CHEMBL5464035` | GCK 抑制 | Inhibition of GCK at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464286` | 无法判断 | % Activity remaining of GCK in the Dundee kinase panel at 1.0 µM |
| `CHEMBL5464534` | GCK 抑制 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5464912` | GCK 抑制 | Inhibition of GCK (h) at 10.0 µM in the Eurofins Kinase panel |
| `CHEMBL5657730` | GCK 抑制 | Inhibition of human GCK assessed as remaining activity at 10 uM in presence of ATP by radiometric kinase assay |
| `CHEMBL5679451` | GCK 抑制 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5679684` | GCK 抑制 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5681798` | GCK 抑制 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5682058` | GCK 抑制 | Inhibition of N-terminal GST- tagged recombinant human GCK (1 to 473 residues) expressed in baculovirus infect |
| `CHEMBL5724151` | 无法判断 | activity of GCK(h) at 1.0 µM in the Eurofins Kinase panel |
| `CHEMBL5733424` | GCK 抑制 | Inhibition Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK, and other k |
| `CHEMBL5735123` | GCK 抑制 | In Vitro Activity Assay: The in vitro activity of the compounds described herein in inhibiting TAK1, HCK and o |
| `CHEMBL6194772` | 无法判断 | Effect of GCK(h) at compound concentration of 1.0 uM using the Cerep Kinase panel |

