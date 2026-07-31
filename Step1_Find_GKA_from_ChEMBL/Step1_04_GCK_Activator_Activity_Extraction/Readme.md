从已确认的 GCK 激活 assay 中提取所有 activity，关联对应小分子，并分别根据：

效力：EC50、AC50 等
效能：% activation、fold activation 等

评价每个分子的 GCK 激活强弱。

分开评价，不先合成一个总分：

效力：EC50/AC50 越小越强；或 pActivity 越大越强。
效能：最大 % activation 或 fold activation 越大越强。
证据可靠性：重复 assay、独立文献更多，且结果一致，证据更强。

最终给每个分子三项结果：效力、效能、证据。

这一步不对activity进行排序。

主要利用activity表。主要应包括：
分子 ID 和结构
assay ID
激活效力：EC50、AC50、Potency
激活效能：% activation、fold activation
relation、单位、pChEMBL
assay 置信度
数据质量标记
文献和来源
同一分子的独立 assay 数、独立文献数

代码：
Step1_04_GCK_Activator_Activity_Extraction.py
