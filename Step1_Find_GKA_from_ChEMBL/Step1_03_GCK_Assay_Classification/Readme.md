逐个判断 228 个 GCK assay 实际测量的内容，并分为：

GCK 激活
GCK 抑制
GCK 结合
GCK–GKRP 相互作用
细胞或表型效应
无法判断
主要依据

target_chembl_id + assay_description + assay_type + BAO format + 该 assay 下的 activity standard_type

查询路径

assay_id → activities.assay_id → 汇总 standard_type、units、relation 和 activity_comment

输出

在 assay 主表中增加：

assay_category
classification_confidence
classification_reason
review_required

这一步完成后，才能确定哪些 assay 应用于搜索真正的 GKA 小分子。

1, 先用规则分类明确记录
例如描述中明确出现 activation、activator、inhibition、binding、GKRP interaction。
代码：Step1_03_GCK_Assay_Classification.py

2, 再让 LLM 处理模糊记录
输入 assay_description + standard_type + units + activity_comment，让它判断实验测的是激活、抑制、结合还是无法确定。
3, LLM 必须输出证据句和置信度
不能只给标签，还要指出依据哪段描述做出判断。
4, 无法明确判断的记录进入人工审核
不允许 LLM 强行分类。

因此，LLM 的角色是：

处理文本语义和模糊 assay，不负责数据库筛选，也不替代明确规则。

对 228 个 assay，建议先规则分类，再只把剩余模糊项交给 LLM。
