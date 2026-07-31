通过已确认的 GCK 靶点 ID，提取全部相关 assay，并汇总每个 assay 的实验类型、描述、物种、靶点可信度、来源及活性数据规模。

查询路径：

target_dictionary.tid → assays.tid → assay_id / assay_chembl_id / assay_description

输出一张 GCK assay 清单，先不搜索小分子。核心目标是从 227 个 assay 中找出真正用于测量 GKA 的实验。

使用python，代码名字：
Step1_2_GCK_Assay_Mapping.py