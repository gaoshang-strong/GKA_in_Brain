通过 UniProt P35557，定位 ChEMBL 中所有对应人葡萄糖激酶 GCK 的靶点记录。

查询策略：
P35557 → component_sequences → component_id → target_components → tid → target_dictionary → chembl_id

输出一张 GCK 靶点映射表，每行代表一个与人 GCK 蛋白组件关联的 ChEMBL target。

至少包含：

accession
component_id
tid
target_chembl_id
pref_name
organism
target_type

使用python，代码名字：
Step1_01_GCK_Target_Anchoring.py
