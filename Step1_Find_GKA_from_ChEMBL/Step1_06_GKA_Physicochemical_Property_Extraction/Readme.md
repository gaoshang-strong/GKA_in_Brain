根据/ShangGaoAIProjects/GKA_in_Brain/Step1_Find_GKA_from_ChEMBL/Step1_05_GCK_Activator_Ranking_and_Candidate_Selection/Step1_05_Followup_Candidates.csv，把里边分子的理化性质全部从ChEMBL里拉出来。

主要用这几张表：

molecule_dictionary：分子主表，拿 molregno、chembl_id、max_phase 等
compound_structures：SMILES、InChI、InChIKey
compound_properties：MW、ALogP、PSA、HBD、HBA、RTB、QED、Ro5 违规数等
molecule_hierarchy：盐型/母体分子归并
molecule_synonyms：名称、研发代号等补充信息

核心路径：

候选 chembl_id → molecule_dictionary.molregno → compound_structures / compound_properties / molecule_hierarchy

