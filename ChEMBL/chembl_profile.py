#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chembl_profile.py — 调查一个 ChEMBL SQLite 数据库，生成 Markdown 概览报告。

设计目标
--------
报告的读者假定为：有生物学 / 生物信息学背景，但**没有药物化学或药理学背景**。
因此报告不只给数字，还会解释每个概念是什么、为什么这么设计、怎么用。

用法
----
    python3 chembl_profile.py /path/to/chembl_37.db
    python3 chembl_profile.py /path/to/chembl_37.db -o report.md
    python3 chembl_profile.py /path/to/chembl_37.db --quick     # 跳过全表扫描，秒级出结果（部分数字为估计值）
    python3 chembl_profile.py /path/to/chembl_37.db --deep      # 额外做 activities×assays×targets 大连接（慢，几分钟）

对 ChEMBL 33+ 的各个版本都应可用：脚本会先探测表 / 列是否存在，缺失的部分自动跳过。
数据库以只读模式打开，不会修改任何内容。
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

# ---------------------------------------------------------------------------
# 知识库：把 ChEMBL 里的行话翻译成生物学家能懂的话
# ---------------------------------------------------------------------------

# 表的一句话说明（覆盖主要表；未收录的表在附录里只给行数）
TABLE_DESC = {
    "activities": "核心事实表：一次测量得到的一个活性数值（如某化合物对某靶点的 IC50）",
    "activity_properties": "某条活性记录的附加参数（如测定时的底物浓度、pH）",
    "activity_supp": "存放补充/原始的多维测定数据（深度数据，通常不需要）",
    "activity_supp_map": "activities 与 activity_supp 的映射",
    "activity_stds_lookup": "标准化字典：哪些活性类型允许被标准化成哪些单位",
    "assays": "实验（assay）描述：这次测量是怎么做的、测的哪个靶点、可信度多高",
    "assay_type": "实验大类字典：B=结合、F=功能、A=ADME、T=毒性、P=理化、U=未分类",
    "assay_classification": "体内实验的疾病/治疗领域分类",
    "assay_class_map": "assays 与 assay_classification 的多对多映射",
    "assay_parameters": "实验的参数（如给药剂量、给药途径、动物品系）",
    "molecule_dictionary": "化合物主表：一个 molregno = 一个唯一化合物（含名称、研发阶段等）",
    "compound_structures": "化学结构：SMILES / InChI / InChIKey / molfile",
    "compound_properties": "计算得到的理化性质（分子量、logP、氢键供受体数等）",
    "compound_records": "文献级记录：同一化合物在不同文献里出现一次就有一行",
    "compound_structural_alerts": "结构警示：命中已知易假阳性/易毒性的子结构",
    "molecule_hierarchy": "化合物家族关系：盐/前药 → 母体活性分子",
    "molecule_synonyms": "化合物别名（商品名、代号、INN 名等）",
    "molecule_atc_classification": "化合物 → WHO ATC 药物分类的映射",
    "atc_classification": "WHO ATC 药物分类字典（按解剖-治疗-化学分级）",
    "target_dictionary": "靶点主表：一个 tid = 一个作用对象（蛋白、蛋白复合物、细胞系、整个生物体…）",
    "target_type": "靶点类型字典（SINGLE PROTEIN / PROTEIN COMPLEX / ORGANISM / CELL-LINE…）",
    "target_components": "靶点 → 组成它的分子组件（蛋白）的映射",
    "target_relations": "靶点之间的关系（如亚基属于复合物、复合物属于家族）",
    "component_sequences": "靶点组件的序列信息（UniProt accession、序列、物种）",
    "component_synonyms": "组件别名（基因名、蛋白名等）— 用基因名找靶点时从这里入手",
    "component_class": "组件 → 蛋白家族分类的映射",
    "component_go": "组件 → GO 注释的映射",
    "component_domains": "组件 → Pfam 结构域的映射",
    "protein_classification": "ChEMBL 自建的蛋白家族树（激酶/GPCR/离子通道…）",
    "binding_sites": "结合位点定义",
    "site_components": "结合位点由哪些组件构成",
    "docs": "数据来源文献/专利/数据集（PubMed ID、DOI、期刊、年份）",
    "source": "数据源字典：这条数据来自文献、PubChem、BindingDB 还是某个捐赠数据集",
    "drug_mechanism": "已知药物的作用机制：药物 × 靶点 × 作用类型（激动/抑制…）",
    "drug_indication": "药物适应症（映射到 MeSH / EFO 疾病本体）",
    "drug_warning": "药物安全性警示（黑框警告、撤市等）",
    "mechanism_refs": "drug_mechanism 的文献出处",
    "indication_refs": "drug_indication 的文献出处",
    "warning_refs": "drug_warning 的文献出处",
    "metabolism": "药物代谢途径（底物 → 代谢产物 → 催化酶）",
    "metabolism_refs": "metabolism 的文献出处",
    "products": "FDA 批准的药品（制剂层面，不是分子层面）",
    "formulations": "药品制剂 → 所含化合物的映射",
    "product_patents": "药品对应的专利信息",
    "defined_daily_dose": "WHO 定义的日剂量",
    "biotherapeutics": "生物药（抗体、多肽等）的额外信息",
    "biotherapeutic_components": "生物药 → 组成序列的映射",
    "bio_component_sequences": "生物药组分的序列",
    "cell_dictionary": "细胞系字典",
    "tissue_dictionary": "组织字典（对齐 UBERON 本体）",
    "organism_class": "物种分类字典（对齐 NCBI taxonomy）",
    "variant_sequences": "突变体序列（如激酶耐药突变）",
    "chembl_id_lookup": "全局 ChEMBL ID 索引：给一个 CHEMBLxxxx，告诉你它是化合物/靶点/文献/实验",
    "chembl_release": "本数据库的版本与发布日期",
    "version": "版本信息（旧字段）",
    "confidence_score_lookup": "confidence_score 0–9 的官方含义字典（靶点粒度 + 指认方式）",
    "relationship_type": "assay 与靶点关系类型字典",
    "data_validity_lookup": "数据可疑标记的含义字典",
    "action_type": "作用类型字典（激动剂、抑制剂…）及其上位归类",
    "bioassay_ontology": "BioAssay Ontology (BAO) 术语字典",
    "structural_alerts": "结构警示子结构定义",
    "structural_alert_sets": "结构警示规则集（PAINS、Dundee 等）",
    "ligand_eff": "配体效率指标（活性按分子大小归一化）",
    "predicted_binding_domains": "预测的结合结构域",
    "domains": "Pfam 结构域字典",
    "go_classification": "GO 分类字典",
    "usan_stems": "USAN 药物命名词干（如 -tinib 表示激酶抑制剂）",
    "patent_use_codes": "专利用途代码字典",
    "pesticide_classification": "农药分类字典",
    "pesticide_class_mapping": "化合物 → 农药分类映射",
    "sqlite_stat1": "SQLite 内部统计表（查询优化器用，不是 ChEMBL 数据）",
}

# 常见 standard_type（测量指标）的解释
ACT_TYPE_DESC = {
    "IC50": "半数抑制浓度：让目标活性下降 50% 所需的化合物浓度。数值越小 = 越强",
    "EC50": "半数效应浓度：产生 50% 最大效应所需浓度。数值越小 = 越强",
    "Ki": "抑制常数：抑制剂与靶点的结合亲和力（热力学量）。越小 = 结合越紧",
    "Kd": "解离常数：配体-靶点结合亲和力。越小 = 结合越紧",
    "AC50": "半数活性浓度（IC50/EC50 的中性叫法，常见于高通量筛选）",
    "Potency": "引发特定效应所需的浓度或剂量；PubChem 高通量筛选数据大量使用这个名字",
    "GI50": "抑制 50% 细胞生长所需浓度（细胞水平）",
    "CC50": "细胞毒性浓度：杀死 50% 细胞所需浓度",
    "TC50": "半数毒性浓度",
    "LD50": "半数致死剂量（整体动物毒性）",
    "MIC": "最低抑菌浓度：抑制细菌生长的最低浓度（抗菌实验）",
    "Inhibition": "在某个固定浓度下的抑制百分比（%），不是浓度值",
    "Activity": "笼统的“活性”，单位/含义随实验而异，通常需要看 assay 描述",
    "Percent Effect": "在固定浓度下的效应百分比",
    "Ratio": "两个测量值的比值",
    "Km": "米氏常数（酶动力学）",
    "Kb": "平衡结合常数（官方定义 equilibrium binding constant，不是拮抗剂解离常数）",
    "T1/2": "半衰期：体内或体外浓度下降一半所需时间（ADME 指标）",
    "Cmax": "峰浓度：给药后达到的最高浓度（可以是血浆，也可以是组织）",
    "AUC": "浓度-时间曲线下面积，反映总暴露量（血浆或组织）",
    "CL": "清除率：单位时间被清除掉的表观体积；也用于体外固有清除率",
    "Papp": "表观渗透系数：分子跨细胞膜（如 Caco-2 单层）的能力",
    "Solubility": "溶解度",
    "LogP": "脂水分配系数：分子的亲脂性（正辛醇/水）",
    "LogD": "在特定 pH 下的分配系数",
    "Fu": "未结合分数：未与蛋白结合的游离比例（血浆或微粒体）",
    "Emax": "最大效应",
    "Selectivity ratio": "选择性比值（对靶点 A 与靶点 B 的活性之比）",
    "kon": "结合速率常数（官方标准写法是 `k_on`，单位 M-1.s-1；数据里 `kon` 是未标准化的写法）",
    "k_off": "解离速率常数：配体从靶点上脱落的快慢；k_off 越小停留时间越长",
    "Residence time": "停留时间 = 1/k_off，药物在靶点上待多久",
    "Z score": "高通量筛选的标准化打分（相对对照组的偏离程度），不是浓度",
    "GI": "生长抑制（百分比或定性）",
    "Ratio IC50": "两个 IC50 的比值，常用于表示选择性",
    "Tissue Severity Score": "组织病理学评分（毒理实验中对病变严重程度的分级）",
    "Drug uptake": "药物摄取量",
    "Stability": "稳定性（如在肝微粒体中的剩余比例）",
    "Bioavailability": "生物利用度：口服后进入血液循环的比例",
    "Vdss": "稳态分布容积（药代动力学）",
    "Survival": "存活率/存活时间（体内实验）",
    "ED50": "半数有效剂量（体内实验）",
    "IZ": "抑菌圈直径（纸片扩散法，单位 mm）",
    "Efficacy": "疗效",
    "Time": "时间类读数",
    "Kinact/Ki": "共价抑制剂的效率参数",
}

# assay_type 字母的含义
ASSAY_TYPE_NOTE = {
    "B": "Binding：直接测化合物与靶点的结合（Ki/Kd/IC50），机制最明确，做 SAR/建模首选",
    "F": "Functional：测功能后果（酶活、细胞信号、表型），靶点归属可能间接",
    "A": "ADME：吸收/分布/代谢/排泄性质（溶解度、渗透性、肝微粒体稳定性…）",
    "T": "Toxicity：毒性相关",
    "P": "Physicochemical：纯理化性质，不涉及生物体系",
    "U": "Unassigned：未分类",
}

# confidence_score 的中文注解。
# 官方英文描述一律从库内的 confidence_score_lookup 表读取（见 §7.2），
# 这里只提供补充说明，避免把官方定义写死在代码里。
CONF_NOTE = {
    9: "实验就是在这个蛋白上做的，靶点实体最明确 —— **建模常用的门槛**",
    8: "实验在同源蛋白上做的（如另一物种的直系同源物），靶点按同源关系映射过来",
    7: "直接指认到某个蛋白复合物的亚基",
    6: "指认到同源蛋白复合物的亚基",
    5: "可能对应多个蛋白，无法唯一确定是哪一个",
    4: "可能对应多个同源蛋白",
    3: "靶点是非蛋白的分子实体（如 DNA、脂质）",
    2: "靶点是亚细胞组分（如膜制备物、微粒体）",
    1: "靶点是非分子实体（如整个细胞、组织、生物体）",
    0: "默认值：靶点未知或尚未指认",
}

MAX_PHASE_NOTE = {
    "4.0": "已获批上市药物",
    "4": "已获批上市药物",
    "3.0": "处于/已完成 III 期临床",
    "3": "处于/已完成 III 期临床",
    "2.0": "II 期临床",
    "2": "II 期临床",
    "1.0": "I 期临床",
    "1": "I 期临床",
    "0.5": "临床前 / 早期临床（ChEMBL 新增的中间档）",
    "0.0": "无临床信息（绝大多数研究化合物）",
    "0": "无临床信息（绝大多数研究化合物）",
    "-1.0": "曾出现在临床相关来源但阶段未知",
    "-1": "曾出现在临床相关来源但阶段未知",
    "None": "未标注",
}

TARGET_TYPE_NOTE = {
    "SINGLE PROTEIN": "单一蛋白 — 可直接对应一个 UniProt / 基因",
    "PROTEIN COMPLEX": "蛋白复合物（多亚基共同构成作用对象）",
    "PROTEIN FAMILY": "蛋白家族（未细分到具体成员）",
    "PROTEIN COMPLEX GROUP": "亚基组成不明确的蛋白复合物（如 GABA-A 受体），不是「家族」",
    "CELL-LINE": "细胞系整体（表型筛选）",
    "ORGANISM": "整个生物体（如某种细菌、寄生虫）",
    "TISSUE": "组织",
    "UNCHECKED": "尚未指认靶点",
    "NO TARGET": "该实验本就不适用靶点概念（如阴性对照、反筛）",
    "PROTEIN-PROTEIN INTERACTION": "以蛋白-蛋白相互作用为对象（官方定义为破坏 PPI；ChEMBL 37 起也用于靶向蛋白降解的效应蛋白/靶蛋白对）",
    "NUCLEIC-ACID": "核酸",
    "SUBCELLULAR": "亚细胞组分／制备物",
    "CHIMERIC PROTEIN": "嵌合蛋白",
    "MACROMOLECULE": "大分子",
    "SMALL MOLECULE": "小分子作为作用对象（如氨基酸、糖、代谢物）",
    "PHENOTYPE": "表型",
    "ADMET": "ADMET 实验，本就不适用靶点概念（如理化性质）",
    "UNKNOWN": "靶点的分子身份未知（药理学上定义的靶点）",
    "SELECTIVITY GROUP": "一对蛋白，用于评估化合物在两者之间的选择性",
    "OLIGOSACCHARIDE": "寡糖",
    "PROTEIN NUCLEIC-ACID COMPLEX": "蛋白-核酸复合物",
    "LIPID": "脂质",
    "3D CELL CULTURE": "三维细胞培养体系（类器官等）",
    "METAL": "金属离子",
    "NON-MOLECULAR": "非分子实体",
}

MOLTYPE_NOTE = {
    "Small molecule": "小分子（绝大多数传统药物与研究化合物）",
    "Protein": "蛋白类（含酶、融合蛋白）",
    "Antibody": "抗体",
    "Oligosaccharide": "寡糖",
    "Oligonucleotide": "寡核苷酸（ASO、siRNA 等）",
    "Cell": "细胞疗法",
    "Enzyme": "酶",
    "Gene": "基因疗法",
    "Unknown": "未标注",
    "Unclassified": "未分类",
    "Antibody drug conjugate": "抗体偶联药物（抗体连上小分子毒素）",
    "Vaccine component": "疫苗组分",
}


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def fmt(n) -> str:
    """千分位；None 显示为 —。"""
    if n is None:
        return "—"
    if isinstance(n, float):
        return f"{n:,.2f}"
    return f"{n:,}"


def pct(part, whole) -> str:
    if not whole:
        return "—"
    return f"{100.0 * part / whole:.1f}%"


def human_bytes(n: int) -> str:
    """十进制单位（与 ls / 官方发布的数字口径一致）。"""
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1000 or unit == "TB":
            return f"{x:,.2f} {unit}" if unit != "B" else f"{int(x)} B"
        x /= 1000
    return f"{x:.2f} TB"


def bar(value: float, vmax: float, width: int = 26) -> str:
    if not vmax or value <= 0:
        return ""
    return "█" * int(round(width * value / vmax))


class DB:
    """只读打开 ChEMBL，并提供带容错的查询封装。"""

    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        uri = f"file:{self.path}?mode=ro"
        self.conn = sqlite3.connect(uri, uri=True)
        self.conn.text_factory = lambda b: b.decode("utf-8", "replace")
        self._tables = {r[0].lower() for r in self.q("SELECT name FROM sqlite_master WHERE type='table'")}
        self._cols: dict[str, set] = {}

    def q(self, sql: str, params=()):
        try:
            return self.conn.execute(sql, params).fetchall()
        except sqlite3.Error as e:
            log(f"  ! 查询失败（已跳过）: {e}")
            return []

    def one(self, sql: str, params=(), default=None):
        rows = self.q(sql, params)
        return rows[0][0] if rows and rows[0] else default

    def has_table(self, t: str) -> bool:
        return t.lower() in self._tables

    def cols(self, t: str) -> set:
        t = t.lower()
        if t not in self._cols:
            self._cols[t] = {r[1].lower() for r in self.q(f"PRAGMA table_info({t})")}
        return self._cols[t]

    def has_col(self, t: str, c: str) -> bool:
        return self.has_table(t) and c.lower() in self.cols(t)

    def tables(self) -> list:
        return sorted(self._tables)

    def count(self, t: str, exact: bool = True):
        """行数。exact=False 时优先用 sqlite_stat1 的估计值（瞬间返回）。"""
        if not self.has_table(t):
            return None
        if not exact and self.has_table("sqlite_stat1"):
            s = self.one("SELECT stat FROM sqlite_stat1 WHERE tbl=? LIMIT 1", (t,))
            if s:
                try:
                    return int(str(s).split()[0])
                except (ValueError, IndexError):
                    pass
        return self.one(f"SELECT COUNT(*) FROM {t}")


# ---------------------------------------------------------------------------
# 报告生成
# ---------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.lines: list[str] = []

    def __call__(self, s: str = "") -> None:
        self.lines.append(s)

    def h(self, level: int, text: str) -> None:
        self("")
        self("#" * level + " " + text)
        self("")

    def table(self, headers: list, rows: list, aligns: list | None = None) -> None:
        if not rows:
            self("*（无数据）*")
            self("")
            return
        self("| " + " | ".join(str(h) for h in headers) + " |")
        if aligns is None:
            aligns = ["---"] * len(headers)
        self("|" + "|".join(f" {a} " for a in aligns) + "|")
        for r in rows:
            cells = ["" if c is None else str(c).replace("|", "\\|").replace("\n", " ") for c in r]
            self("| " + " | ".join(cells) + " |")
        self("")

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


# ---------------------------------------------------------------------------

def sec_intro(r: Report, db: DB, args) -> dict:
    """概览 + 版本信息。返回一些后续要复用的数字。"""
    log("章节 1/10：概览与版本")
    size = os.path.getsize(db.path)

    rel_name, rel_date, rel_comment = None, None, None
    if db.has_table("chembl_release"):
        rows = db.q("SELECT chembl_release, creation_date FROM chembl_release ORDER BY chembl_release_id DESC LIMIT 1")
        if rows:
            rel_name, rel_date = rows[0][0], rows[0][1]
    if db.has_table("version"):
        # version 表存的是「本版依赖的各种外部资源版本」，其中一行才是 ChEMBL 自身
        rows = db.q("SELECT name, creation_date, comments FROM version WHERE name LIKE 'ChEMBL/_%' ESCAPE '/' LIMIT 1")
        if rows:
            rel_name = rel_name or rows[0][0]
            rel_date = rel_date or rows[0][1]
            rel_comment = rows[0][2]

    n_mol = db.count("molecule_dictionary", not args.quick)
    n_act = db.count("activities", not args.quick)
    n_assay = db.count("assays", not args.quick)
    n_tgt = db.count("target_dictionary", not args.quick)
    n_doc = db.count("docs", not args.quick)

    r("# ChEMBL 数据库结构与内容概览报告")
    r("")
    r(f"> 自动生成于 {datetime.now().strftime('%Y-%m-%d %H:%M')}　|　"
      f"数据库文件：`{db.path}`")
    r(">")
    r("> 本报告面向**有生物学 / 生信背景、但没有药物化学背景**的读者。"
      "每一节先解释「这是什么」，再给出这个数据库里的实际统计。")
    r("")

    r.h(2, "0. 三十秒速览")
    r("")
    r("**ChEMBL 是什么？** 一个由 EMBL-EBI 维护的人工整理（manually curated）的生物活性数据库。"
      "它把散落在几万篇论文和专利里的一句句话——「化合物 X 抑制蛋白 Y，IC50 = 12 nM」"
      "——抽取成结构化的表格。你可以把它理解成**化合物-靶点相互作用领域的 UniProt + GEO**：")
    r("")
    r("- 像 UniProt 一样，它有稳定 ID（`CHEMBL25` 就是阿司匹林）、人工审编、版本化发布；")
    r("- 像 GEO 一样，它的一条记录必须绑定「在什么实验条件下测的」，脱离实验条件的数值没有意义。")
    r("")
    r("**一句话的数据模型**：")
    r("")
    r("```")
    r("某个化合物  在某个实验里  针对某个靶点  测到了某个数值  出处是某篇文献")
    r("(compound)    (assay)       (target)      (activity)      (document)")
    r("```")
    r("")
    r("这五个实体就是 ChEMBL 的骨架，其余 60 多张表都是围绕它们的注释、字典和映射。")
    r("")

    r.table(
        ["条目", "值"],
        [
            ["数据库版本", rel_name or "（未标注）"],
            ["版本创建日期", rel_date or "（未标注）"],
            ["版本说明", rel_comment or "—"],
            ["文件大小", human_bytes(size)],
            ["表数量", fmt(len(db.tables()))],
            ["化合物数（molecule_dictionary）", fmt(n_mol)],
            ["活性数据点数（activities）", fmt(n_act)],
            ["实验数（assays）", fmt(n_assay)],
            ["靶点数（target_dictionary）", fmt(n_tgt)],
            ["来源文献数（docs）", fmt(n_doc)],
            ["SQLite 引擎（运行本脚本的）", sqlite3.sqlite_version],
        ],
    )
    if args.quick:
        r("> ⚠️ `--quick` 模式：以上行数为 SQLite 统计表给出的**估计值**，可能与精确值有出入。")
        r("")

    # 外部资源版本 —— 对可重复性很重要
    if db.has_table("version"):
        rows = db.q("SELECT name, comments FROM version ORDER BY name")
        rows = [x for x in rows if not str(x[0]).lower().startswith("chembl_3")]
        if rows:
            r.h(3, "0.1 本版所依赖的外部资源版本")
            r("")
            r("ChEMBL 的注释是**站在其他数据库肩膀上**的。下表来自 `version` 表，"
              "记录了本次发布所用的外部资源版本。做可重复性分析或跨库整合时，"
              "这决定了你的 UniProt / GO / MeSH 注释应该对齐到哪个版本。")
            r("")
            r.table(["资源", "用途"], [[f"**{x[0]}**", x[1]] for x in rows])

    # 发布历史
    if db.has_table("chembl_release"):
        rows = db.q("SELECT chembl_release, creation_date FROM chembl_release ORDER BY chembl_release_id DESC LIMIT 8")
        if rows:
            r.h(3, "0.2 最近的发布历史")
            r("")
            r("> ChEMBL 大约每 6–12 个月发一版。**每一版都会回溯性地修订旧数据**"
              "（合并重复化合物、重新指认靶点、调整临床阶段标注），"
              "所以任何分析都必须写明用的是哪一版。")
            r("")
            r.table(["版本", "创建日期"], [[x[0], str(x[1])[:10]] for x in rows])

    # ChEMBL ID 分布
    if db.has_table("chembl_id_lookup") and not args.quick:
        rows = db.q("SELECT entity_type, COUNT(*) c FROM chembl_id_lookup GROUP BY 1 ORDER BY c DESC")
        if rows:
            r.h(3, "0.3 所有 CHEMBLxxxxx 标识符的构成")
            r("")
            r("拿到一个陌生的 `CHEMBLxxxxx`，可以在 `chembl_id_lookup` 里查它到底是什么实体。")
            r("")
            r.table(["entity_type", "数量", "是什么"],
                    [[f"`{x[0]}`", fmt(x[1]),
                      {"COMPOUND": "化合物", "ASSAY": "实验", "TARGET": "靶点",
                       "DOCUMENT": "文献/数据集", "CELL": "细胞系", "TISSUE": "组织"}.get(x[0], "")]
                     for x in rows],
                    aligns=["---", "---:", "---"])

    return {"n_mol": n_mol, "n_act": n_act, "n_assay": n_assay, "n_tgt": n_tgt, "n_doc": n_doc}


def sec_glossary(r: Report, db: DB) -> None:
    log("章节 2/10：概念词表")
    r.h(2, "1. 先搞懂这些词")
    r("")
    r("如果你是从生物学过来的，下面这几个概念是读懂 ChEMBL 的全部前提。")
    r("")

    r.h(3, "1.1 化合物（compound / molecule）")
    r("")
    r("一个小分子药物或研究用化合物。它在数据库里的身份证有两套：")
    r("")
    r("- **`molregno`** — 内部整数主键，所有表之间用它连接（相当于 NCBI 的 GeneID）；")
    r("- **`chembl_id`** — 对外稳定 ID，形如 `CHEMBL25`（相当于基因的 Ensembl ID）。")
    r("")
    r("化学结构本身以文本形式存放，最常用的是 **SMILES**：一行 ASCII 字符串就编码了一个分子的"
      "原子连接关系，例如阿司匹林是 `CC(=O)Oc1ccccc1C(=O)O`。你不需要会读 SMILES —— "
      "RDKit 之类的库可以把它解析成分子对象、算描述符、算相似度。另外还有 **InChIKey**，"
      "一个 27 位的定长哈希（如 `BSYNRYMUTXBXSQ-UHFFFAOYSA-N`），"
      "作用完全等同于序列的 MD5：用来跨数据库精确匹配同一个分子。")
    r("")

    r.h(3, "1.2 靶点（target）")
    r("")
    r("化合物作用的对象。**注意它不总是一个蛋白**：可能是单一蛋白、蛋白复合物、"
      "一条蛋白家族、一个细胞系，甚至一整个物种（抗菌实验里靶点就是「大肠杆菌」）。"
      "`target_dictionary.target_type` 记录了是哪一种。")
    r("")
    r("要从**基因名**找到靶点，路径是："
      "`component_synonyms`（基因名/蛋白别名）→ `component_sequences`（UniProt accession）"
      "→ `target_components` → `target_dictionary`。"
      "`component_sequences.accession` 就是 UniProt 号，这是 ChEMBL 与你熟悉的生物学世界之间最重要的一座桥。")
    r("")

    r.h(3, "1.3 实验（assay）")
    r("")
    r("一次具体的测量设定：用什么体系、测什么读数、针对哪个靶点。"
      "**同一个化合物-靶点对，在不同实验里测出的数值可能差一个数量级**，这不是数据错误，"
      "而是实验条件不同（ATP 浓度、细胞类型、孵育时间…）。所以 ChEMBL 从不把活性直接挂在"
      "「化合物-靶点」上，而是必须经过 assay 这一层。这一点和 GEO 里样本必须绑定实验设计是一个道理。")
    r("")
    r("两个字段决定了这条数据能不能用：")
    r("")
    r("- **`assay_type`** — 实验大类（见 §5）；")
    r("- **`confidence_score`（0–9）** — 描述这个实验被归因到了**什么粒度的靶点**、"
      "以及归因是直接的还是靠同源关系映射的。它**不是数据质量分**（详见 §7.2）。"
      "做定量建模一般要求 **≥ 8**，即靶点是一个单一蛋白。")
    r("")

    r.h(3, "1.4 活性（activity）——最容易踩坑的一张表")
    r("")
    r("`activities` 是核心事实表，一行 = 一个数值。关键在于它有**两套值**：")
    r("")
    r.table(
        ["字段", "含义"],
        [
            ["`type` / `value` / `units` / `relation`", "**原文照抄**：论文里怎么写就怎么存（可能是 `Ic-50`、单位 µM）"],
            ["`standard_type` / `standard_value` / `standard_units` / `standard_relation`", "**标准化后**：类型名统一（`IC50`）、单位统一（浓度一律 nM）"],
        ],
    )
    r("**永远用 `standard_*` 这一套做分析**，`type/value` 只在追溯原文时用。")
    r("")
    r("`standard_relation` 是个容易被忽略的坑：它可能是 `=`，也可能是 `>` 或 `<`。"
      "`> 10000 nM` 的意思是「测到最高浓度也没打到 50%，实际值未知，只知道比这大」——"
      "这是**删失数据（censored data）**，直接当成 10000 会污染回归模型。生存分析里的删失是同一个概念。")
    r("")

    r.h(3, "1.5 pChEMBL —— 你真正应该用的那个数")
    r("")
    r("`pchembl_value = −log10(标准化浓度，单位 M)`，只对浓度型指标"
      "（IC50/EC50/XC50/AC50/Ki/Kd/Potency）且 `relation = '='` 时计算。")
    r("")
    r("- IC50 = 1 nM  → pChEMBL = 9")
    r("- IC50 = 100 nM → pChEMBL = 7")
    r("- IC50 = 10 µM  → pChEMBL = 5")
    r("")
    r("为什么要取负对数？因为**结合自由能与浓度的对数成正比**（ΔG = −RT·ln K），"
      "所以 pChEMBL 才是那个物理上线性、分布近似正态、适合做回归的量。"
      "这跟你对表达量取 log2 是完全一样的动机。**数值越大 = 活性越强**（方向和 IC50 相反，注意别搞反）。")
    r("")
    r("经验尺度：pChEMBL ≥ 6（IC50 ≤ 1 µM）算「有活性」，≥ 8（≤ 10 nM）算「强效」。")
    r("")

    r.h(3, "1.6 max_phase —— 这个化合物走到哪一步了")
    r("")
    r("`molecule_dictionary.max_phase` 记录化合物达到过的最高临床阶段："
      "4 = 已上市药物，3/2/1 = 对应临床期，0.5 = 早期，0 或空 = 纯研究化合物。"
      "想快速捞出「所有已知药物」，条件就是 `max_phase = 4`。")
    r("")
    r("> 注意：不同 ChEMBL 版本之间 `max_phase` 会变动（官方在持续复核「什么算获批」），"
      "跨版本比较药物数量时务必注明版本号。")
    r("")


def sec_schema(r: Report, db: DB, args) -> None:
    log("章节 3/10：数据模型与连接关系")
    r.h(2, "2. 数据模型：表是怎么连起来的")
    r("")
    r("```")
    r("                    compound_properties   compound_structures")
    r("                          │ molregno            │ molregno")
    r("                          └──────┬──────────────┘")
    r("                                 │")
    r("                        molecule_dictionary          docs ──────┐")
    r("                                 │ molregno       (文献/专利)   │ doc_id")
    r("                                 │                              │")
    r("        compound_records ────────┤                              │")
    r("         (文献级记录)  record_id │                              │")
    r("                                 │                              │")
    r("                          ┌──────┴──────────────────────────────┴──┐")
    r("                          │            ACTIVITIES                  │  ← 核心事实表")
    r("                          │  molregno / assay_id / doc_id / src_id │")
    r("                          └──────────────────┬─────────────────────┘")
    r("                                             │ assay_id")
    r("                                          ASSAYS")
    r("                                     (实验设定 + 可信度)")
    r("                                             │ tid")
    r("                                     target_dictionary")
    r("                                             │ tid")
    r("                                      target_components")
    r("                                             │ component_id")
    r("                              component_sequences  (UniProt accession)")
    r("                                             │")
    r("                        component_synonyms / component_go / component_class")
    r("```")
    r("")
    r("**读法**：从 `activities` 出发，向左连到「测的是什么分子」，向下连到「怎么测的、测的什么靶点」，"
      "向上连到「出处是哪篇文献」。几乎所有实用查询都是这个骨架的变形。")
    r("")

    # 从 PRAGMA 抽外键，展示核心表的实际连接
    r.h(3, "2.1 核心表的外键（从数据库实际读取）")
    r("")
    core = ["activities", "assays", "target_dictionary", "target_components",
            "molecule_dictionary", "compound_records", "compound_structures"]
    rows = []
    for t in core:
        if not db.has_table(t):
            continue
        for fk in db.q(f"PRAGMA foreign_key_list({t})"):
            # (id, seq, table, from, to, on_update, on_delete, match)
            rows.append([f"`{t}`", f"`{fk[3]}`", f"`{fk[2]}`", f"`{fk[4]}`"])
    r.table(["表", "外键列", "指向表", "指向列"], rows)

    r.h(3, "2.2 主键命名约定")
    r("")
    r.table(
        ["主键", "属于", "说明"],
        [
            ["`molregno`", "化合物", "molecule registration number"],
            ["`tid`", "靶点", "target id"],
            ["`assay_id`", "实验", ""],
            ["`activity_id`", "活性数据点", ""],
            ["`doc_id`", "文献", ""],
            ["`record_id`", "文献级化合物记录", "同一化合物出现在 N 篇文献 → N 个 record_id，但只有 1 个 molregno"],
            ["`component_id`", "靶点组件（蛋白）", ""],
            ["`src_id`", "数据来源", ""],
        ],
    )
    r("> `chembl_id_lookup` 表可以反查任意 `CHEMBLxxxxx` 属于哪一类实体——"
      "拿到一个陌生 ChEMBL ID 时先查它。")
    r("")


def sec_tables(r: Report, db: DB, args) -> None:
    log("章节 4/10：表清单与行数（可能较慢）")
    r.h(2, "3. 表清单与规模")
    r("")
    exact = not args.quick
    rows = []
    for t in db.tables():
        n = db.count(t, exact)
        rows.append((t, n if n is not None else -1))
    rows.sort(key=lambda x: -x[1])

    r(f"共 **{len(rows)}** 张表，按行数排序（{'精确计数' if exact else '估计值'}）。"
      "「说明」一栏是本脚本内置的中文注解。")
    r("")
    r.table(
        ["表名", "行数", "说明"],
        [[f"`{t}`", fmt(n if n >= 0 else None), TABLE_DESC.get(t, "")] for t, n in rows],
        aligns=["---", "---:", "---"],
    )


def sec_sources(r: Report, db: DB, args) -> None:
    log("章节 5/10：数据来源")
    r.h(2, "4. 数据从哪来")
    r("")
    r("ChEMBL 不是单一来源。`source` 表列出了所有数据源，`activities.src_id` 指向它。"
      "了解来源很重要，因为**不同来源的数据质量与稠密度差别很大**："
      "文献数据经过人工审编但稀疏，高通量筛选数据量大但多为单浓度、阴性居多。")
    r("")
    if not db.has_table("source"):
        r("*（本库无 source 表）*")
        return

    if db.has_col("activities", "src_id") and not args.quick:
        rows = db.q("""
            SELECT s.src_id, s.src_short_name, s.src_description, COUNT(a.activity_id)
            FROM source s LEFT JOIN activities a ON a.src_id = s.src_id
            GROUP BY s.src_id, s.src_short_name, s.src_description
            ORDER BY COUNT(a.activity_id) DESC
        """)
        vmax = max([x[3] for x in rows], default=0)
        r.table(
            ["src_id", "简称", "描述", "activities 数", ""],
            [[x[0], f"`{x[1]}`", x[2], fmt(x[3]), bar(x[3], vmax, 20)] for x in rows],
            aligns=["---:", "---", "---", "---:", "---"],
        )
    else:
        rows = db.q("SELECT src_id, src_short_name, src_description FROM source ORDER BY src_id")
        r.table(["src_id", "简称", "描述"], [[x[0], f"`{x[1]}`", x[2]] for x in rows])


def sec_compounds(r: Report, db: DB, args) -> None:
    log("章节 6/10：化合物")
    r.h(2, "5. 化合物层面")
    r("")
    n_mol = db.count("molecule_dictionary", not args.quick)
    n_struct = db.count("compound_structures", not args.quick)
    n_prop = db.count("compound_properties", not args.quick)
    r(f"- 唯一化合物：**{fmt(n_mol)}**")
    r(f"- 其中有化学结构（SMILES 等）：**{fmt(n_struct)}**"
      f"（{pct(n_struct or 0, n_mol or 0)}）—— 没有结构的多为抗体、细胞疗法等大分子")
    r(f"- 有计算理化性质：**{fmt(n_prop)}**（{pct(n_prop or 0, n_mol or 0)}）")
    r("")

    if db.has_col("molecule_dictionary", "molecule_type"):
        rows = db.q("SELECT COALESCE(molecule_type,'(空)'), COUNT(*) c FROM molecule_dictionary GROUP BY 1 ORDER BY c DESC")
        vmax = max([x[1] for x in rows], default=0)
        r.h(3, "5.1 化合物类型")
        r.table(
            ["molecule_type", "数量", "占比", "", "说明"],
            [[x[0], fmt(x[1]), pct(x[1], n_mol or 0), bar(x[1], vmax, 18), MOLTYPE_NOTE.get(x[0], "")] for x in rows],
            aligns=["---", "---:", "---:", "---", "---"],
        )

    if db.has_col("molecule_dictionary", "max_phase"):
        rows = db.q("SELECT COALESCE(CAST(max_phase AS TEXT),'None'), COUNT(*) c FROM molecule_dictionary GROUP BY 1 ORDER BY c DESC")
        r.h(3, "5.2 研发阶段 max_phase")
        r.table(
            ["max_phase", "化合物数", "含义"],
            [[x[0], fmt(x[1]), MAX_PHASE_NOTE.get(str(x[0]), "")] for x in rows],
            aligns=["---", "---:", "---"],
        )
        n4 = db.one("SELECT COUNT(*) FROM molecule_dictionary WHERE max_phase = 4")
        r(f"> 想拿「已上市药物」集合：`SELECT * FROM molecule_dictionary WHERE max_phase = 4`，"
          f"本库共 **{fmt(n4)}** 个。")
        r("")

    # 若干标记位
    flags = [
        ("therapeutic_flag", "被标记为治疗用药"),
        ("natural_product", "天然产物来源"),
        ("oral", "可口服"),
        ("parenteral", "注射给药"),
        ("topical", "外用"),
        ("black_box_warning", "有黑框警告（严重安全性问题）"),
        ("withdrawn_flag", "已撤市"),
        ("prodrug", "前药（体内代谢后才有活性）"),
        ("chemical_probe", "化学探针（工具分子，选择性经过验证）"),
        ("inorganic_flag", "无机物"),
        ("polymer_flag", "聚合物"),
    ]
    rows = []
    for col, note in flags:
        if db.has_col("molecule_dictionary", col):
            n = db.one(f"SELECT COUNT(*) FROM molecule_dictionary WHERE {col} = 1")
            rows.append([f"`{col}`", fmt(n), note])
    if rows:
        r.h(3, "5.3 化合物标记位（值为 1 的数量）")
        r.table(["字段", "数量", "含义"], rows, aligns=["---", "---:", "---"])

    # 理化性质
    if db.has_table("compound_properties"):
        r.h(3, "5.4 计算理化性质的分布")
        r("")
        r("这些是**用软件从结构算出来的**（不是实验测的），常用于快速过滤化合物。"
          "对非化学背景读者，最需要知道的是所谓 **Lipinski 类药五规则（Rule of Five）**："
          "分子量 ≤ 500、logP ≤ 5、氢键供体 ≤ 5、氢键受体 ≤ 10 —— "
          "满足这些的分子更可能有良好口服吸收。`num_ro5_violations` 就是违反了几条。")
        r("")
        prop_notes = {
            "mw_freebase": ("分子量（游离碱形式）", "Da"),
            "alogp": ("计算 logP，亲脂性；越大越亲脂", ""),
            "hba": ("氢键受体数", ""),
            "hbd": ("氢键供体数", ""),
            "psa": ("极性表面积；与膜通透性、口服吸收相关", "Å²"),
            "rtb": ("可旋转键数；反映分子柔性", ""),
            "num_ro5_violations": ("违反 Lipinski 五规则的条数", ""),
            "aromatic_rings": ("芳香环数", ""),
            "heavy_atoms": ("重原子（非氢原子）数", ""),
            "qed_weighted": ("QED 类药性综合评分，0–1，越大越「像药」", ""),
            "cx_logd": ("pH 7.4 下的分配系数", ""),
            "np_likeness_score": ("天然产物相似度评分", ""),
        }
        rows = []
        for col, (note, unit) in prop_notes.items():
            if not db.has_col("compound_properties", col):
                continue
            res = db.q(f"SELECT COUNT({col}), MIN({col}), AVG({col}), MAX({col}) FROM compound_properties WHERE {col} IS NOT NULL")
            if not res or not res[0][0]:
                continue
            cnt, mn, avg, mx = res[0]
            med = db.one(f"SELECT {col} FROM compound_properties WHERE {col} IS NOT NULL ORDER BY {col} LIMIT 1 OFFSET {cnt // 2}")
            rows.append([f"`{col}`", fmt(cnt), f"{mn:,.2f}", f"{med:,.2f}" if med is not None else "—",
                         f"{avg:,.2f}", f"{mx:,.2f}", unit, note])
        r.table(["字段", "有值数", "最小", "中位", "均值", "最大", "单位", "含义"], rows,
                aligns=["---", "---:", "---:", "---:", "---:", "---:", "---", "---"])


def sec_targets(r: Report, db: DB, args) -> None:
    log("章节 7/10：靶点")
    r.h(2, "6. 靶点层面")
    r("")
    n_tgt = db.count("target_dictionary", not args.quick)
    r(f"靶点总数：**{fmt(n_tgt)}**。再强调一次：靶点不等于蛋白，见下表分布。")
    r("")

    if db.has_col("target_dictionary", "target_type"):
        rows = db.q("SELECT COALESCE(target_type,'(空)'), COUNT(*) c FROM target_dictionary GROUP BY 1 ORDER BY c DESC")
        vmax = max([x[1] for x in rows], default=0)
        tt_off = {x[0]: x[1] for x in db.q("SELECT target_type, target_desc FROM target_type")}
        r.h(3, "6.1 靶点类型")
        r("")
        r("「官方描述」一列读自库内的 `target_type` 字典表，中文一列是本脚本的补充注解。")
        r("")
        r.table(
            ["target_type", "数量", "", "官方描述（原文）", "中文说明"],
            [[x[0], fmt(x[1]), bar(x[1], vmax, 12), tt_off.get(x[0], "—"),
              TARGET_TYPE_NOTE.get(x[0], "")] for x in rows],
            aligns=["---", "---:", "---", "---", "---"],
        )

    if db.has_col("target_dictionary", "organism"):
        rows = db.q(f"SELECT COALESCE(organism,'(空)'), COUNT(*) c FROM target_dictionary GROUP BY 1 ORDER BY c DESC LIMIT {args.top}")
        r.h(3, f"6.2 靶点物种分布（Top {args.top}）")
        r.table(["物种", "靶点数"], [[x[0], fmt(x[1])] for x in rows], aligns=["---", "---:"])

    # 人源单一蛋白靶点 —— 生信读者最关心的子集
    if db.has_table("target_dictionary"):
        n_human_sp = db.one(
            "SELECT COUNT(*) FROM target_dictionary WHERE target_type='SINGLE PROTEIN' AND organism='Homo sapiens'")
        r(f"> **人源单一蛋白靶点**共 **{fmt(n_human_sp)}** 个 —— 这通常是做人类靶点分析时的起点子集。")
        r("")

    # 靶点 → UniProt 桥梁
    if db.has_table("component_sequences"):
        r.h(3, "6.3 与 UniProt / 基因的对接")
        n_comp = db.count("component_sequences", not args.quick)
        n_acc = db.one("SELECT COUNT(DISTINCT accession) FROM component_sequences WHERE accession IS NOT NULL")
        r(f"- `component_sequences` 共 {fmt(n_comp)} 行，其中 **{fmt(n_acc)} 个不同的 UniProt accession**。")
        if db.has_col("component_sequences", "component_type"):
            rows = db.q("SELECT COALESCE(component_type,'(空)'), COUNT(*) c FROM component_sequences GROUP BY 1 ORDER BY c DESC")
            r("- 组件类型：" + "，".join(f"{x[0]} {fmt(x[1])}" for x in rows))
        r("")
        r("从基因名查靶点的标准写法：")
        r("")
        r("```sql")
        r("SELECT DISTINCT td.chembl_id, td.pref_name, td.organism, td.target_type")
        r("FROM component_synonyms cs")
        r("JOIN component_sequences seq ON cs.component_id = seq.component_id")
        r("JOIN target_components tc    ON tc.component_id = seq.component_id")
        r("JOIN target_dictionary td    ON td.tid = tc.tid")
        r("WHERE cs.component_synonym = 'GCK'        -- 基因名")
        r("  AND cs.syn_type = 'GENE_SYMBOL';")
        r("```")
        r("")

    # 蛋白家族分类
    if db.has_table("protein_classification") and db.has_col("protein_classification", "pref_name"):
        rows = db.q("""
            SELECT pc.pref_name, COUNT(DISTINCT cc.component_id) c
            FROM protein_classification pc
            JOIN component_class cc ON cc.protein_class_id = pc.protein_class_id
            GROUP BY pc.pref_name ORDER BY c DESC LIMIT 15
        """)
        if rows:
            r.h(3, "6.4 蛋白家族分类（Top 15，按组件数）")
            r.table(["蛋白家族", "组件数"], [[x[0], fmt(x[1])] for x in rows], aligns=["---", "---:"])


def sec_assays(r: Report, db: DB, args) -> None:
    log("章节 8/10：实验")
    r.h(2, "7. 实验（assay）层面")
    r("")
    n_assay = db.count("assays", not args.quick)
    r(f"实验总数：**{fmt(n_assay)}**。")
    r("")

    if db.has_col("assays", "assay_type"):
        rows = db.q("SELECT COALESCE(assay_type,'(空)'), COUNT(*) c FROM assays GROUP BY 1 ORDER BY c DESC")
        vmax = max([x[1] for x in rows], default=0)
        r.h(3, "7.1 实验类型 assay_type")
        r.table(
            ["类型", "实验数", "", "含义"],
            [[x[0], fmt(x[1]), bar(x[1], vmax, 18), ASSAY_TYPE_NOTE.get(x[0], "")] for x in rows],
            aligns=["---", "---:", "---", "---"],
        )

    if db.has_col("assays", "confidence_score"):
        rows = db.q("SELECT COALESCE(confidence_score,-1), COUNT(*) c FROM assays GROUP BY 1 ORDER BY 1 DESC")
        vmax = max([x[1] for x in rows], default=0)
        # 官方定义直接从库内字典表读，不写死在代码里
        lookup = {x[0]: (x[1], x[2]) for x in
                  db.q("SELECT confidence_score, description, target_mapping FROM confidence_score_lookup")}
        r.h(3, "7.2 confidence_score：靶点指认的类型与直接程度 ⭐")
        r("")
        r("**这是筛选数据时最重要的字段之一，也是最容易被误解的一个。**")
        r("")
        r("它常被叫作「可信度分数」，但它**不是一个单纯的数据质量分**。"
          "它同时编码了两件事：")
        r("")
        r("1. **靶点实体有多确定** —— 单一蛋白 > 蛋白复合物 > 多个候选蛋白 > 非蛋白分子 > 亚细胞组分 > 非分子实体；")
        r("2. **指认是直接的还是靠同源关系映射的** —— 这正是 9 与 8、7 与 6、5 与 4 之间的区别："
          "奇数档是 direct（实验就在该靶点上做的），偶数档是 homologous"
          "（实验在同源物上做的，例如用大鼠蛋白测的活性被映射到人的直系同源蛋白）。")
        r("")
        r("所以低分**不等于**数据质量差 —— 一个 `confidence_score = 1` 的抗菌实验可能做得非常严谨，"
          "只是它的靶点是「整个细菌」这种非分子实体。这个字段回答的是"
          "「**这条数据能归因到什么粒度的靶点**」，而不是「这条数据可不可信」。")
        r("")
        r("下表的「官方描述」一列直接读自本库的 `confidence_score_lookup` 表：")
        r("")
        r.table(
            ["分数", "实验数", "", "官方描述（原文）", "target_mapping", "中文说明"],
            [[x[0] if x[0] != -1 else "(空)", fmt(x[1]), bar(x[1], vmax, 12),
              lookup.get(x[0], ("—", "—"))[0], lookup.get(x[0], ("—", "—"))[1],
              CONF_NOTE.get(x[0], "")] for x in rows],
            aligns=["---:", "---:", "---", "---", "---", "---"],
        )
        n_hi = db.one("SELECT COUNT(*) FROM assays WHERE confidence_score >= 8")
        r(f"> **常用过滤条件 `confidence_score >= 8`**（本库 **{fmt(n_hi)}** 个实验，"
          f"{pct(n_hi or 0, n_assay or 0)}）的真实含义是"
          "「靶点是一个单一蛋白，无论直接指认还是同源映射」，而不是「数据质量达到 8 分」。"
          "做 QSAR / 机器学习建模时这几乎是标配，因为建模需要每条数据对应一个明确的蛋白。"
          "如果你连同源映射也不想要（例如只认人源蛋白上实测的数据），就用 `= 9`。")
        r("")

    if db.has_col("assays", "assay_organism"):
        rows = db.q(f"SELECT COALESCE(assay_organism,'(空)'), COUNT(*) c FROM assays GROUP BY 1 ORDER BY c DESC LIMIT {args.top}")
        r.h(3, f"7.3 实验体系物种（Top {args.top}）")
        r("")
        r("> 注意区分 `assay_organism`（实验体系来自哪个物种，例如用大鼠肝微粒体）与 "
          "`target_dictionary.organism`（靶点蛋白本身的物种）—— 两者可以不同。")
        r("")
        r.table(["assay_organism", "实验数"], [[x[0], fmt(x[1])] for x in rows], aligns=["---", "---:"])

    # 每个靶点有多少实验（快，assays 表内即可）
    if db.has_col("assays", "tid") and db.has_table("target_dictionary"):
        rows = db.q(f"""
            SELECT td.chembl_id, td.pref_name, td.organism, COUNT(*) c
            FROM assays a JOIN target_dictionary td ON a.tid = td.tid
            GROUP BY td.chembl_id, td.pref_name, td.organism
            ORDER BY c DESC LIMIT {args.top}
        """)
        r.h(3, f"7.4 实验数最多的靶点（Top {args.top}）")
        r.table(["ChEMBL ID", "靶点名", "物种", "实验数"],
                [[f"`{x[0]}`", x[1], x[2], fmt(x[3])] for x in rows],
                aligns=["---", "---", "---", "---:"])


def sec_activities(r: Report, db: DB, args) -> None:
    log("章节 9/10：活性数据（全表扫描，较慢）")
    r.h(2, "8. 活性数据层面 ⭐")
    r("")
    n_act = db.count("activities", not args.quick)
    r(f"活性数据点总数：**{fmt(n_act)}**。这是整个数据库的重心，也是最需要小心处理的部分。")
    r("")

    if args.quick:
        r("> ⚠️ `--quick` 模式跳过了本节的全表统计。去掉该参数可获得完整分布。")
        r("")
        return

    # 单次扫描算出所有计数
    has_mod = db.has_col("activities", "modality")
    sel = [
        "COUNT(*)",
        "SUM(CASE WHEN standard_value IS NOT NULL THEN 1 ELSE 0 END)",
        "SUM(CASE WHEN pchembl_value IS NOT NULL THEN 1 ELSE 0 END)",
        "SUM(CASE WHEN standard_relation = '=' THEN 1 ELSE 0 END)",
        "SUM(CASE WHEN data_validity_comment IS NOT NULL THEN 1 ELSE 0 END)",
        "SUM(CASE WHEN potential_duplicate = 1 THEN 1 ELSE 0 END)",
        "SUM(CASE WHEN standard_flag = 1 THEN 1 ELSE 0 END)",
    ]
    if has_mod:
        sel.append("SUM(CASE WHEN modality IS NOT NULL THEN 1 ELSE 0 END)")
    log("  · 单次扫描统计各标记位…")
    res = db.q(f"SELECT {', '.join(sel)} FROM activities")
    if res:
        v = list(res[0])
        total = v[0]
        rows = [
            ["有标准化数值 `standard_value`", fmt(v[1]), pct(v[1], total), "没有数值的多为定性结论（Active/Inactive）"],
            ["有 `pchembl_value`", fmt(v[2]), pct(v[2], total), "**做定量分析的可用子集**"],
            ["`standard_relation = '='`", fmt(v[3]), pct(v[3], total), "精确值；其余为 > / < 的删失数据"],
            ["有 `data_validity_comment` 标记", fmt(v[4]), pct(v[4], total), "多数是可疑标记（建议剔除），但 `Manually validated` 是正面标记，见 §8.6"],
            ["疑似重复引用 `potential_duplicate = 1`", fmt(v[5]), pct(v[5], total), "同一数值被多篇文献转述，建议剔除"],
            ["已人工标准化 `standard_flag = 1`", fmt(v[6]), pct(v[6], total), ""],
        ]
        if has_mod:
            rows.append(["有 `modality` 标注", fmt(v[7]), pct(v[7], total), "ChEMBL 37 起新增，目前主要标注靶向蛋白降解"])
        r.h(3, "8.1 数据完整性与质量标记")
        r.table(["指标", "数量", "占比", "备注"], rows, aligns=["---", "---:", "---:", "---"])

    # standard_type 分布
    log("  · 统计 standard_type 分布…")
    rows = db.q(f"SELECT COALESCE(standard_type,'(空)'), COUNT(*) c FROM activities GROUP BY 1 ORDER BY c DESC LIMIT {args.top}")
    if rows:
        vmax = max(x[1] for x in rows)
        # 官方定义读自 activity_stds_lookup（同一 type 可能有多行，只是允许的单位不同）
        std_off: dict[str, str] = {}
        std_units: dict[str, set] = {}
        for st, defn, un in db.q("SELECT standard_type, definition, standard_units FROM activity_stds_lookup"):
            std_off[st] = defn
            std_units.setdefault(st, set()).add(un)

        r.h(3, f"8.2 测量指标 standard_type（Top {args.top}）")
        r("")
        r("每一行是一种「测的是什么量」。**「官方定义」一列读自库内的 `activity_stds_lookup` 表**，"
          "中文一列是本脚本的补充注解。")
        r("")
        r("这张表还有个额外用处：**能不能在 `activity_stds_lookup` 里查到，本身就是一个质量信号**。"
          "查得到，说明 ChEMBL 为它定义了标准单位和合理取值范围，会做单位换算与越界检查；"
          "查不到（下表中官方定义为「—」的那些，如 `Activity`、`Percent Effect`），"
          "说明它是原样收录的自由文本，**含义随实验而异，不能跨实验汇总**。")
        r("")
        r.table(
            ["standard_type", "数量", "", "官方定义（原文）", "标准单位", "中文说明"],
            [[f"`{x[0]}`", fmt(x[1]), bar(x[1], vmax, 10),
              std_off.get(x[0], "—"),
              "／".join(sorted(u for u in std_units.get(x[0], set()) if u)) or "—",
              ACT_TYPE_DESC.get(x[0], "")] for x in rows],
            aligns=["---", "---:", "---", "---", "---", "---"],
        )
        n_std = sum(1 for x in rows if x[0] in std_off)
        r(f"> 上面 {len(rows)} 个最常见的类型里，只有 **{n_std}** 个有官方标准化规则。"
          "越靠后的类型越可能是未标准化的自由文本。")
        r("")
        # 大小写变体是个真实存在的坑，自动检测并提醒
        allt = db.q("SELECT standard_type, COUNT(*) FROM activities WHERE standard_type IS NOT NULL GROUP BY 1")
        groups: dict[str, list] = {}
        for name, cnt in allt:
            groups.setdefault(str(name).lower(), []).append((name, cnt))
        dup = sorted([g for g in groups.values() if len(g) > 1],
                     key=lambda g: -sum(c for _, c in g))[:8]
        if dup:
            r("> ⚠️ **陷阱：同一指标存在大小写/写法不同的变体**，它们在数据库里是不同的字符串，"
              "`GROUP BY standard_type` 会把它们拆开。本库中最主要的几组：")
            r(">")
            for g in dup:
                r("> - " + "、".join(f"`{n}`（{fmt(c)}）" for n, c in sorted(g, key=lambda x: -x[1])))
            r(">")
            r("> 分析前建议先 `UPPER(standard_type)` 归一，或明确列出你要的写法。")
            r("")

    # 单位
    log("  · 统计单位分布…")
    rows = db.q("SELECT COALESCE(standard_units,'(空)'), COUNT(*) c FROM activities GROUP BY 1 ORDER BY c DESC LIMIT 15")
    if rows:
        r.h(3, "8.3 标准化单位 standard_units（Top 15）")
        r("")
        r("**「浓度一律换算成 nM」是个流传很广但不准确的说法。** 实际规则是："
          "ChEMBL 为每个标准类型规定了一组允许的标准单位，浓度型指标以 **nM** 为主，"
          "但 **`ug.mL-1` 同样是官方认可的标准单位**——"
          "当样品分子量未知时（天然产物提取物、抗菌 MIC 等）只能用质量浓度。")
        r("")
        um = db.one("SELECT COUNT(*) FROM activities WHERE standard_units = 'uM'") or 0
        um_std = db.one("SELECT COUNT(*) FROM activities WHERE standard_units = 'uM' AND standard_flag = 1") or 0
        ug_std = db.one("SELECT COUNT(*) FROM activities WHERE standard_units = 'ug.mL-1' AND standard_flag = 1") or 0
        r(f"本库的实测情况：`ug.mL-1` 中有 **{fmt(ug_std)}** 条 `standard_flag = 1`（确实是标准化过的）；"
          f"而 `uM` 共 {fmt(um)} 条，其中 `standard_flag = 1` 的只有 **{fmt(um_std)}** 条——"
          "**看到 `uM` 基本就意味着这条记录压根没被标准化**，用之前要自己换算并核对。")
        r("")
        r("> 看到 `%` 说明是百分比读数（如抑制率），`(空)` 多半是无量纲比值或定性结论。"
          "**不同单位的数值绝不能混在一起做统计。**")
        r("")
        r.table(["单位", "数量"], [[f"`{x[0]}`", fmt(x[1])] for x in rows], aligns=["---", "---:"])

    # relation
    rows = db.q("SELECT COALESCE(standard_relation,'(空)'), COUNT(*) c FROM activities GROUP BY 1 ORDER BY c DESC LIMIT 10")
    if rows:
        r.h(3, "8.4 关系符 standard_relation")
        r.table(["关系", "数量", "含义"],
                [[f"`{x[0]}`", fmt(x[1]),
                  {"=": "精确值", ">": "大于（未达到效应，实际值更大 → 活性更弱）",
                   "<": "小于", ">=": "大于等于", "<=": "小于等于", "~": "约等于"}.get(x[0], "")]
                 for x in rows],
                aligns=["---", "---:", "---"])

    # pChEMBL 直方图
    log("  · 统计 pChEMBL 分布…")
    rows = db.q("""
        SELECT CAST(pchembl_value AS INT) b, COUNT(*) c
        FROM activities WHERE pchembl_value IS NOT NULL
        GROUP BY b ORDER BY b
    """)
    if rows:
        vmax = max(x[1] for x in rows)
        r.h(3, "8.5 pChEMBL 分布")
        r("")
        r("**读法**：pChEMBL 9 = 1 nM（很强），7 = 100 nM（不错），5 = 10 µM（弱）。"
          "分布通常在 5–8 之间呈钟形，这既反映真实的活性分布，也反映发表偏倚"
          "（太弱的化合物没人报道）。")
        r("")
        r.table(
            ["pChEMBL 区间", "对应浓度", "数量", ""],
            [[f"[{x[0]}, {x[0]+1})", _conc_label(x[0]), fmt(x[1]), bar(x[1], vmax, 30)] for x in rows],
            aligns=["---", "---", "---:", "---"],
        )

    # data_validity_comment 明细
    rows = db.q("SELECT data_validity_comment, COUNT(*) c FROM activities WHERE data_validity_comment IS NOT NULL GROUP BY 1 ORDER BY c DESC")
    if rows:
        dv_off = {x[0]: x[1] for x in db.q("SELECT data_validity_comment, description FROM data_validity_lookup")}
        r.h(3, "8.6 被标记为可疑的数据")
        r("")
        r("> 注意 `Manually validated` 是**正面**标记（已对照原文核实无误），"
          "不要因为它出现在这一列就一并剔除。官方描述读自 `data_validity_lookup` 表。")
        r("")
        r.table(["data_validity_comment", "数量", "官方描述（原文）"],
                [[x[0], fmt(x[1]), dv_off.get(x[0], "—")] for x in rows],
                aligns=["---", "---:", "---"])

    if has_mod:
        rows = db.q("SELECT modality, COUNT(*) c FROM activities WHERE modality IS NOT NULL GROUP BY 1 ORDER BY c DESC")
        if rows:
            r.h(3, "8.7 modality（作用模态）")
            r("")
            r("ChEMBL 37 新增字段，标注化合物的设计模态。目前主要值是"
              "「靶向蛋白降解」（PROTAC、分子胶这类不是抑制靶点、而是诱导其被降解的分子）。"
              "注意它与 `action_type` 不同：标了 modality 的化合物**未必**有活性。")
            r("")
            r.table(["modality", "数量"], [[x[0], fmt(x[1])] for x in rows], aligns=["---", "---:"])

    if args.deep:
        log("  · [--deep] 活性数最多的靶点（大连接）…")
        rows = db.q(f"""
            SELECT td.chembl_id, td.pref_name, td.organism, COUNT(*) c
            FROM activities act
            JOIN assays a ON act.assay_id = a.assay_id
            JOIN target_dictionary td ON a.tid = td.tid
            GROUP BY td.chembl_id, td.pref_name, td.organism
            ORDER BY c DESC LIMIT {args.top}
        """)
        if rows:
            r.h(3, f"8.8 活性数据点最多的靶点（Top {args.top}）")
            r("")
            r("**这张表本身就是一堂课。** 排在最前面的往往不是什么热门蛋白，而是 "
              "`Unchecked`、`ADMET`、细胞系（HepG2、MCF7）和整个物种（疟原虫、金黄色葡萄球菌）——"
              "因为大规模高通量筛选和表型筛选贡献了海量数据点，但它们的靶点归属是模糊的。"
              "**数据量大 ≠ 数据可用**。要找某个具体蛋白的可建模数据，"
              "必须叠加 `target_type = 'SINGLE PROTEIN'` 和 `confidence_score >= 8`。")
            r("")
            r.table(["ChEMBL ID", "靶点名", "物种", "活性数"],
                    [[f"`{x[0]}`", x[1], x[2], fmt(x[3])] for x in rows],
                    aligns=["---", "---", "---", "---:"])

            log("  · [--deep] 高质量子集下的靶点排行…")
            rows2 = db.q(f"""
                SELECT td.chembl_id, td.pref_name, td.organism, COUNT(*) c
                FROM activities act
                JOIN assays a ON act.assay_id = a.assay_id
                JOIN target_dictionary td ON a.tid = td.tid
                WHERE td.target_type = 'SINGLE PROTEIN'
                  AND a.confidence_score >= 8
                  AND act.pchembl_value IS NOT NULL
                  AND act.standard_relation = '='
                  AND act.data_validity_comment IS NULL
                  AND act.potential_duplicate = 0
                GROUP BY td.chembl_id, td.pref_name, td.organism
                ORDER BY c DESC LIMIT {args.top}
            """)
            if rows2:
                r.h(3, f"8.9 加上质量过滤后的排行（Top {args.top}）")
                r("")
                r("过滤条件：`target_type = 'SINGLE PROTEIN'` + `confidence_score ≥ 8` + 有 pChEMBL + "
                  "`standard_relation = '='` + 排除可疑与重复。")
                r("")
                r.table(["ChEMBL ID", "靶点名", "物种", "可用活性数"],
                        [[f"`{x[0]}`", x[1], x[2], fmt(x[3])] for x in rows2],
                        aligns=["---", "---", "---", "---:"])
                r(f"**注意榜单可能还是不太对劲**：占据前列的（{_names(rows2)} ……）"
                  "未必是最重要的药物靶点，往往只是因为 PubChem 上的大规模高通量筛选"
                  "（一次几十万化合物）恰好打了这些靶点。质量标记只能保证「这条数据本身可信」，"
                  "**保证不了「这批数据代表了该靶点的研究现状」**。")
                r("")

                # 再叠加数据来源过滤，展示文献数据的样子
                log("  · [--deep] 仅文献来源的靶点排行…")
                rows3 = db.q(f"""
                    SELECT td.chembl_id, td.pref_name, td.organism, COUNT(*) c
                    FROM activities act
                    JOIN assays a ON act.assay_id = a.assay_id
                    JOIN target_dictionary td ON a.tid = td.tid
                    WHERE td.target_type = 'SINGLE PROTEIN'
                      AND a.confidence_score >= 8
                      AND act.pchembl_value IS NOT NULL
                      AND act.standard_relation = '='
                      AND act.data_validity_comment IS NULL
                      AND act.potential_duplicate = 0
                      AND act.src_id = 1
                    GROUP BY td.chembl_id, td.pref_name, td.organism
                    ORDER BY c DESC LIMIT {args.top}
                """)
                if rows3:
                    r.h(3, f"8.10 再叠加「仅科学文献来源」（`src_id = 1`，Top {args.top}）")
                    r("")
                    r(f"再限定到人工审编的文献数据后，榜首变成了 {_names(rows3)} 等"
                      "被药物化学界长期反复研究的靶点。与上一张表对比即可看出，"
                      "**「哪个靶点数据最多」这个问题的答案，取决于你是否把高通量筛选数据算进来**。")
                    r("")
                    r.table(["ChEMBL ID", "靶点名", "物种", "可用活性数"],
                            [[f"`{x[0]}`", x[1], x[2], fmt(x[3])] for x in rows3],
                            aligns=["---", "---", "---", "---:"])
                    r("> **这三张表（8.8 → 8.9 → 8.10）是本报告最重要的一节。** "
                      "同一个数据库，只是换了过滤条件，「最热门靶点」的答案就完全不同。"
                      "用 ChEMBL 时，**先想清楚你的科学问题决定了哪种过滤**，再写 SQL。")
                    r("")


def _names(rows: list, n: int = 3) -> str:
    """从排行结果里取前 n 个靶点名，拼成一句可嵌入正文的话。

    正文里不写死具体靶点，避免换版本/换库后叙述与下方表格矛盾。
    """
    picked = [str(x[1]) for x in rows[:n] if x[1]]
    return "、".join(picked) if picked else "（见下表）"


def _conc_label(b: int) -> str:
    """把 pChEMBL 整数档翻译成浓度范围。"""
    table = {
        11: "≤ 10 pM", 10: "10–100 pM", 9: "0.1–1 nM", 8: "1–10 nM",
        7: "10–100 nM", 6: "0.1–1 µM", 5: "1–10 µM", 4: "10–100 µM",
        3: "0.1–1 mM", 2: "1–10 mM",
    }
    return table.get(b, "—")


def sec_docs_drugs(r: Report, db: DB, args) -> None:
    log("章节 10/10：文献与药物注释")
    r.h(2, "9. 文献与药物注释")
    r("")

    if db.has_table("docs"):
        n_doc = db.count("docs", not args.quick)
        n_pmid = db.one("SELECT COUNT(*) FROM docs WHERE pubmed_id IS NOT NULL") if db.has_col("docs", "pubmed_id") else None
        n_doi = db.one("SELECT COUNT(*) FROM docs WHERE doi IS NOT NULL") if db.has_col("docs", "doi") else None
        r(f"- 文献/数据集条目：**{fmt(n_doc)}**；有 PubMed ID：{fmt(n_pmid)}；有 DOI：{fmt(n_doi)}")
        r("")
        if db.has_col("docs", "doc_type"):
            rows = db.q("SELECT COALESCE(doc_type,'(空)'), COUNT(*) c FROM docs GROUP BY 1 ORDER BY c DESC")
            r("- 文献类型：" + "；".join(f"{x[0]} {fmt(x[1])}" for x in rows))
            r("")
        if db.has_col("docs", "journal"):
            rows = db.q(f"SELECT journal, COUNT(*) c FROM docs WHERE journal IS NOT NULL GROUP BY 1 ORDER BY c DESC LIMIT {args.top}")
            r.h(3, f"9.1 主要期刊（Top {args.top}）")
            r.table(["期刊", "文献数"], [[x[0], fmt(x[1])] for x in rows], aligns=["---", "---:"])
        if db.has_col("docs", "year"):
            rows = db.q("SELECT year, COUNT(*) c FROM docs WHERE year IS NOT NULL GROUP BY 1 ORDER BY 1")
            if rows:
                vmax = max(x[1] for x in rows)
                recent = [x for x in rows if x[0] and int(x[0]) >= 1990]
                r.h(3, "9.2 文献年份分布（1990 起）")
                r.table(["年份", "文献数", ""], [[x[0], fmt(x[1]), bar(x[1], vmax, 24)] for x in recent],
                        aligns=["---:", "---:", "---"])

    # 药物相关表
    r.h(3, "9.3 药物层面的注释表")
    r("")
    r("这几张表只覆盖**已知药物**（不是全部化合物），但信息密度很高，适合做药物重定位、"
      "靶点-疾病关联一类的分析。")
    r("")
    rows = []
    for t, note in [
        ("drug_mechanism", "药物 × 靶点 × 作用类型（抑制剂/激动剂…），人工审编"),
        ("drug_indication", "药物 × 适应症，疾病映射到 MeSH / EFO"),
        ("drug_warning", "安全性警示（黑框警告、撤市）"),
        ("metabolism", "代谢途径：底物 → 产物 → 催化酶"),
        ("molecule_atc_classification", "药物 × WHO ATC 分类"),
        ("products", "FDA 批准的药品（制剂级）"),
        ("formulations", "药品 → 成分化合物"),
    ]:
        if db.has_table(t):
            rows.append([f"`{t}`", fmt(db.count(t, not args.quick)), note])
    r.table(["表", "行数", "说明"], rows, aligns=["---", "---:", "---"])

    if db.has_table("drug_mechanism") and db.has_col("drug_mechanism", "action_type"):
        rows = db.q("SELECT COALESCE(action_type,'(空)'), COUNT(*) c FROM drug_mechanism GROUP BY 1 ORDER BY c DESC LIMIT 15")
        at_off = {x[0]: (x[1], x[2]) for x in db.q("SELECT action_type, description, parent_type FROM action_type")}
        r.h(3, "9.4 药物作用类型 action_type（Top 15）")
        r("")
        r("> 描述与上位归类读自库内的 `action_type` 表。`parent_type` 把细分类型归成"
          "正向调节 / 负向调节 / 其他三大类，做粗粒度分析时用它比用 `action_type` 更稳。")
        r("")
        r.table(["action_type", "数量", "parent_type", "官方描述（原文）"],
                [[f"`{x[0]}`", fmt(x[1]), at_off.get(x[0], ("", "—"))[1] or "—",
                  at_off.get(x[0], ("—", ""))[0]] for x in rows],
                aligns=["---", "---:", "---", "---"])


def sec_recipes(r: Report, db: DB, args) -> None:
    r.h(2, "10. 上手：几个可直接运行的查询")
    r("")
    r("以下 SQL 可直接在本库上运行（`python3 -c` 或 `sqlite3` 均可）。")
    r("")

    r.h(3, "10.1 已知某个基因，取它的高质量活性数据")
    r("")
    r("这是最常用的一条流水线：基因名 → 靶点 → 高可信实验 → 精确的定量活性。")
    r("")
    r("```sql")
    r("SELECT md.chembl_id                AS compound,")
    r("       cs.canonical_smiles         AS smiles,")
    r("       act.standard_type,")
    r("       act.standard_value, act.standard_units,")
    r("       act.pchembl_value,")
    r("       td.pref_name                AS target,")
    r("       d.year, d.pubmed_id")
    r("FROM activities act")
    r("JOIN assays a               ON act.assay_id = a.assay_id")
    r("JOIN target_dictionary td   ON a.tid        = td.tid")
    r("JOIN target_components tc   ON tc.tid       = td.tid")
    r("JOIN component_sequences seq ON seq.component_id = tc.component_id")
    r("JOIN molecule_dictionary md ON act.molregno = md.molregno")
    r("JOIN compound_structures cs ON cs.molregno  = md.molregno")
    r("LEFT JOIN docs d            ON act.doc_id   = d.doc_id")
    r("WHERE seq.accession = 'P35557'          -- UniProt：人葡萄糖激酶 GCK")
    r("  AND td.target_type = 'SINGLE PROTEIN'")
    r("  AND a.confidence_score >= 8           -- 靶点是单一蛋白（9=直接指认, 8=同源映射）")
    r("  AND act.pchembl_value IS NOT NULL     -- 只要能定量的")
    r("  AND act.standard_relation = '='       -- 排除删失数据")
    r("  AND act.data_validity_comment IS NULL -- 排除可疑数据")
    r("  AND act.potential_duplicate = 0       -- 排除重复引用")
    r("ORDER BY act.pchembl_value DESC;")
    r("```")
    r("")

    r.h(3, "10.2 同一化合物-靶点对有多条记录时怎么办")
    r("")
    r("很常见：不同实验室、不同实验条件重复测过。标准做法是**取中位数**，"
      "并丢弃离散度过大的对（例如 max−min > 1 个 log 单位说明数据打架）。")
    r("")
    r("```sql")
    r("SELECT md.chembl_id, COUNT(*) n,")
    r("       ROUND(AVG(act.pchembl_value), 2)  AS mean_p,")
    r("       ROUND(MAX(act.pchembl_value) - MIN(act.pchembl_value), 2) AS spread")
    r("FROM activities act")
    r("JOIN assays a ON act.assay_id = a.assay_id")
    r("JOIN molecule_dictionary md ON act.molregno = md.molregno")
    r("WHERE a.tid = ? AND act.pchembl_value IS NOT NULL")
    r("GROUP BY md.chembl_id")
    r("HAVING n > 1")
    r("ORDER BY spread DESC;")
    r("```")
    r("")

    r.h(3, "10.3 查一个陌生的 CHEMBL ID 是什么")
    r("")
    r("```sql")
    r("SELECT entity_type, entity_id FROM chembl_id_lookup WHERE chembl_id = 'CHEMBL25';")
    r("```")
    r("")

    r.h(3, "10.4 Python + pandas 的读法")
    r("")
    r("```python")
    r("import sqlite3, pandas as pd")
    r(f"con = sqlite3.connect('file:{db.path}?mode=ro', uri=True)   # 只读，避免误写")
    r("df = pd.read_sql_query(open('query.sql').read(), con)")
    r("```")
    r("")
    r("> 30 GB 的库不要 `SELECT *` 全表读进内存。先用 SQL 过滤好，再交给 pandas。")
    r("")

    r.h(2, "11. 使用这份数据前必须知道的坑")
    r("")
    r.table(
        ["坑", "后果", "对策"],
        [
            ["用了 `value` 而不是 `standard_value`", "单位混杂（µM 与 nM 混在一起），结论全错", "永远用 `standard_*` 列"],
            ["忽略 `standard_relation`", "把 `>10 µM` 当成 10 µM，把「无活性」当成「弱活性」", "定量建模时限定 `= '='`；或按删失数据处理"],
            ["不过滤 `confidence_score`", "把归因到细胞、组织、整个生物体的数据当成单一蛋白上的活性", "建模用 `>= 8`（单一蛋白）；只要直接实测则用 `= 9`"],
            ["把 `confidence_score` 当成数据质量分", "误以为低分 = 实验做得差；实际它描述的是靶点粒度与指认方式", "见 §7.2；低分数据在表型筛选场景下完全可用"],
            ["混用不同 `standard_type`", "IC50 与 Ki 与 %抑制率不可比", "分开处理；至少分开 IC50/EC50 与 Ki/Kd"],
            ["忽略 `potential_duplicate` / `data_validity_comment`", "同一数值被重复计数；纳入已知错误值", "两者都加进过滤条件；但注意 `data_validity_comment IS NULL` 会连 `Manually validated` 的记录一起丢掉"],
            ["把化合物记录数当成化合物数", "`compound_records` 是文献级的，数量远大于唯一化合物", "唯一化合物看 `molecule_dictionary`"],
            ["把盐和母体当成两个分子", "同一药物的不同盐型被算作不同化合物", "用 `molecule_hierarchy` 归并到 parent_molregno"],
            ["假设「没有数据 = 没有活性」", "ChEMBL 存在强烈发表偏倚，阴性结果严重缺失", "做机器学习时需谨慎构造负样本"],
            ["跨版本直接比较数字", "官方会回溯性地修订（合并重复、重新指认靶点、改 max_phase）", "论文里注明具体版本号"],
        ],
    )

    r.h(2, "12. 延伸资料")
    r("")
    r("- 官方 schema 文档：https://chembl.gitbook.io/chembl-interface-documentation/db-schema-description")
    r("- 同目录下的 `schema_documentation.txt`（逐字段说明）与 `chembl_*_schema.pdf`（ER 图）")
    r("- 官方 release notes：`chembl_*_release_notes.txt`")
    r("- Web 界面（适合抽查单个化合物/靶点）：https://www.ebi.ac.uk/chembl/")
    r("- 数据许可 CC BY-SA；发表时需引用 ChEMBL 论文并注明 release 号（见 `REQUIRED.ATTRIBUTION`）")
    r("")
    r("---")
    r("")
    r("*本报告由 `chembl_profile.py` 自动生成。表中的中文解释为脚本内置的领域注解，"
      "统计数字均实时查询自上述数据库文件。*")


# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description="调查 ChEMBL SQLite 数据库，生成面向非药化背景读者的 Markdown 报告。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("db_path", help="ChEMBL SQLite 文件路径，如 chembl_37.db")
    p.add_argument("-o", "--out", default=None,
                   help="输出 Markdown 路径（默认：与数据库同目录的 chembl_profile_report.md）")
    p.add_argument("--quick", action="store_true",
                   help="快速模式：用统计表估算行数、跳过 activities 全表扫描（秒级）")
    p.add_argument("--deep", action="store_true",
                   help="深度模式：额外做 activities×assays×targets 大连接（慢，可能数分钟）")
    p.add_argument("--top", type=int, default=20, help="各排行榜取前 N 项（默认 20）")
    p.add_argument("--stdout", action="store_true", help="同时把报告打印到标准输出")
    args = p.parse_args()

    if not os.path.isfile(args.db_path):
        print(f"错误：找不到文件 {args.db_path}", file=sys.stderr)
        return 1

    t0 = time.time()
    log(f"打开数据库（只读）：{args.db_path}")
    try:
        db = DB(args.db_path)
    except sqlite3.Error as e:
        print(f"错误：无法打开数据库 — {e}", file=sys.stderr)
        return 1

    if not db.has_table("activities") or not db.has_table("molecule_dictionary"):
        print("错误：这不像是一个 ChEMBL 数据库（缺少 activities / molecule_dictionary 表）。", file=sys.stderr)
        return 1

    r = Report()
    sec_intro(r, db, args)
    sec_glossary(r, db)
    sec_schema(r, db, args)
    sec_tables(r, db, args)
    sec_sources(r, db, args)
    sec_compounds(r, db, args)
    sec_targets(r, db, args)
    sec_assays(r, db, args)
    sec_activities(r, db, args)
    sec_docs_drugs(r, db, args)
    sec_recipes(r, db, args)

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.db_path)), "chembl_profile_report.md")
    text = r.text()
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)

    if args.stdout:
        print(text)

    log(f"完成，用时 {time.time() - t0:.1f}s")
    print(f"报告已写入：{out}")
    print(f"  行数 {len(r.lines):,}　大小 {human_bytes(len(text.encode('utf-8')))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
