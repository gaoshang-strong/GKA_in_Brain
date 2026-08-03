#!/usr/bin/env python
"""
Step3_01a：构建统一输入表

把两批分子合成一张表，作为 Step3 后续所有计算（RDKit、SwissADME、ADMETlab）的
**唯一输入**：

  1. GKA 候选   787 个   ← Step1_GKA_Candidates_with_Properties.csv
  2. 入脑对照   487 个   ← Step3_00_BBB_Control_Set.csv
                          （445 B3DB 分类 + 42 Fridén 定量）

合成一张表而不是分开跑，是为了保证**候选与对照在完全相同的条件下走完全程**——
分批提交给网页工具时，对照必须与候选同批，否则对照标定的结论套不到候选上。

⚠ 只搬运，不改结构：本步不做任何标准化，`input_smiles` 逐字来自各自上游。
   去盐、中性化在 Step3_01b 做，两步分开以便核对改了什么。

用法：
  micromamba run -n GKA_in_Brain python Step3_01a_Build_Input_Table.py
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
GKA_CSV = (HERE.parents[1] / "Step1_Find_GKA_from_ChEMBL"
           / "Step1_GKA_Candidates_with_Properties.csv")
CTRL_CSV = (HERE.parent / "Step3_00_BBB_Control_Set"
            / "Step3_00_BBB_Control_Set.csv")

# 从上游带过来的注解列。只带下游要用的，其余保持在上游表里按 ID 关联，
# 避免这张表越长越宽、失去「输入清单」的性质。
GKA_CARRY = ["molecule_pref_name", "max_phase", "source", "curated_direction",
             "priority", "potency_band", "pactivity_median",
             "is_positive_control", "control_name", "parent_chembl_id", "is_salt"]


def build_gka(path: Path) -> pd.DataFrame:
    g = pd.read_csv(path)
    out = pd.DataFrame({
        "set": "gka_candidate",
        "source_id": g.molecule_chembl_id,
        "compound_name": g.molecule_pref_name,
        "input_smiles": g.canonical_smiles,
    })
    for c in GKA_CARRY:
        out[f"gka_{c}"] = g[c] if c in g.columns else pd.NA
    # 对照专用列留空，保证两批合表后列对齐
    for c in ["control_class", "measure_basis", "logbb", "kpuu_brain",
              "b3db_label", "b3db_group"]:
        out[c] = pd.NA
    return out


def build_controls(path: Path) -> pd.DataFrame:
    c = pd.read_csv(path)
    out = pd.DataFrame({
        "set": c.control_set.map({
            "b3db_classification": "bbb_control_b3db",
            "friden_quantitative": "bbb_control_friden",
        }).fillna(c.control_set),
        # B3DB 侧没有 ChEMBL ID，用 control_id；Fridén 侧优先用 ChEMBL ID
        "source_id": c.molecule_chembl_id.fillna(c.control_id),
        "compound_name": c.compound_name,
        "input_smiles": c.smiles,
    })
    for col in GKA_CARRY:
        out[f"gka_{col}"] = pd.NA
    for col in ["control_class", "measure_basis", "logbb", "kpuu_brain",
                "b3db_label", "b3db_group"]:
        out[col] = c[col] if col in c.columns else pd.NA
    out["control_id"] = c.control_id
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gka", type=Path, default=GKA_CSV)
    ap.add_argument("--controls", type=Path, default=CTRL_CSV)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    started = datetime.now()
    for f in (args.gka, args.controls):
        if not f.exists():
            sys.exit(f"[FATAL] 找不到输入：{f}")

    print("[1/3] 读入两批分子 …")
    gka = build_gka(args.gka)
    ctrl = build_controls(args.controls)
    print(f"      GKA 候选 {len(gka)}；入脑对照 {len(ctrl)}")

    print("[2/3] 合表并编号 …")
    df = pd.concat([gka, ctrl], ignore_index=True)
    if "control_id" not in df.columns:
        df["control_id"] = pd.NA
    prefix = {"gka_candidate": "GKA", "bbb_control_b3db": "B3D",
              "bbb_control_friden": "FRI"}
    ids, counters = [], {}
    for s in df["set"]:
        counters[s] = counters.get(s, 0) + 1
        ids.append(f"{prefix.get(s, 'MOL')}_{counters[s]:04d}")
    df.insert(0, "mol_id", ids)

    # SMILES 缺失是事实，如实标出，不丢行
    df["input_smiles_missing"] = df.input_smiles.isna()
    n_missing = int(df.input_smiles_missing.sum())
    if n_missing:
        print(f"      ⚠ {n_missing} 行没有 SMILES，已标 input_smiles_missing")

    front = ["mol_id", "set", "source_id", "compound_name", "input_smiles"]
    df = df[front + [c for c in df.columns if c not in front]]

    print("[3/3] 写出 …")
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / "Step3_01_Molecule_Input.csv"
    df.to_csv(out, index=False)
    print(f"      {out.name}  ({df.shape[0]} × {df.shape[1]})")

    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {"gka": str(args.gka), "controls": str(args.controls)},
        "counts": df["set"].value_counts().to_dict(),
        "total": int(len(df)),
        "smiles_missing": n_missing,
    }
    (args.outdir / "Step3_01a_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\n组成：")
    for k, v in df["set"].value_counts().items():
        print(f"  {k:22s} {v}")
    print(f"  {'合计':22s} {len(df)}")


if __name__ == "__main__":
    main()
