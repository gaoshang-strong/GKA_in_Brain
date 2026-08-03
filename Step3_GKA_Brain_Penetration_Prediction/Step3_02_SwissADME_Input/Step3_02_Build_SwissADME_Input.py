#!/usr/bin/env python
"""
Step3_02：生成 SwissADME 提交文件

读 Step3_01_RDKit_Processed.csv（1,274 个分子，已去盐并中性化），
切成若干批提交文件 + 一份回填清单。

SwissADME 的硬约束（取自其 FAQ / Help 页，2026-08-02 核对）：
  · 输入是文本框，一行一个条目，格式 "SMILES<空格>名称"，名称可省略
  · "Do NOT exceed 200 entries per list."
  · "Do NOT launch several calculations at the same time, please wait for the
     whole calculation to be terminated before running the next batch."  → 必须串行
  · 约 1–5 秒/分子；结果面板下方红色 CSV 图标导出

三个设计要点：

1. **名称一律用 mol_id，不用化合物名。**
   名称与 SMILES 以空格分隔，而 1,274 个里有 93 个化合物名本身含空格
   （如 "MK-0941 FREE BASE"），直接写进去会被截断，结果对不回来。

2. **每批都放同一组锚点分子（重复提交）。**
   SwissADME 是网页模型、可能静默更新，7 批不是同时跑完的。
   同一个锚点在不同批次给出不同结果，就是批次漂移的直接证据；
   没有锚点则漂移不可见，而漂移会被误读成候选与对照的真实差异。

3. **各批的集合构成必须均匀混合。**
   若把 GKA 候选和对照分别装进不同批次，一旦发生漂移，
   就无法把它与「候选 vs 对照」的真实差异区分开——两者完全混杂。
   因此按 set 分层后轮转分配。

用法：
  micromamba run -n GKA_in_Brain python Step3_02_Build_SwissADME_Input.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_IN = (HERE.parent / "Step3_01_Structure_Standardization_and_RDKit_Properties"
              / "Step3_01_RDKit_Processed.csv")

BATCH_SIZE = 200          # SwissADME FAQ 明文上限，不要调大
N_ANCHORS = 20            # 每批重复提交的锚点数
N_ANCHOR_FRIDEN = 12      # 锚点里取自 Fridén 定量集的个数（覆盖 Kp,uu 全跨度）


def pick_anchors(df: pd.DataFrame) -> pd.Index:
    """选批次锚点：覆盖预测值可能的全跨度，且结构上不扎堆。

    取自对照集而非 GKA 候选——锚点的作用是监测工具本身的稳定性，
    用带实测值的分子才能同时看出「有没有漂」和「漂向哪边」。
    """
    fr = df[(df.set == "bbb_control_friden") & df.kpuu_brain.notna()]
    fr = fr.sort_values("kpuu_brain")
    # 沿 Kp,uu 等距取点，两端必取（Kp,uu 0.006 与 2.0 是跨度的两极）
    idx = np.linspace(0, len(fr) - 1, min(N_ANCHOR_FRIDEN, len(fr))).round().astype(int)
    anchors = list(fr.index[np.unique(idx)])

    # 其余名额给 B3DB 两端，让锚点也覆盖 logBB 判据的极值
    need = N_ANCHORS - len(anchors)
    if need > 0:
        b = df[(df.set == "bbb_control_b3db") & df.logbb.notna()].sort_values("logbb")
        half = need // 2
        anchors += list(b.index[:half]) + list(b.index[-(need - half):])
    return pd.Index(anchors[:N_ANCHORS])


def assign_batches(df: pd.DataFrame, anchor_idx: pd.Index,
                   n_batches: int) -> pd.Series:
    """把非锚点分子按 set 分层后轮转分配到各批，保证各批构成相近。"""
    rest = df.drop(index=anchor_idx)
    # 分层：同一 set 内按 mol_id 稳定排序后依次发牌，避免任何一批被某个 set 占满
    order = rest.sort_values(["set", "mol_id"]).index
    batch = pd.Series(index=order, dtype="int64")
    for i, ix in enumerate(order):
        batch[ix] = i % n_batches + 1
    return batch


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_IN)
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    args = ap.parse_args()

    started = datetime.now()
    if not args.input.exists():
        sys.exit(f"[FATAL] 找不到输入：{args.input}（先跑 Step3_01b）")

    print("[1/5] 读入并挑出可提交的分子 …")
    df = pd.read_csv(args.input)
    ok = df.std_ok & df.std_smiles.notna()
    print(f"      {len(df)} 行；标准化成功 {int(ok.sum())}")
    if (~ok).any():
        print(f"      ⚠ {int((~ok).sum())} 行无法提交，已在清单中标注原因")

    # 同一结构只提交一次，结果按 std_smiles 回填给所有共享它的行。
    # 这不是删行——1,274 行原样保留在清单里，只是不重复占用提交名额。
    sub = df[ok].drop_duplicates("std_smiles").copy()
    print(f"      唯一结构 {len(sub)}（{int(ok.sum()) - len(sub)} 个重复结构不重复提交）")

    print("[2/5] 选锚点 …")
    anchor_idx = pick_anchors(sub)
    anchors = sub.loc[anchor_idx]
    print(f"      {len(anchors)} 个：Fridén {int((anchors.set=='bbb_control_friden').sum())}"
          f" / B3DB {int((anchors.set=='bbb_control_b3db').sum())}")

    print("[3/5] 分批 …")
    capacity = args.batch_size - len(anchors)
    n_batches = math.ceil((len(sub) - len(anchors)) / capacity)
    batch = assign_batches(sub, anchor_idx, n_batches)
    sub["batch"] = batch
    sub.loc[anchor_idx, "batch"] = 0          # 0 = 锚点，出现在每一批
    sub["is_batch_anchor"] = sub.index.isin(anchor_idx)
    print(f"      {n_batches} 批 × 每批 {capacity} 新分子 + {len(anchors)} 锚点 "
          f"= 每批 ≤ {capacity + len(anchors)}")

    print("[4/5] 写提交文件 …")
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for b in range(1, n_batches + 1):
        part = pd.concat([sub[sub.batch == b], anchors])
        lines = [f"{r.std_smiles} {r.mol_id}" for r in part.itertuples()]
        assert len(lines) <= args.batch_size, f"批 {b} 超过 {args.batch_size} 条"
        p = args.outdir / f"Step3_02_SwissADME_batch{b:02d}.smi"
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"      {p.name}  {len(lines)} 条"
              f"（新 {int((sub.batch == b).sum())} + 锚点 {len(anchors)}）")
        for pos, r in enumerate(part.itertuples(), start=1):
            manifest_rows.append({
                "batch": b, "line_no": pos, "submitted_name": r.mol_id,
                "std_smiles": r.std_smiles, "set": r.set,
                "is_batch_anchor": bool(r.is_batch_anchor),
            })
    manifest = pd.DataFrame(manifest_rows)

    # 回填清单：1,274 行全在，标出各自结构由哪一批提交
    smi2batch = (sub[~sub.is_batch_anchor].set_index("std_smiles").batch.to_dict())
    back = df[["mol_id", "set", "source_id", "compound_name", "std_smiles",
               "std_ok", "inchikey"]].copy()
    back["submitted_as"] = back.std_smiles.map(
        sub.set_index("std_smiles").mol_id.to_dict())
    back["is_batch_anchor"] = back.std_smiles.map(
        sub.set_index("std_smiles").is_batch_anchor.to_dict()).fillna(False)
    back["batch"] = back.std_smiles.map(smi2batch)
    back.loc[back.is_batch_anchor, "batch"] = 0
    back["not_submitted_reason"] = np.where(
        back.std_ok, "", "标准化失败，无可提交结构")
    back.loc[back.std_ok & (back.mol_id != back.submitted_as),
             "not_submitted_reason"] = "与 submitted_as 结构相同，结果按结构回填"

    p_man = args.outdir / "Step3_02_Submission_Manifest.csv"
    p_back = args.outdir / "Step3_02_Result_Join_Map.csv"
    manifest.to_csv(p_man, index=False)
    back.to_csv(p_back, index=False)
    print(f"      {p_man.name}  ({manifest.shape[0]} 条提交记录)")
    print(f"      {p_back.name}  ({back.shape[0]} 行，1,274 行全保留)")

    print("[5/5] 自检 …")
    problems = []
    over = manifest.groupby("batch").size()
    if (over > args.batch_size).any():
        problems.append(f"有批次超过 {args.batch_size} 条：{over[over > args.batch_size].to_dict()}")
    if manifest.submitted_name.str.contains(" ").any():
        problems.append("提交名称含空格，SwissADME 会截断")
    missing = set(sub.mol_id) - set(manifest.submitted_name)
    if missing:
        problems.append(f"{len(missing)} 个唯一结构没有出现在任何批次中")
    n_anchor_per_batch = manifest[manifest.is_batch_anchor].groupby("batch").size()
    if n_anchor_per_batch.nunique() != 1:
        problems.append(f"各批锚点数不一致：{n_anchor_per_batch.to_dict()}")
    # 各批的集合构成应当相近，否则漂移与真实差异分不开
    mix = manifest[~manifest.is_batch_anchor].pivot_table(
        index="batch", columns="set", aggfunc="size", fill_value=0)
    if len(mix) > 1 and (mix.max() - mix.min() > 2).any():
        problems.append(f"各批 set 构成差异过大：\n{mix}")
    if problems:
        for x in problems:
            print(f"  ⚠ {x}")
    else:
        print("  ✓ 全部通过")
        print(f"  ✓ 各批构成均匀：\n{mix.to_string()}")

    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "input": str(args.input),
        "tool": "SwissADME (https://www.swissadme.ch/)",
        "constraints_checked_on": "2026-08-02",
        "batch_size_limit": args.batch_size,
        "n_rows_in": int(len(df)),
        "n_unique_structures": int(len(sub)),
        "n_batches": int(n_batches),
        "n_anchors_per_batch": int(len(anchors)),
        "total_submissions": int(len(manifest)),
        "anchor_mol_ids": list(anchors.mol_id),
        "selfcheck_problems": problems,
    }
    (args.outdir / "Step3_02_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n{n_batches} 个提交文件，共 {len(manifest)} 条"
          f"（{len(sub)} 唯一结构 + 锚点重复 {len(manifest) - len(sub)}）")


if __name__ == "__main__":
    main()
