#!/usr/bin/env python
"""
Step3_03：生成 ADMETlab 3.0 提交文件

读 Step3_01_RDKit_Processed.csv（1,274 个分子，已去盐并中性化），
切批写出提交文件 + 行序清单。

ADMETlab 3.0 的约束（2026-08-02 核对其 screening 页与 example 文件）：
  · 批量上限页面明写 "Submit SMILES(MAX 1000)"
  · 接受 SDF / TXT / CSV
      example.csv = 单列，表头 "SMILES"，一行一个
      example.txt = 裸 SMILES，无表头
  · 有 API（/api/openapi.json），但实测 /api/single/admet 返回服务端 500，
    /api/washmol 可用。本步仍走网页上传，API 留待其修复。

⚠⚠ 与 SwissADME 最大的不同：**输入里没有名称/ID 字段。**
   SwissADME 可以写 "SMILES mol_id"，这里只能给 SMILES。
   因此结果只能靠**行序**对回去——一旦工具静默丢行或重排，
   该行之后的全部结果都会整体错位，而且不会报错。

   三重防护，缺一不可：
   1. 提交清单逐行记录 (batch, line_no, mol_id, std_smiles)；
   2. 同一批内 SMILES 保证唯一（已按结构去重），使「按结构回填」成为
      行序之外的独立校验路径；
   3. 结果回来后必须先核对行数与结构，对不上就整批作废重跑，
      绝不允许按行序硬对。

⚠ 人工排除项从 Step3_02_Manual_Exclusions.csv 读入。
   B3D_0012（MW 1802.7）提交 SwissADME 时已被人工删除；
   此处默认同样排除——它是全表唯一 MW > 1000 的分子，
   若 ADMETlab 也静默丢弃它，会直接破坏行序对齐。
   要包含它用 --include-excluded。

用法：
  micromamba run -n GKA_in_Brain python Step3_03_Build_ADMETlab_Input.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_IN = (HERE.parent / "Step3_01_Structure_Standardization_and_RDKit_Properties"
              / "Step3_01_RDKit_Processed.csv")
EXCLUSIONS = HERE.parent / "Step3_02_SwissADME_Input" / "Step3_02_Manual_Exclusions.csv"
SWISS_SUMMARY = HERE.parent / "Step3_02_SwissADME_Input" / "Step3_02_summary.json"

BATCH_LIMIT = 1000        # 页面明写上限
# ⚠ 但页面上限不等于实际能跑通的量：实测 646 个/批时网页返回 504 Gateway Time-out，
#   站点本身 HTTP 200，是计算后端超时。所以真正的约束是**服务端算得完**，不是那个 1000。
#   默认切到 100，宁可多提交几次也别整批超时重来。
BATCH_SIZE = 100
N_ANCHORS_TARGET = 20     # 与 Step3_02 同一组锚点；批次多了漂移风险更大，值得保留


def load_anchors(sub: pd.DataFrame) -> pd.Index:
    """沿用 Step3_02 选定的同一组锚点，两个工具用同一把尺子才好互相印证。"""
    if SWISS_SUMMARY.exists():
        ids = json.loads(SWISS_SUMMARY.read_text(encoding="utf-8"))["anchor_mol_ids"]
        hit = sub[sub.mol_id.isin(ids)]
        if len(hit):
            return hit.index
    # 拿不到就退回：Fridén 全取（42 个里覆盖 Kp,uu 全跨度）
    return sub[sub.set == "bbb_control_friden"].index


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_IN)
    ap.add_argument("--exclusions", type=Path, default=EXCLUSIONS)
    ap.add_argument("--include-excluded", action="store_true",
                    help="连人工排除项一并提交")
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE,
                    help="每批分子数（含锚点）。实测 646 会 504，默认 100")
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    started = datetime.now()
    if not args.input.exists():
        sys.exit(f"[FATAL] 找不到输入：{args.input}（先跑 Step3_01b）")

    print("[1/5] 读入 …")
    df = pd.read_csv(args.input)
    ok = df.std_ok & df.std_smiles.notna()
    print(f"      {len(df)} 行；标准化成功 {int(ok.sum())}")

    excluded: list[str] = []
    if args.exclusions.exists() and not args.include_excluded:
        ex = pd.read_csv(args.exclusions)
        excluded = list(ex.mol_id)
        print(f"      人工排除 {len(excluded)} 个：{', '.join(excluded)}")
        for _, r in ex.iterrows():
            print(f"        {r.mol_id} {r.compound_name} — {r.reason}（MW {r.mw}）")

    # 同一结构只提交一次；1,274 行全部保留在回填表里
    sub = df[ok & ~df.mol_id.isin(excluded)].drop_duplicates("std_smiles").copy()
    print(f"      待提交唯一结构 {len(sub)}")

    print("[2/5] 取锚点（与 Step3_02 同一组）…")
    anchor_idx = load_anchors(sub)
    anchors = sub.loc[anchor_idx]
    print(f"      {len(anchors)} 个")

    print("[3/5] 分批（等分，按 set 分层轮转）…")
    rest = sub.drop(index=anchor_idx).sort_values(["set", "mol_id"])
    capacity = args.batch_size - len(anchors)
    if capacity <= 0:
        sys.exit(f"[FATAL] 批大小 {args.batch_size} 容不下 {len(anchors)} 个锚点")
    n_batches = math.ceil(len(rest) / capacity)
    batch = {ix: i % n_batches + 1 for i, ix in enumerate(rest.index)}
    sub["batch"] = sub.index.map(lambda i: 0 if i in set(anchor_idx) else batch[i])
    sub["is_batch_anchor"] = sub.index.isin(anchor_idx)
    print(f"      {n_batches} 批 × 每批 ≤{capacity} 个新分子 + {len(anchors)} 锚点 "
          f"= 每批 ≤{args.batch_size}")

    print("[4/5] 写提交文件 …")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for b in range(1, n_batches + 1):
        part = pd.concat([sub[sub.batch == b], anchors])
        if len(part) > args.batch_size:
            sys.exit(f"[FATAL] 批 {b} 有 {len(part)} 条，超过设定的 {args.batch_size}")
        # CSV：单列 + 表头 SMILES，与官方 example.csv 一致
        pc = args.outdir / f"Step3_03_ADMETlab_batch{b:02d}.csv"
        pc.write_text("SMILES\n" + "\n".join(part.std_smiles) + "\n", encoding="utf-8")
        # TXT：裸 SMILES 无表头，与官方 example.txt 一致，作为备选格式
        pt = args.outdir / f"Step3_03_ADMETlab_batch{b:02d}.txt"
        pt.write_text("\n".join(part.std_smiles) + "\n", encoding="utf-8")
        print(f"      {pc.name} / {pt.name}  {len(part)} 条"
              f"（新 {int((sub.batch == b).sum())} + 锚点 {len(anchors)}）")
        for pos, r in enumerate(part.itertuples(), start=1):
            rows.append({"batch": b, "line_no": pos, "mol_id": r.mol_id,
                         "std_smiles": r.std_smiles, "inchikey": r.inchikey,
                         "set": r.set, "is_batch_anchor": bool(r.is_batch_anchor)})
    man = pd.DataFrame(rows)

    back = df[["mol_id", "set", "source_id", "compound_name",
               "std_smiles", "std_ok", "inchikey"]].copy()
    b2 = sub.set_index("std_smiles")
    back["submitted_as"] = back.std_smiles.map(b2.mol_id.to_dict())
    back["batch"] = back.std_smiles.map(b2.batch.to_dict())
    back["is_batch_anchor"] = back.std_smiles.map(
        b2.is_batch_anchor.to_dict()).fillna(False)
    back["not_submitted_reason"] = ""
    back.loc[back.mol_id.isin(excluded), "not_submitted_reason"] = "人工排除（见 Step3_02_Manual_Exclusions.csv）"
    back.loc[~back.std_ok, "not_submitted_reason"] = "标准化失败"
    back.loc[back.std_ok & back.submitted_as.notna()
             & (back.mol_id != back.submitted_as), "not_submitted_reason"] = (
        "与 submitted_as 结构相同，结果按结构回填")

    p_man = args.outdir / "Step3_03_Submission_Manifest.csv"
    p_back = args.outdir / "Step3_03_Result_Join_Map.csv"
    man.to_csv(p_man, index=False)
    back.to_csv(p_back, index=False)
    print(f"      {p_man.name}  ({len(man)} 条，含 line_no 行序)")
    print(f"      {p_back.name}  ({len(back)} 行，1,274 行全保留)")

    print("[5/5] 自检 …")
    problems = []
    sz = man.groupby("batch").size()
    if (sz > args.batch_size).any():
        problems.append(f"有批次超过 {args.batch_size}：{sz.to_dict()}")
    # 行序回填的前提：同一批内 SMILES 必须唯一，否则按结构校验时无法定位
    for b, g in man.groupby("batch"):
        if g.std_smiles.duplicated().any():
            problems.append(f"批 {b} 内有重复 SMILES，按结构回填会有歧义")
    if man.line_no.min() != 1:
        problems.append("line_no 未从 1 开始")
    for b, g in man.groupby("batch"):
        if list(g.line_no) != list(range(1, len(g) + 1)):
            problems.append(f"批 {b} 的 line_no 不连续")
    missing = set(sub.mol_id) - set(man.mol_id)
    if missing:
        problems.append(f"{len(missing)} 个结构未出现在任何批次")
    na = man[man.is_batch_anchor].groupby("batch").size()
    if na.nunique() != 1:
        problems.append(f"各批锚点数不一致：{na.to_dict()}")
    mix = man[~man.is_batch_anchor].pivot_table(
        index="batch", columns="set", aggfunc="size", fill_value=0)
    if len(mix) > 1 and (mix.max() - mix.min() > 2).any():
        problems.append(f"各批构成差异过大：\n{mix}")
    if problems:
        for x in problems:
            print(f"  ⚠ {x}")
    else:
        print("  ✓ 全部通过")
        print(f"  ✓ 各批构成：\n{mix.to_string()}")

    summary = {
        "generated_at": started.strftime("%Y-%m-%d %H:%M:%S"),
        "input": str(args.input),
        "tool": "ADMETlab 3.0 (https://admetlab3.scbdd.com/server/screening)",
        "constraints_checked_on": "2026-08-02",
        "batch_limit_stated_by_site": BATCH_LIMIT,
        "batch_size_used": int(args.batch_size),
        "why_smaller": "实测 646/批返回 504 Gateway Time-out（站点 HTTP 200，计算后端超时）",
        "input_format": "CSV 单列表头 SMILES / TXT 裸 SMILES；**不支持名称字段**",
        "join_strategy": "行序 (line_no) 为主，SMILES/InChIKey 为独立校验",
        "n_rows_in": int(len(df)),
        "manually_excluded": excluded,
        "n_unique_structures_submitted": int(len(sub)),
        "n_batches": int(n_batches),
        "n_anchors_per_batch": int(len(anchors)),
        "total_submissions": int(len(man)),
        "anchor_mol_ids": list(anchors.mol_id),
        "selfcheck_problems": problems,
    }
    (args.outdir / "Step3_03_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n{n_batches} 批，共 {len(man)} 条"
          f"（{len(sub)} 唯一结构 + 锚点重复 {len(man) - len(sub)}）")


if __name__ == "__main__":
    main()
