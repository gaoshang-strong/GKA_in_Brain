#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step1_03_LLM_Assay_Classification（第 2 阶段：LLM 处理模糊记录）
===============================================================

只处理规则阶段判不明确的 assay（`review_required = TRUE`），调用 DeepSeek
判断实验测的是激活 / 抑制 / 结合 / GKRP 相互作用 / 细胞或表型效应 / 无法判断。

对应 README 的第 2–4 步，因此严格遵守三条约束：

1. **LLM 不负责数据库筛选，也不覆盖规则的明确判定。** 默认只跑规则标记为
   需复核的子集；规则给出 high 置信度的记录原样保留。
2. **必须输出证据句与置信度**，不允许只给标签。脚本会校验 LLM 引用的证据句
   是否**逐字出现在原始描述里**；对不上就判为无效（`llm_evidence_verified=FALSE`），
   该条强制转人工，不采纳其标签。这是防编造的硬约束。
3. **不允许强行分类。** LLM 返回「无法判断」或低置信度时，保持
   `review_required = TRUE`，交人工。

模型
----
`deepseek-v4-pro`（DeepSeek 的推理模型）。密钥读自环境变量 `DEEPSEEK_API_KEY`，
脚本不会打印或写出密钥。注意该模型会先产出 reasoning tokens，
`max_tokens` 给小了会导致正文为空，因此默认给到 3000。

输入 / 输出
-----------
输入：`Step1_03_GCK_Assay_Classification.csv`（规则阶段的完整输出）
输出：
    Step1_03_LLM_Assay_Classification.csv    仅 LLM 处理过的记录及其证据
    Step1_03_Assay_Classification_final.csv  规则 + LLM 合并后的最终分类表
    Step1_03_LLM_Assay_Classification.md     报告
    .llm_cache.json                          按输入哈希缓存，重跑不重复计费

用法
----
    python3 Step1_03_LLM_Assay_Classification.py --dry-run   # 只打印提示词，不调用
    python3 Step1_03_LLM_Assay_Classification.py
    python3 Step1_03_LLM_Assay_Classification.py --limit 5   # 先小批量试
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

try:
    import requests
except ImportError:
    print("错误：需要 requests（pip install requests）", file=sys.stderr)
    raise SystemExit(1)

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / "Step1_03_GCK_Assay_Classification.csv"
CACHE_PATH = HERE / ".llm_cache.json"

API_URL = "https://api.deepseek.com/chat/completions"
MODEL = "deepseek-v4-pro"

CATEGORIES = ["GCK 激活", "GCK 抑制", "GCK 结合", "GCK–GKRP 相互作用",
              "细胞或表型效应", "无法判断"]

SYSTEM_PROMPT = """你是一名协助整理 ChEMBL 生物活性数据的专家。

任务：根据给定的 assay 信息，判断这个实验**实际测量的内容**属于以下哪一类：

- "GCK 激活"：测量化合物增强葡萄糖激酶(GCK/hexokinase-4)酶活性的能力
- "GCK 抑制"：测量化合物降低 GCK 酶活性的能力
- "GCK 结合"：只测结合/亲和力(Kd、SPR、竞争置换)，不测功能方向
- "GCK–GKRP 相互作用"：测量 GCK 与其调节蛋白 GKRP 之间相互作用的破坏或形成
- "细胞或表型效应"：读数是细胞或整体表型(胰岛素分泌、葡萄糖摄取、核质转位、
  血糖变化)，而非 GCK 酶活本身
- "无法判断"：信息不足或自相矛盾

严格要求：
1. 必须给出 evidence：从"实验描述"中**逐字复制**的一小段原文，作为你判断的依据。
   不得改写、翻译或自行造句。若描述中找不到可支撑判断的原文，evidence 填空字符串
   并把 category 设为"无法判断"。
2. confidence 取 "high" / "medium" / "low"。信息不足时不要硬给 high。
3. 宁可返回"无法判断"，也不要猜测。这些数据将用于科研，错误分类的代价高于漏判。
4. 注意常见陷阱：
   - "expressed in CHO/Sf9 cells" 指表达宿主，不代表这是细胞实验
   - "filter-binding assay"、"binding assay" 可能只是方法名，不代表测的是结合
   - "enzyme affinity for glucose" 指的是底物亲和力，属于酶动力学表征而非化合物结合
   - "in absence of GKRP" 是实验条件，不等于在测 GCK–GKRP 相互作用

只输出 JSON，格式：
{"category": "...", "confidence": "...", "evidence": "...", "reasoning": "一到两句中文说明"}"""


def build_user_prompt(row: dict) -> str:
    parts = [
        f"靶点 ChEMBL ID：{row.get('target_chembl_id', '')}"
        f"（{row.get('target_pref_name', '')}，{row.get('target_type', '')}）",
        f"实验描述：{row.get('assay_description', '')}",
        f"assay_type：{row.get('assay_type', '')}"
        f"（{row.get('assay_type_desc', '')}）",
        f"BAO 实验格式：{row.get('bao_label', '')}",
        f"该 assay 下实测的 standard_type 分布：{row.get('standard_types', '') or '无'}",
    ]
    comments = row.get("activity_comments", "")
    if comments:
        try:
            items = json.loads(comments)
            txt = "；".join(f"{i['comment']}({i['n']}条)" for i in items[:8])
        except (json.JSONDecodeError, KeyError, TypeError):
            txt = comments[:300]
        parts.append(f"activity_comment 备注：{txt}")
    parts.append(f"规则阶段的暂定结论：{row.get('assay_category', '')}"
                 f"（置信度 {row.get('classification_confidence', '')}），"
                 f"理由：{row.get('classification_reason', '')}")
    parts.append("\n请独立判断，不必迁就规则阶段的结论。")
    return "\n".join(parts)


def _norm(s: str) -> str:
    """比对证据句时忽略大小写与空白差异。"""
    return re.sub(r"\s+", " ", (s or "").lower()).strip()


def verify_evidence(evidence: str, row: dict) -> bool:
    """证据必须逐字出自实验描述（或 standard_types / 备注）。"""
    ev = _norm(evidence)
    if not ev:
        return False
    haystack = _norm(" || ".join([
        row.get("assay_description", ""),
        row.get("standard_types", ""),
        row.get("activity_comments", ""),
    ]))
    return ev in haystack


def load_cache() -> dict:
    if CACHE_PATH.is_file():
        try:
            return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}
    return {}


def save_cache(cache: dict) -> None:
    CACHE_PATH.write_text(json.dumps(cache, ensure_ascii=False, indent=1),
                          encoding="utf-8")


def call_llm(prompt: str, api_key: str, max_tokens: int,
             retries: int = 3) -> dict:
    payload = {
        "model": MODEL,
        "messages": [{"role": "system", "content": SYSTEM_PROMPT},
                     {"role": "user", "content": prompt}],
        "response_format": {"type": "json_object"},
        "temperature": 0,
        "max_tokens": max_tokens,
    }
    last = None
    for attempt in range(retries):
        try:
            r = requests.post(
                API_URL,
                headers={"Authorization": f"Bearer {api_key}",
                         "Content-Type": "application/json"},
                json=payload, timeout=180)
            if r.status_code == 200:
                data = r.json()
                content = data["choices"][0]["message"]["content"]
                if not (content or "").strip():
                    last = "模型返回空正文（推理 token 可能耗尽，可调大 --max-tokens）"
                else:
                    return {"raw": content, "usage": data.get("usage", {})}
            elif r.status_code in (429, 500, 502, 503, 529):
                last = f"HTTP {r.status_code}"
            else:
                return {"error": f"HTTP {r.status_code}: {r.text[:200]}"}
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(2 * (attempt + 1))
    return {"error": f"重试 {retries} 次仍失败：{last}"}


def parse_response(raw: str) -> dict:
    """解析模型返回的 JSON；容忍被 ``` 包裹的情况。"""
    txt = raw.strip()
    if txt.startswith("```"):
        txt = re.sub(r"^```[a-z]*\s*|\s*```$", "", txt, flags=re.S)
    try:
        obj = json.loads(txt)
    except json.JSONDecodeError:
        m = re.search(r"\{.*\}", txt, re.S)
        if not m:
            return {"error": "返回内容不是 JSON"}
        try:
            obj = json.loads(m.group(0))
        except json.JSONDecodeError as e:
            return {"error": f"JSON 解析失败：{e}"}
    if not isinstance(obj, dict):
        return {"error": "返回的 JSON 不是对象"}
    return obj


def main() -> int:
    p = argparse.ArgumentParser(
        description="用 DeepSeek 处理规则阶段判不明确的 GCK assay（第 2 阶段）。")
    p.add_argument("--in-csv", type=Path, default=DEFAULT_IN)
    p.add_argument("--outdir", type=Path, default=HERE)
    p.add_argument("--all", action="store_true",
                   help="对全部 assay 调用 LLM（默认只跑 review_required=TRUE 的）")
    p.add_argument("--limit", type=int, default=None, help="只处理前 N 条，便于试跑")
    p.add_argument("--max-tokens", type=int, default=3000,
                   help="deepseek-v4-pro 会先耗 reasoning token，给小了正文会空")
    p.add_argument("--dry-run", action="store_true", help="只打印提示词，不调用 API")
    p.add_argument("--no-cache", action="store_true")
    args = p.parse_args()

    if not args.in_csv.is_file():
        print(f"错误：找不到 {args.in_csv}，请先运行 Step1_03_GCK_Assay_Classification.py",
              file=sys.stderr)
        return 1

    with args.in_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        cols = list(reader.fieldnames or [])
        rows = list(reader)

    targets = rows if args.all else [r for r in rows if r.get("review_required") == "TRUE"]
    if args.limit:
        targets = targets[:args.limit]

    print(f"输入：{args.in_csv}")
    print(f"总记录 {len(rows)}，本次交给 LLM 的 {len(targets)} 条"
          f"（{'全部' if args.all else 'review_required=TRUE'}）")
    print(f"模型：{MODEL}\n")

    if args.dry_run:
        for r in targets[:3]:
            print("=" * 70)
            print(f"# {r['assay_chembl_id']}")
            print(build_user_prompt(r))
        print(f"\n（--dry-run，未调用 API。共 {len(targets)} 条待处理）")
        return 0

    api_key = os.environ.get("DEEPSEEK_API_KEY")
    if not api_key:
        print("错误：环境变量 DEEPSEEK_API_KEY 未设置。", file=sys.stderr)
        return 1

    cache = {} if args.no_cache else load_cache()
    results: dict = {}
    n_cached = n_called = 0
    tok_in = tok_out = 0

    for i, r in enumerate(targets, 1):
        prompt = build_user_prompt(r)
        key = hashlib.sha256((MODEL + "\x00" + SYSTEM_PROMPT + "\x00" + prompt)
                             .encode("utf-8")).hexdigest()[:32]
        if key in cache:
            resp, n_cached = cache[key], n_cached + 1
        else:
            resp = call_llm(prompt, api_key, args.max_tokens)
            n_called += 1
            if "error" not in resp:
                cache[key] = resp
                if not args.no_cache:
                    save_cache(cache)
        u = resp.get("usage") or {}
        tok_in += u.get("prompt_tokens", 0)
        tok_out += u.get("completion_tokens", 0)

        out = {"llm_category": "", "llm_confidence": "", "llm_evidence": "",
               "llm_reasoning": "", "llm_evidence_verified": "", "llm_error": ""}
        if "error" in resp:
            out["llm_error"] = resp["error"]
        else:
            obj = parse_response(resp["raw"])
            if "error" in obj:
                out["llm_error"] = obj["error"]
            else:
                cat = str(obj.get("category", "")).strip()
                out["llm_category"] = cat if cat in CATEGORIES else ""
                if cat and cat not in CATEGORIES:
                    out["llm_error"] = f"返回了未定义的类别：{cat}"
                out["llm_confidence"] = str(obj.get("confidence", "")).strip().lower()
                out["llm_evidence"] = str(obj.get("evidence", "")).strip()
                out["llm_reasoning"] = str(obj.get("reasoning", "")).strip()
                out["llm_evidence_verified"] = (
                    "TRUE" if verify_evidence(out["llm_evidence"], r) else "FALSE")
        results[r["assay_chembl_id"]] = out
        status = out["llm_category"] or out["llm_error"][:40] or "?"
        print(f"  [{i}/{len(targets)}] {r['assay_chembl_id']}: {status} "
              f"({out['llm_confidence']}, 证据校验 {out['llm_evidence_verified']})")

    if not args.no_cache:
        save_cache(cache)

    # ---- 合并：规则 high 的保留；LLM 结果只有在证据可核验且非低置信时才采纳 ----
    llm_cols = ["llm_category", "llm_confidence", "llm_evidence", "llm_reasoning",
                "llm_evidence_verified", "llm_error"]
    final_cols = ["final_category", "final_confidence", "final_method",
                  "final_review_required"]
    for r in rows:
        res = results.get(r["assay_chembl_id"])
        for c in llm_cols:
            r[c] = res[c] if res else ""
        if res and res["llm_category"] and res["llm_evidence_verified"] == "TRUE" \
                and res["llm_confidence"] in ("high", "medium") \
                and res["llm_category"] != "无法判断":
            r["final_category"] = res["llm_category"]
            r["final_confidence"] = res["llm_confidence"]
            r["final_method"] = "llm"
            r["final_review_required"] = "FALSE" if res["llm_confidence"] == "high" else "TRUE"
        elif res:
            # LLM 也没能确定，或证据核验不过 → 保留规则暂定值，转人工
            r["final_category"] = r["assay_category"]
            r["final_confidence"] = "low"
            r["final_method"] = "manual_review_required"
            r["final_review_required"] = "TRUE"
        else:
            r["final_category"] = r["assay_category"]
            r["final_confidence"] = r["classification_confidence"]
            r["final_method"] = r.get("classification_method", "rule")
            r["final_review_required"] = r["review_required"]

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_all = args.outdir / "Step1_03_Assay_Classification_final.csv"
    with out_all.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols + llm_cols + final_cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols + llm_cols + final_cols})

    llm_only = [r for r in rows if r["assay_chembl_id"] in results]
    sub_cols = (["assay_chembl_id", "assay_description", "standard_types",
                 "assay_category", "classification_confidence"]
                + llm_cols + final_cols)
    out_llm = args.outdir / "Step1_03_LLM_Assay_Classification.csv"
    with out_llm.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=sub_cols)
        w.writeheader()
        for r in llm_only:
            w.writerow({k: r.get(k, "") for k in sub_cols})

    write_report(rows, llm_only, args.outdir / "Step1_03_LLM_Assay_Classification.md",
                 n_called, n_cached, tok_in, tok_out)

    print(f"\nAPI 调用 {n_called} 次，命中缓存 {n_cached} 次，"
          f"token 输入 {tok_in:,} / 输出 {tok_out:,}")
    print(f"仍需人工审核：{sum(1 for r in rows if r['final_review_required'] == 'TRUE')} / {len(rows)}")
    print(f"\nLLM 结果：{out_llm}")
    print(f"最终分类：{out_all}")
    return 0


def write_report(rows: list, llm_only: list, path: Path,
                 n_called: int, n_cached: int, tok_in: int, tok_out: int) -> None:
    L: list = []
    L.append("# Step1_03 LLM 分类结果（第 2 阶段）")
    L.append("")
    L.append(f"- 模型：`{MODEL}`")
    L.append(f"- 运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- 交给 LLM 的记录：**{len(llm_only)}** / {len(rows)}")
    L.append(f"- API 调用 {n_called} 次（缓存命中 {n_cached} 次），"
             f"token 输入 {tok_in:,} / 输出 {tok_out:,}")
    L.append("")
    L.append("> LLM 只处理规则阶段判不明确的记录，不覆盖规则的 high 置信度结论。"
             "每条结论都要求给出**逐字引自原描述**的证据句；脚本会校验该证据是否真实存在，"
             "核验不通过的一律不采纳标签并转人工。")
    L.append("")

    ver = Counter(r["llm_evidence_verified"] for r in llm_only)
    L.append("## 证据校验")
    L.append("")
    L.append(f"- 通过：**{ver.get('TRUE', 0)}**")
    L.append(f"- 未通过（证据句无法在原文中找到，标签不予采纳）：**{ver.get('FALSE', 0)}**")
    L.append("")

    L.append("## 最终分类分布")
    L.append("")
    L.append("| 类别 | assay 数 | 其中仍需人工审核 |")
    L.append("| --- | ---: | ---: |")
    for c, n in Counter(r["final_category"] for r in rows).most_common():
        need = sum(1 for r in rows
                   if r["final_category"] == c and r["final_review_required"] == "TRUE")
        L.append(f"| {c} | {n} | {need} |")
    L.append("")

    L.append("## 判定方式")
    L.append("")
    L.append("| final_method | 数量 | 含义 |")
    L.append("| --- | ---: | --- |")
    meaning = {"rule": "规则明确判定，未经 LLM",
               "llm": "规则不明确，采纳了 LLM 的判定（证据已核验）",
               "manual_review_required": "规则与 LLM 都未能确定，需人工审核"}
    for m, n in Counter(r["final_method"] for r in rows).most_common():
        L.append(f"| {m} | {n} | {meaning.get(m, '')} |")
    L.append("")

    L.append("## LLM 逐条结果")
    L.append("")
    L.append("| assay | 规则暂定 | LLM 判定 | 置信度 | 证据核验 | 证据句 | LLM 理由 |")
    L.append("| --- | --- | --- | --- | --- | --- | --- |")
    for r in llm_only:
        ev = (r["llm_evidence"] or "")[:90].replace("|", "\\|")
        rs = (r["llm_reasoning"] or "")[:90].replace("|", "\\|")
        err = f" ⚠{r['llm_error'][:40]}" if r["llm_error"] else ""
        L.append(f"| `{r['assay_chembl_id']}` | {r['assay_category']} | "
                 f"{r['llm_category'] or '—'}{err} | {r['llm_confidence'] or '—'} | "
                 f"{r['llm_evidence_verified']} | {ev} | {rs} |")
    L.append("")

    still = [r for r in rows if r["final_review_required"] == "TRUE"]
    if still:
        L.append("## 仍需人工审核")
        L.append("")
        L.append(f"共 **{len(still)}** 条。按 README 第 4 步，这些不允许由 LLM 强行分类。")
        L.append("")
        L.append("| assay | 暂定类别 | 描述 |")
        L.append("| --- | --- | --- |")
        for r in still:
            d = (r["assay_description"] or "")[:110].replace("|", "\\|")
            L.append(f"| `{r['assay_chembl_id']}` | {r['final_category']} | {d} |")
        L.append("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
