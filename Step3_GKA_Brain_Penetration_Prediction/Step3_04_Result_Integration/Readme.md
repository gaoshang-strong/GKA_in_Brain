# Step3_04：结果合并

## 任务

把 Step3_02（SwissADME，7 批）与 Step3_03（ADMETlab 3.0，16 批）的原始返回结果，
校验后并进 Step3_01 那张 1,274 行的骨架表，产出本轮的终点：一张候选与对照
在同条件下跑完、各项预测值并排的表。

**本步只取数，不下判断。** 阈值、排序、流程验收都属 Step3_05。

## 输入

```text
../Step3_01_Structure_Standardization_and_RDKit_Properties/Step3_01_RDKit_Processed.csv
../Step3_02_SwissADME_Input/    Step3_02_Submission_Manifest.csv
                                Step3_02_Result_Join_Map.csv
                                Step3_02_Manual_Exclusions.csv
                                swissadme_batch1..7.csv          ← 网页导出的原始结果
../Step3_03_ADMETlab_Input/     Step3_03_Submission_Manifest.csv
                                Step3_03_Result_Join_Map.csv
                                ADMetlab_batch1..16.csv          ← 网页导出的原始结果
```

原始结果文件**不得修改**。脚本只读。

## 产物

```text
Step3_04_Integrated_Brain_Penetration_Results.csv   1,274 × 261   主产物
Step3_04_CNS_MPO.py                  CNS MPO 打分函数（拐点 + 原文算例自检）
Step3_04_Batch_Provenance.csv        每批的行数 / mtime / SHA256 / 校验计数
Step3_04_Verification_Failures.csv   结构核对不通过、被剔除的结果行
Step3_04_Anchor_Drift.csv            20 个锚点在各批之间的分歧量
Step3_04_summary.json
Step3_04_Result_Integration.md       报告
```

## 四条不能松的规则

1. **回填前逐行比对「提交的结构」与「返回的结构」的 InChIKey，不符的整行剔除并记录。**
   SwissADME 出过返回樟脑的事故，名称是对的、只有结构错了——按名称或行序对是发现不了的。
2. **加列不删行。** 1,274 行进、1,274 行出；没结果的写 `*_ok = False` + `*_missing_reason`，
   不填默认值。**人工排除与工具失败必须分得开**（读 `Step3_02_Manual_Exclusions.csv`）。
3. **锚点重复要合并并留痕**：数值取中位数、文本取众数，分歧量写进 `Step3_04_Anchor_Drift.csv`。
4. **CNS MPO 的拐点只能来自原文，且自检写死在脚本里。** 拐点取自 `../cn100008c.pdf`
   的 Table 1 + Figure 4，公式取自 Methods eq 1/2；`Step3_04_CNS_MPO.py` 用原文
   Table 4 的算例（输入未舍入）与 Table 3 的三个候选回算，对不上直接退出、拒绝出分。
   **算分不等于套阈值**——`≥4` 这类判据属 Step3_05。

## 运行

```bash
/home/sgao30/micromamba/bin/micromamba run -n GKA_in_Brain python Step3_04_Merge_Results.py
```

约 33 秒（大部分花在 2,961 次 InChIKey 计算上）。
自检不通过会打印 `⚠` 并写进 `summary.json` 的 `selfcheck_problems`。
