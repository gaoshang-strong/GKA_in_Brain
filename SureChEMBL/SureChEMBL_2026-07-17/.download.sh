#!/usr/bin/env bash
# SureChEMBL 2.0 全量下载（断点续传）。重跑安全：已完整的文件会被跳过。
BASE=https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/2026-07-17
DEST=/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17
cd "$DEST" || exit 1
# 小文件在前，大文件在后：先把结构性字典拿到手
for f in LICENCE fields.parquet biomedical_types.parquet biomedical_entities.parquet \
         fpsim2_fingerprints.h5 biomedical_locations.parquet compounds.parquet \
         patent_compound_map.parquet patents.parquet; do
  want=$(awk -F'\t' -v n="$f" '$1==n{print $2}' .server_manifest.tsv)
  have=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$have" = "$want" ]; then
    echo "[skip] $f 已完整 ($have)"
    continue
  fi
  echo "[get ] $f  ($have / $want)"
  curl -fsS -C - -o "$f" "$BASE/$f" || { echo "[FAIL] $f"; exit 1; }
done
echo "=== 下载结束，核对大小 ==="
fail=0
while IFS=$'\t' read -r f want _; do
  [ "$f" = "file" ] && continue
  have=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$have" = "$want" ]; then echo "OK   $f  $have"
  else echo "BAD  $f  期望 $want 实得 $have"; fail=1; fi
done < .server_manifest.tsv
exit $fail
