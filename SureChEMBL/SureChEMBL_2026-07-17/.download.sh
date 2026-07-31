#!/usr/bin/env bash
# SureChEMBL 2.0 全量下载（断点续传 + 多轮补传）。
# 重跑安全：已完整的文件跳过，不完整的从断点续。
#
# EBI 的连接会中途掉（实测 "OpenSSL SSL_read: No route to host, errno 113"），
# 所以这里做了三层容错：
#   1. curl 自身重试
#   2. 单个文件失败不中断整轮，继续下一个
#   3. 外层多轮扫描，直到全部齐了或达到 MAX_PASS
#
# ⚠ 本机 curl 是 7.68.0，**不支持 --retry-all-errors**（7.71 才加入）。
#   用了会直接报 "option is unknown" 且每轮空转——只用下面这几个老选项。
BASE=https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/2026-07-17
DEST=/ShangGaoAIProjects/GKA_in_Brain/SureChEMBL/SureChEMBL_2026-07-17
MAX_PASS=${MAX_PASS:-30}
cd "$DEST" || exit 1

# 小文件在前，大文件在后：先把结构性字典拿到手
FILES="LICENCE fields.parquet biomedical_types.parquet biomedical_entities.parquet
       fpsim2_fingerprints.h5 biomedical_locations.parquet compounds.parquet
       patent_compound_map.parquet patents.parquet"

want_of() { awk -F'\t' -v n="$1" '$1==n{print $2}' .server_manifest.tsv; }
have_of() { stat -c%s "$1" 2>/dev/null || echo 0; }

for pass in $(seq 1 "$MAX_PASS"); do
  missing=0
  for f in $FILES; do
    want=$(want_of "$f"); have=$(have_of "$f")
    [ "$have" = "$want" ] && continue
    missing=$((missing + 1))
    echo "[pass $pass] $f  $have / $want"
    # --speed-limit/--speed-time：60 秒内低于 10 KB/s 视为卡死，断开由外层重来
    curl -fsS -C - -o "$f" \
         --retry 8 --retry-delay 10 --retry-connrefused \
         --speed-limit 10240 --speed-time 60 \
         "$BASE/$f" || echo "[warn] $f 本轮未完成，下一轮继续"
    now=$(have_of "$f")
    echo "[pass $pass] $f  → $now / $want"
  done
  [ "$missing" = 0 ] && break
  sleep 5
done

echo "=== 核对字节数 ==="
fail=0
for f in $FILES; do
  want=$(want_of "$f"); have=$(have_of "$f")
  if [ "$have" = "$want" ]; then
    printf 'OK   %-32s %s\n' "$f" "$have"
  else
    printf 'BAD  %-32s 期望 %s 实得 %s\n' "$f" "$want" "$have"; fail=1
  fi
done
[ "$fail" = 0 ] && echo "全部完整" || echo "仍有缺口，再跑一次本脚本即可续传"
exit $fail
