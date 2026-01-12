#!/bin/bash

# 输入文件
ID_FILE="/home/gao/data/hsmilk_cp/vfdb/vfdb_high_copy_genes.tsv"
FASTA="/home/gao/databases/VFDB/VFDB_setA_nt.fas"
OUT="/home/gao/data/hsmilk_cp/vfdb/vfdb_high_copy_genes.fa"

# 提取第一列的基因ID (去掉表头)
tail -n +2 "$ID_FILE" | cut -f1 > gene_ids.txt

# 构建匹配ID的正则（以 | 分隔供 awk 使用）
IDS=$(paste -sd'|' gene_ids.txt)

# 用 awk 逐条扫描 fasta，把匹配到ID的序列输出
awk -v ids="$IDS" '
  BEGIN {split(ids, arr, "|"); for (i in arr) keep[arr[i]]=1}
  /^>/ {
    # 取 > 和 (gb| 之间的部分作为ID
    header=$0
    match(header, /^>([^ ]+)\(gb\|/, m)
    if (m[1] in keep) {printit=1; print header}
    else {printit=0}
    next
  }
  printit {print}
' "$FASTA" > "$OUT"

echo "提取完成，结果保存在 $OUT"
