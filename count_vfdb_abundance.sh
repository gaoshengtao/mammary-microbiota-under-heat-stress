#!/bin/bash
set -euo pipefail

# 输入文件
GENOME_FA="/home/gao/data/hsmilk_cp/vfdb/vfdb_high_copy_genes.fa"
HOSTFREE_DIR="/home/gao/data/hsmilk_cp/hostfree"
OUT_DIR="/home/gao/data/hsmilk_cp/vfdb/mapping_results"
mkdir -p "$OUT_DIR"

# 1. 建立索引 (只需要运行一次)
if [ ! -f "${GENOME_FA}.bwt" ]; then
    echo "Building bwa index..."
    bwa index "$GENOME_FA"
fi

# 2. 循环所有样本
for SAMPLE in $(ls -d ${HOSTFREE_DIR}/*/); do
    SAMPLE_NAME=$(basename $SAMPLE)
    R1="${SAMPLE}/${SAMPLE_NAME}_filtered_1.fastq.gz"
    R2="${SAMPLE}/${SAMPLE_NAME}_filtered_2.fastq.gz"

    if [ ! -f "$R1" ]; then
        echo "Skip $SAMPLE_NAME, no fastq found"
        continue
    fi

    echo "Processing $SAMPLE_NAME ..."

    # 2.1 比对
    BAM="$OUT_DIR/${SAMPLE_NAME}.bam"
    bwa mem -t 16 "$GENOME_FA" "$R1" "$R2" | samtools view -bS - | samtools sort -o "$BAM"
    samtools index "$BAM"

    # 2.2 统计每个基因的 reads 数
    samtools idxstats "$BAM" > "$OUT_DIR/${SAMPLE_NAME}_idxstats.txt"

    # 2.3 统计总 reads 数
    TOTAL=$(($(zcat "$R1" | wc -l)/4 + $(zcat "$R2" | wc -l)/4))

    # 2.4 计算相对丰度
    awk -v total=$TOTAL -v sample=$SAMPLE_NAME 'BEGIN{OFS="\t"} 
        {if($1 != "*" ){gene=$1; count=$3; rel=(count/total); print sample, gene, count, total, rel}}
    ' "$OUT_DIR/${SAMPLE_NAME}_idxstats.txt" > "$OUT_DIR/${SAMPLE_NAME}_abundance.tsv"
done

