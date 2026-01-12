#!/bin/bash
# RNA-seq 转录组分析流程
# 工具链: HISAT2 + StringTie + featureCounts
# 作者：ChatGPT（为高胜涛定制）
# 使用前提：
#   1. 已安装 hisat2, samtools, stringtie, subread (featureCounts)
#   2. 已构建好参考基因组索引 *.ht2
#   3. 有对应的 GTF 注释文件

# 路径设置（请修改为自己的）
REF_INDEX="/home/gao/ref/genome"          # HISAT2 索引前缀（不含 .1.ht2）
GTF_FILE="/home/gao/ref/genes.gtf"        # 注释文件
RAW_DIR="/home/gao/data/clean_data"       # 原始测序数据目录
OUT_DIR="/home/gao/data/rnaseq_out"       # 输出目录
THREADS=16

mkdir -p "$OUT_DIR"

# 遍历所有样本（假设命名为 sample1_1.fq.gz, sample1_2.fq.gz）
for R1 in "$RAW_DIR"/*_1.fq.gz; do
    sample=$(basename "$R1" _1.fq.gz)
    R2="${RAW_DIR}/${sample}_2.fq.gz"
    echo "=== 处理样本: $sample ==="

    # 1️⃣ 比对
    echo "[1/3] HISAT2 比对..."
    hisat2 -x "$REF_INDEX" -1 "$R1" -2 "$R2" -p $THREADS -S "$OUT_DIR/${sample}.sam"

    # 2️⃣ SAM 转 BAM 并排序
    echo "[2/3] SAMtools 排序..."
    samtools sort -@ $THREADS -o "$OUT_DIR/${sample}.sorted.bam" "$OUT_DIR/${sample}.sam"
    samtools index "$OUT_DIR/${sample}.sorted.bam"
    rm "$OUT_DIR/${sample}.sam"

    # 3️⃣ StringTie 定量
    echo "[3/3] StringTie 定量..."
    mkdir -p "$OUT_DIR/${sample}_stringtie"
    stringtie "$OUT_DIR/${sample}.sorted.bam" -G "$GTF_FILE" -o "$OUT_DIR/${sample}_stringtie/${sample}.gtf" -p $THREADS -e -B

    echo "=== 样本 $sample 完成 ==="
    echo
done

# 🔢 4️⃣ 汇总 featureCounts 结果（用于差异表达分析）
echo "[4] 运行 featureCounts 汇总基因表达矩阵..."
featureCounts -T $THREADS -a "$GTF_FILE" -o "$OUT_DIR/gene_counts.txt" "$OUT_DIR"/*.sorted.bam

echo "✅ 全部样本处理完成！"

