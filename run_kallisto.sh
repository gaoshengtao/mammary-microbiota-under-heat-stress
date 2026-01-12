#!/bin/bash

# 定义路径和参数
clean_data_dir="/home/gao/data/hsmilk_cp/hostfree"
output_base="/home/gao/data/hsmilk_cp/kallisto"
index_file="/home/gao/data/hsmilk_cp/megahit_out/prokka_annotation_90/metagG-index"
threads=64

# 遍历 cleandata 目录下的所有样本目录
for sample_dir in "$clean_data_dir"/*/; do
    # 提取样本名（目录名）
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"

    # 动态匹配带样本名的 FASTQ 文件（允许任意前缀）
    r1_files=("$sample_dir"/*"${sample}_filtered_1.fastq.gz")
    r2_files=("$sample_dir"/*"${sample}_filtered_1.fastq.gz")

    # 检查文件匹配情况
    if [[ ${#r1_files[@]} -eq 0 || ${#r2_files[@]} -eq 0 ]]; then
        echo "Error: No FASTQ files found for $sample"
        continue
    elif [[ ${#r1_files[@]} -gt 1 || ${#r2_files[@]} -gt 1 ]]; then
        echo "Error: Multiple R1/R2 files detected for $sample"
        echo "R1 candidates: ${r1_files[*]}"
        echo "R2 candidates: ${r2_files[*]}"
        continue
    fi

    # 创建输出目录
    output_dir="$output_base/$sample"
    mkdir -p "$output_dir"

    # 运行 kallisto quant
    kallisto quant \
        --threads "$threads" \
        -i "$index_file" \
        -o "$output_dir" \
        "${r1_files[0]}" \
        "${r2_files[0]}"

    echo "Finished processing $sample"
done

echo "All samples processed!"
