#!/bin/bash

# 配置参数
input_dir="/home/gao/backup/heatstress_milk_microbe/MJ20241223222-40milk/rawdata/"
metaphlan_DIR="/home/gao/data/hsmilk_cp/metaphlan_rawdata"
METAPHLAN_DB="/home/gao/data/database/metaphlan"


# 智能提取样本ID函数
extract_sample_id() {
    local filename=$(basename "$1")
    # 匹配模式：最后一个连字符后到第一个点前的内容（如TAH01N_2HZ）
    local sample_id=$(echo "$filename" | sed -E 's/.*-([^.]+)\..*/\1/')
    # 二次校验格式：TAH/TAP开头 + 数字 + 字母 + _ + 数字 + HZ
    if [[ "$sample_id" =~ ^TA[HP][0-9]+[A-Z]_[0-9]+HZ$ ]]; then
        echo "$sample_id"
    else
        echo "INVALID"
    fi
}

# 遍历所有R1文件
find "${input_dir}" -name "*.R1.raw.fastq.gz" | while read r1_file; do
    # 提取并验证样本ID
    sample_id=$(extract_sample_id "${r1_file}")
    if [ "$sample_id" == "INVALID" ]; then
        echo "[ERROR] 文件名格式错误: ${r1_file}"
        continue
    fi

    # 构建R2文件路径
    r2_file=$(echo "${r1_file}" | sed 's/\.R1\.raw\.fastq\.gz/.R2.raw.fastq.gz/')
    
    # 检查R2文件是否存在
    if [ ! -f "${r2_file}" ]; then
        echo "[ERROR] 找不到R2文件: ${r2_file}，跳过样本 ${sample_id}"
        continue
    fi

    # 创建输出目录
    output_dir="${metaphlan_DIR}/${sample_id}"
    mkdir -p "${output_dir}"

  #运行MetaPhlAn
  metaphlan "${r1_file},${r2_file}" \
  -t rel_ab_w_read_stats \
  --bowtie2out "${output_dir}/${sample_id}.bowtie2.bz2" \
  --nproc 360 \
  --input_type fastq \
  -o "${output_dir}/${sample_id}_rawdata.txt" \
  --bowtie2db "${METAPHLAN_DB}"

    echo "[$(date)] 样本 ${sample_id} 处理完成"
done

echo "所有样本处理完成！日志文件保存在：${output_base}/logs"
