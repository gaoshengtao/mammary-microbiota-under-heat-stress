#!/bin/bash

# 配置路径
BOWTIE2_PATH="$HOME/apps/bowtie2-2.5.4-linux-x86_64/bowtie2"
HOST_INDEX="/home/gao/backup/database/kneaddata/cow/bowtie2/Bos_taurus"
CLEANDATA_DIR="/home/gao/data/hsmilk_cp/cleandata"
HOSTFREE_DIR="/home/gao/data/hsmilk_cp/hostfree"

# 遍历所有样本目录（根据cleandata下的子目录名）
for ID_DIR in ${CLEANDATA_DIR}/T*; do
  # 提取样本ID（目录名）
  ID=$(basename "${ID_DIR}")

  # 定义输入输出路径
  INPUT_DIR="${CLEANDATA_DIR}/${ID}"
  OUTPUT_DIR="${HOSTFREE_DIR}/${ID}"
  mkdir -p "${OUTPUT_DIR}"

  # 动态匹配输入文件名（支持任意前缀，但需符合 *-${ID}.R1.paired.fastq.gz 格式）
  R1_FILE=$(ls "${INPUT_DIR}/${ID}.R1.paired.fastq.gz" 2>/dev/null)
  R2_FILE=$(ls "${INPUT_DIR}/${ID}.R2.paired.fastq.gz" 2>/dev/null)

  # 检查文件是否存在且唯一
  if [ ! -f "${R1_FILE}" ] || [ ! -f "${R2_FILE}" ]; then
    echo "[ERROR] 输入文件缺失: ${ID}" >&2
    continue
  elif [ $(echo "${R1_FILE}" | wc -l) -gt 1 ] || [ $(echo "${R2_FILE}" | wc -l) -gt 1 ]; then
    echo "[ERROR] 存在多个匹配文件: ${ID}" >&2
    continue
  fi

  # 运行Bowtie2
  ${BOWTIE2_PATH} \
  -x "${HOST_INDEX}" \
  -1 "${R1_FILE}" \
  -2 "${R2_FILE}" \
  --threads 64 \
  --al-conc-gz "${OUTPUT_DIR}/${ID}_aligned_%.fastq.gz" \
  --un-conc-gz "${OUTPUT_DIR}/${ID}_filtered_%.fastq.gz" \
  -S /dev/null \
  --score-min L,0,-0.6 \
  --mp 8,3 \
  --rdg 15,5 \
  --rfg 15,5 \
  -L 20 \
  -D 10 -R 2

  echo "[SUCCESS] 已完成样本: ${ID}"
done
