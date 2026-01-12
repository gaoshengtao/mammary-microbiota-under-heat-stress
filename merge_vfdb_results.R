#!/usr/bin/env Rscript

library(dplyr)
library(readr)

# 输入和输出路径
in_dir <- "/home/gao/data/hsmilk_cp/vfdb/mapping_results"
out_file <- "/home/gao/data/hsmilk_cp/vfdb/vfdb_abundance_matrix.tsv"

# 找到所有 abundance 文件
files <- list.files(in_dir, pattern = "_abundance.tsv$", full.names = TRUE)

# 读入并合并
abundance_list <- lapply(files, function(f) {
  df <- read_tsv(f,
                 col_names = c("Sample", "GeneID", "Mapped_Reads", "Total_Reads", "Relative_Abundance"),
                 col_types = cols(
                   Sample = col_character(),
                   GeneID = col_character(),
                   Mapped_Reads = col_double(),
                   Total_Reads = col_double(),
                   Relative_Abundance = col_double()
                 )
  )
  sample_name <- unique(df$Sample)
  df %>%
    select(GeneID, Relative_Abundance) %>%
    rename(!!sample_name := Relative_Abundance)
})

# 合并成矩阵
abundance_matrix <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), abundance_list)

# NA 填 0
abundance_matrix[is.na(abundance_matrix)] <- 0

# 保存结果
write_tsv(abundance_matrix, out_file)

cat("✅ 合并完成，结果输出到: ", out_file, "\n")
