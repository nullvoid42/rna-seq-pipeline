#!/usr/bin/env Rscript
# ============================================================================
# download_data.R
# 从 Bioconductor 下载 airway 数据集
# Dataset: 人类气道平滑肌细胞，DEX 处理 vs 对照
# ============================================================================

suppressPackageStartupMessages({
  library(BiocManager)
  library(airway)
  library(DESeq2)
})

# 输出目录
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  data_dir <- "data"
} else {
  data_dir <- args[1]
}

cat("=== 下载 airway 数据集 ===\n")

# airway 数据集已经包含在包里面，不需要额外下载
# 它是一个 RangedSummarizedExperiment 对象
data("airway")
cat("airway 数据集加载成功！\n")

# 提取 count matrix 和 coldata
count_matrix <- assay(airway)
coldata <- colData(airway)

cat("Count matrix 维度:", dim(count_matrix), "\n")
cat("样本数:", ncol(count_matrix), "\n")
cat("基因数:", nrow(count_matrix), "\n")
cat("\n样本信息:\n")
print(coldata)

# 保存为 RDS
count_matrix_path <- file.path(data_dir, "airway_counts.rds")
coldata_path <- file.path(data_dir, "airway_coldata.rds")

saveRDS(count_matrix, count_matrix_path)
saveRDS(as.data.frame(coldata), coldata_path)

cat("\n=== 数据保存完成 ===\n")
cat("Count matrix:", count_matrix_path, "\n")
cat("Sample info:", coldata_path, "\n")
