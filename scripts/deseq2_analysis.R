#!/usr/bin/env Rscript
# ============================================================================
# deseq2_analysis.R
# 使用 DESeq2 进行差异表达分析
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(tidyverse)
})

cat("=== DESeq2 差异表达分析 ===\n")

# 命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("需要输入: count_matrix.rds coldata.rds [output_dir]")
}

count_matrix_path <- args[1]
coldata_path <- args[2]
output_dir <- if (length(args) >= 3) args[3] else "results/tables"

# 加载数据
cat("加载数据...\n")
count_matrix <- readRDS(count_matrix_path)
coldata <- readRDS(coldata_path)

# 确保 coldata 的行名和 count matrix 的列名一致
coldata <- coldata[match(colnames(count_matrix), rownames(coldata)), ]
rownames(coldata) <- colnames(count_matrix)

cat("构建 DESeq2 数据集...\n")
# 创庺 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ dex
)

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep, ]
cat("过滤后基因数:", nrow(dds), "\n")

# 运行 DESeq2
cat("运行 DESeq2 分析...\n")
dds <- DESeq(dds)

# 获取结果
res <- results(dds, contrast = c("dex", "trt", "untrt"))
res <- na.omit(res)

# 添加基因 symbol
res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df <- as.data.frame(res) %>%
  arrange(padj)

# 保存完整结果
results_path <- file.path(output_dir, "deseq2_results.csv")
write.csv(res_df, results_path, quote = FALSE)
cat("结果已保存:", results_path, "\n")

# 筛选显著差异基因
sig_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
sig_path <- file.path(output_dir, "significant_genes.csv")
write.csv(sig_genes, sig_path, quote = FALSE)
cat("显著差异基因数:", nrow(sig_genes), "\n")

# 保存 dds 对象
dds_path <- file.path(dirname(count_matrix_path), "dds_object.rds")
saveRDS(dds, dds_path)

cat("\n=== 分析完成 ===\n")
cat("上调基因:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("下调基因:", sum(sig_genes$log2FoldChange < 0), "\n")
