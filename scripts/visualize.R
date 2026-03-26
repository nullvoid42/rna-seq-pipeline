#!/usr/bin/env Rscript
# ============================================================================
# visualize.R
# 生成 PCA、火山图、热图
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(tidyverse)
})

cat("=== 生成可视化图表 ===\n")

# 命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("需要输入: dds.rds deseq2_results.csv [figures_dir]")
}

dds_path <- args[1]
results_path <- args[2]
figures_dir <- if (length(args) >= 3) args[3] else "results/figures"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# 加载数据
cat("加载数据...\n")
dds <- readRDS(dds_path)
res_df <- read.csv(results_path, row.names = 1)

# ============================================================================
# 1. PCA 图
# ============================================================================
cat("生成 PCA 图...\n")

# 需要 vst 变换
vst <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vst, intgroup = "dex", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = dex)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), max.overlaps = 15) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA - Airway RNA-seq", color = "Treatment") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("steelblue", "firebrick"))

ggsave(
  file.path(figures_dir, "pca_plot.png"),
  pca_plot,
  width = 8, height = 6, dpi = 300
)

cat("PCA 图已保存\n")

# ============================================================================
# 2. 火山图
# ============================================================================
cat("生成火山图...\n")

res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not significant"
    ),
    label = ifelse(!is.na(symbol) & padj < 0.05 & abs(log2FoldChange) > 2,
                   symbol, NA)
  )

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_text_repel(aes(label = label), max.overlaps = 20, na.rm = TRUE) +
  scale_color_manual(
    values = c("Up" = "firebrick", "Down" = "steelblue", "Not significant" = "gray60")
  ) +
  labs(
    title = "Volcano Plot - DEX vs Control",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~P~value),
    color = "Significance"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(
  file.path(figures_dir, "volcano_plot.png"),
  volcano_plot,
  width = 9, height = 7, dpi = 300
)

cat("火山图已保存\n")

# ============================================================================
# 3. 热图 (Top 20 差异基因)
# ============================================================================
cat("生成热图...\n")

# 选择 top 20 padj 的显著基因
top_genes <- res_df %>%
  filter(!is.na(symbol), padj < 0.05) %>%
  arrange(padj) %>%
  head(20) %>%
  rownames()

# 提取 counts 并标准化
counts <- counts(dds, normalized = TRUE)[top_genes, ]
rownames(counts) <- res_df[top_genes, "symbol"]

# 注释数据
annot_col <- data.frame(
  Treatment = ifelse(colData(dds)$dex == "trt", "DEX", "Control"),
  row.names = colnames(dds)
)
annot_colors <- list(
  Treatment = c(DEX = "firebrick", Control = "steelblue")
)

# Z-score 标准化
counts_z <- t(scale(t(counts)))

pheatmap(
  counts_z,
  annotation_col = annot_col,
  annotation_colors = annot_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
  main = "Top 20 Differentially Expressed Genes",
  filename = file.path(figures_dir, "heatmap_plot.png"),
  width = 10,
  height = 8
)

cat("热图已保存\n")
cat("\n=== 可视化完成 ===\n")
cat("图片保存在:", figures_dir, "\n")
