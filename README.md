# RNA-seq 差异表达分析全流程

使用 DESeq2 对 Bioconductor 的 airway 数据集进行差异表达分析。

## 数据集

**airway**: 人类气道平滑肌细胞 RNA-seq 数据
- 8 个样本：4 个 DEX 处理 vs 4 个对照
- 差异表达分析：DEX 处理 vs 无处理

## 工具链

- **R 4.3+** + Bioconductor
- **DESeq2**: 差异表达分析
- **ggplot2/ggrepel**: 可视化
- **pheatmap**: 热图
- **Snakemake**: 流程管理（可选）

## 目录结构

```
rna-seq-pipeline/
├── Snakefile              # Snakemake 流程定义
├── config.yaml            # 配置文件
├── envs/
│   └── deseq2.yaml        # Conda 环境
├── scripts/
│   ├── download_data.R    # 数据下载
│   ├── deseq2_analysis.R  # 差异表达分析
│   └── visualize.R        # 可视化
└── results/
    ├── figures/           # 图片输出
    └── tables/            # 表格输出
```

## 快速开始

### 方法 1: 使用 conda 环境（推荐）

```bash
# 1. 安装 miniforge (如果没有 conda)
brew install --cask miniforge
source ~/.zshrc
mamba init

# 2. 创建环境
cd ~/bioinformatics/rna-seq-pipeline
mamba env create -f envs/deseq2.yaml

# 3. 激活环境
mamba activate deseq2_pipeline

# 4. 下载数据
Rscript scripts/download_data.R data

# 5. 运行分析
Rscript scripts/deseq2_analysis.R data/airway_counts.rds data/airway_coldata.rds results/tables

# 6. 生成可视化
Rscript scripts/visualize.R data/dds_object.rds results/tables/deseq2_results.csv results/figures
```

### 方法 2: 直接安装 R 包

```bash
# 安装 R
brew install --cask r-base

# 启动 R 并安装包
R

install.packages(c("tidyverse", "ggrepel", "pheatmap"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "airway", "org.Hs.eg.db"))
```

### 方法 3: 使用 Docker

```bash
# 拉取 Bioconductor Docker 镜像
docker pull bioconductor/bioconductor_docker:latest

# 运行分析
docker run -v $(pwd):/workspace bioconductor/bioconductor_docker:latest Rscript scripts/download_data.R /workspace/data
```

## 输出结果

### 表格
- `deseq2_results.csv`: 所有基因的差异表达结果
- `significant_genes.csv`: 显著差异基因（padj < 0.05, |log2FC| > 1）

### 图片
- `pca_plot.png`: PCA 图
- `volcano_plot.png`: 火山图
- `heatmap_plot.png`: 差异基因热图

## 结果解读

### PCA 图
- 展示样本间的整体差异
- 如果处理组和对照组在 PCA 空间分离，说明有系统性差异

### 火山图
- X轴: log2 Fold Change（表达量变化倍数）
- Y轴: -log10(adjusted p-value)（统计显著性）
- 右上/左上: 显著上调/下调基因

### 热图
- 展示 Top 20 差异基因的表达模式
- 红色=高表达，蓝色=低表达

## 常见问题

**Q: 包安装失败？**
A: 确保 R 版本 >= 4.3，使用 BiocManager 安装 Bioconductor 包

**Q: 内存不足？**
A: airway 数据集较小，一般不需要特殊配置

**Q: 没有 conda？**
A: 使用 Homebrew 安装: `brew install --cask miniforge`
