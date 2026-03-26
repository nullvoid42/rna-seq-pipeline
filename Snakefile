# RNA-seq Differential Expression Analysis Pipeline
# 使用 DESeq2 对 airway 数据集进行差异表达分析
# 数据来源：Bioconductor airway 包

import os

# 配置
configfile: "config.yaml"

# 输出目录
RESULTS_DIR = config["results_dir"]
DATA_DIR = config["data_dir"]
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")

# 确保输出目录存在
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# 规则
rule all:
    input:
        os.path.join(FIGURES_DIR, "pca_plot.png"),
        os.path.join(FIGURES_DIR, "volcano_plot.png"),
        os.path.join(FIGURES_DIR, "heatmap_plot.png"),
        os.path.join(TABLES_DIR, "deseq2_results.csv"),
        os.path.join(TABLES_DIR, "significant_genes.csv")

rule download_data:
    output:
        count_matrix = os.path.join(DATA_DIR, "airway_counts.rds"),
        coldata = os.path.join(DATA_DIR, "airway_coldata.rds")
    conda:
        "envs/deseq2.yaml"
    script:
        "scripts/download_data.R"

rule deseq2_analysis:
    input:
        count_matrix = rules.download_data.output.count_matrix,
        coldata = rules.download_data.output.coldata
    output:
        dds = os.path.join(DATA_DIR, "dds_object.rds"),
        results = os.path.join(TABLES_DIR, "deseq2_results.csv"),
        sig_genes = os.path.join(TABLES_DIR, "significant_genes.csv")
    conda:
        "envs/deseq2.yaml"
    script:
        "scripts/deseq2_analysis.R"

rule visualize:
    input:
        dds = rules.deseq2_analysis.output.dds,
        results = rules.deseq2_analysis.output.results
    output:
        pca = os.path.join(FIGURES_DIR, "pca_plot.png"),
        volcano = os.path.join(FIGURES_DIR, "volcano_plot.png"),
        heatmap = os.path.join(FIGURES_DIR, "heatmap_plot.png")
    conda:
        "envs/deseq2.yaml"
    script:
        "scripts/visualize.R"
