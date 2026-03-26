#!/usr/bin/env python3
"""
download_data.py
生成模拟的 airway RNA-seq count matrix
8样本：4 control, 4 DEX treatment
"""

import os
import pandas as pd
import numpy as np

# 输出目录
script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(script_dir, "data")
os.makedirs(data_dir, exist_ok=True)

print("=== 生成模拟 airway RNA-seq 数据 ===")

np.random.seed(42)
n_genes = 15000
n_samples = 8

# 基因 ID
gene_ids = [f"ENSG{str(i).zfill(10)}" for i in range(1, n_genes + 1)]

# 样本 ID
sample_ids = [f"SRR{i}" for i in range(1, 9)]

# 生成基础表达量 (log-normal)
base_means = np.random.lognormal(4, 1.5, n_genes)

# 生成 counts (negative binomial)
counts = np.zeros((n_genes, n_samples), dtype=int)
for i in range(n_genes):
    for j in range(n_samples):
        mu = base_means[i]
        # Convert to nb params: var = mu + mu^2/theta
        theta = max(mu, 5)  # dispersion
        var = mu + mu**2 / theta
        p = mu / var
        n = mu**2 / (var - mu)
        counts[i, j] = np.random.negative_binomial(n, p)

# 转为 DataFrame
counts_df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)

# 差异表达基因: 前200个基因有差异
# 前100上调，后100下调
for i in range(100):
    counts_df.iloc[i, 4:] = (counts_df.iloc[i, 4:] * 2.5).astype(int)
for i in range(100, 200):
    counts_df.iloc[i, 4:] = (counts_df.iloc[i, 4:] * 0.35).astype(int)

# 样本信息
coldata = pd.DataFrame({
    "sample": sample_ids,
    "condition": ["control"] * 4 + ["dex"] * 4
}, index=sample_ids)

# 保存
counts_path = os.path.join(data_dir, "airway_counts.tsv")
coldata_path = os.path.join(data_dir, "airway_coldata.tsv")
counts_df.to_csv(counts_path, sep="\t")
coldata.to_csv(coldata_path, sep="\t")

print(f"Count matrix: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
print(f"\n样本信息:")
print(coldata)
print(f"\n数据已保存:")
print(f"  {counts_path}")
print(f"  {coldata_path}")
print("=== 完成 ===")
