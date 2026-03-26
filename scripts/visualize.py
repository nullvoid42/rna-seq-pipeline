#!/usr/bin/env python3
"""
visualize.py
生成 PCA、火山图、热图
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    results_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    figures_dir = os.path.join(results_dir, "figures")
    tables_dir = os.path.join(results_dir, "tables")
    
    os.makedirs(figures_dir, exist_ok=True)
    
    print("=== 生成可视化 ===")
    
    # 加载数据
    with open(os.path.join(data_dir, "dds_object.pkl"), "rb") as f:
        obj = pickle.load(f)
    
    counts = obj["counts"]
    coldata = obj["coldata"]
    results = obj["results"]
    
    # log2 标准化 (VST-like)
    log_counts = np.log2(counts + 1)
    
    # =========================================================================
    # 1. PCA 图
    # =========================================================================
    print("生成 PCA 图...")
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(log_counts.T)
    
    colors = ["steelblue" if c == "control" else "firebrick" for c in coldata["condition"]]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, (x, y) in enumerate(pca_result):
        ax.scatter(x, y, c=colors[i], s=100, zorder=3)
        ax.annotate(coldata.index[i], (x, y), xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)", fontsize=12)
    ax.set_title("PCA - Airway RNA-seq (DEX vs Control)", fontsize=14)
    ax.legend(handles=[Patch(color="steelblue", label="Control"), Patch(color="firebrick", label="DEX")], loc="best")
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    pca_path = os.path.join(figures_dir, "pca_plot.png")
    plt.savefig(pca_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"PCA 图已保存: {pca_path}")
    
    # =========================================================================
    # 2. 火山图
    # =========================================================================
    print("生成火山图...")
    
    colors = []
    for _, row in results.iterrows():
        if row["padj"] < 0.05 and row["log2FoldChange"] > 1:
            colors.append("firebrick")
        elif row["padj"] < 0.05 and row["log2FoldChange"] < -1:
            colors.append("steelblue")
        else:
            colors.append("#999999")
    
    fig, ax = plt.subplots(figsize=(9, 7))
    
    # 预定义颜色
    color_map = {"firebrick": "#DC143C", "steelblue": "#4682B4", "#999999": "#999999"}
    rgba_colors = [color_map.get(c, c) for c in colors]
    
    ax.scatter(results["log2FoldChange"].values, -np.log10(results["padj"].values + 1e-300), 
               c=rgba_colors, alpha=0.6, s=15)
    
    ax.axhline(-np.log10(0.05), color="#666666", linestyle="--", linewidth=1)
    ax.axvline(1, color="#666666", linestyle="--", linewidth=1)
    ax.axvline(-1, color="#666666", linestyle="--", linewidth=1)
    
    # 标注 top genes
    top_genes = results.nsmallest(10, "padj")
    for _, row in top_genes.iterrows():
        if row["padj"] < 0.05 and abs(row["log2FoldChange"]) > 1.5:
            ax.annotate(row["symbol"], 
                       (row["log2FoldChange"], -np.log10(row["padj"] + 1e-300)),
                       fontsize=7, alpha=0.8)
    
    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10 (adjusted P-value)", fontsize=12)
    ax.set_title("Volcano Plot - DEX vs Control", fontsize=14)
    n_up = sum((results['padj'] < 0.05) & (results['log2FoldChange'] > 1))
    n_down = sum((results['padj'] < 0.05) & (results['log2FoldChange'] < -1))
    ax.legend(handles=[
        Patch(color="firebrick", label=f"Up ({n_up})"),
        Patch(color="steelblue", label=f"Down ({n_down})"),
        Patch(color="#999999", label="Not significant")
    ], loc="upper right")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    volcano_path = os.path.join(figures_dir, "volcano_plot.png")
    plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"火山图已保存: {volcano_path}")
    
    # =========================================================================
    # 3. 热图 (Top 20 差异基因)
    # =========================================================================
    print("生成热图...")
    
    sig = results[(results["padj"] < 0.05) & (abs(results["log2FoldChange"]) > 1)]
    top20 = sig.head(20)
    
    if len(top20) > 0:
        top20_counts = log_counts.loc[top20.index]
        
        # Z-score 标准化
        top20_z = (top20_counts - top20_counts.mean(axis=1).values.reshape(-1,1)) / top20_counts.std(axis=1).values.reshape(-1,1)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(top20_z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
        
        ax.set_xticks(range(len(top20_counts.columns)))
        ax.set_xticklabels(top20_counts.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20["symbol"].values, fontsize=8)
        
        # 列注释颜色
        for i, cond in enumerate(coldata["condition"]):
            ax.get_xticklabels()[i].set_color("firebrick" if cond == "dex" else "steelblue")
        
        plt.colorbar(im, ax=ax, label="Z-score")
        ax.set_title("Top 20 Differentially Expressed Genes", fontsize=14)
        
        plt.tight_layout()
        heatmap_path = os.path.join(figures_dir, "heatmap_plot.png")
        plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"热图已保存: {heatmap_path}")
    else:
        print("没有显著差异基因，跳过热图")
        heatmap_path = None
    
    print("\n=== 可视化完成 ===")
    print(f"PCA: {pca_path}")
    print(f"火山图: {volcano_path}")
    if heatmap_path:
        print(f"热图: {heatmap_path}")

if __name__ == "__main__":
    main()
