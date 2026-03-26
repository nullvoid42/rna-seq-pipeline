#!/usr/bin/env python3
"""
deseq2_analysis.py
使用 Python 进行差异表达分析
模拟 DESeq2 的统计方法
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats

def deseq2_like_analysis(counts, condition):
    """
    简化版 DESeq2 分析流程
    使用对数变换 + t-test 作为差异表达检验
    """
    control = counts.loc[:, condition == "control"]
    treatment = counts.loc[:, condition == "dex"]
    
    # 计算 mean 和 log2FC
    control_mean = control.mean(axis=1) + 1
    treatment_mean = treatment.mean(axis=1) + 1
    log2fc = np.log2(treatment_mean / control_mean)
    
    # 计算 p-value (Welch's t-test)
    pvals = []
    for gene in counts.index:
        c_vals = control.loc[gene].values + 1
        t_vals = treatment.loc[gene].values + 1
        if np.var(c_vals) == 0 and np.var(t_vals) == 0:
            pvals.append(1.0)
        else:
            _, p = stats.ttest_ind(np.log2(t_vals), np.log2(c_vals), equal_var=False)
            pvals.append(p)
    
    pvals = np.array(pvals)
    
    # BH 校正
    padj = stats.false_discovery_control(pvals, method='bh')
    
    results = pd.DataFrame({
        "baseMean": (control_mean + treatment_mean) / 2,
        "log2FoldChange": log2fc,
        "pvalue": pvals,
        "padj": padj
    }, index=counts.index)
    
    return results

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.dirname(__file__)), "results", "tables")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== 差异表达分析 ===")
    
    # 加载数据
    counts = pd.read_csv(os.path.join(data_dir, "airway_counts.tsv"), sep="\t", index_col=0)
    coldata = pd.read_csv(os.path.join(data_dir, "airway_coldata.tsv"), sep="\t", index_col=0)
    
    condition = coldata["condition"].values
    
    print(f"样本数: {counts.shape[1]}, 基因数: {counts.shape[0]}")
    
    # 运行分析
    results = deseq2_like_analysis(counts, condition)
    
    # 添加基因符号（模拟）
    results["symbol"] = [f"GENE_{i}" for i in range(len(results))]
    
    # 按 padj 排序
    results = results.sort_values("padj")
    
    # 保存完整结果
    results_file = os.path.join(output_dir, "deseq2_results.csv")
    results.to_csv(results_file)
    print(f"结果已保存: {results_file}")
    
    # 显著基因
    sig = results[(results["padj"] < 0.05) & (np.abs(results["log2FoldChange"]) > 1)]
    sig_file = os.path.join(output_dir, "significant_genes.csv")
    sig.to_csv(sig_file)
    print(f"显著差异基因数: {len(sig)} (上调: {sum(sig['log2FoldChange'] > 0)}, 下调: {sum(sig['log2FoldChange'] < 0)})")
    
    # 保存 dds-like object (即 counts + coldata, 用于可视化)
    import pickle
    dds_file = os.path.join(data_dir, "dds_object.pkl")
    with open(dds_file, "wb") as f:
        pickle.dump({"counts": counts, "coldata": coldata, "results": results}, f)
    print(f"对象已保存: {dds_file}")
    
    print("\n=== 分析完成 ===")

if __name__ == "__main__":
    main()
