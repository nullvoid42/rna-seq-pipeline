# 🤝 Contributing to RNA-seq Pipeline

欢迎贡献代码！请遵循以下指南。

## 如何贡献

### 1. Fork & Clone
```bash
git clone https://github.com/YOUR_USERNAME/rna-seq-pipeline.git
cd rna-seq-pipeline
```

### 2. 创建分支
```bash
git checkout -b feature/your-feature-name
```

### 3. 开发
- 添加新功能或修复 bug
- 确保代码有注释
- 添加测试（如果有）

### 4. 提交
```bash
git commit -m "Add: your feature description"
git push origin feature/your-feature-name
```

### 5. Pull Request
在 GitHub 上创建 PR，描述你的改动。

## 项目结构
```
scripts/
  download_data.py    # 数据下载/生成
  deseq2_analysis.py  # 差异表达分析
  visualize.py        # 可视化
```

## 规范
- Python 代码遵循 PEP 8
- 所有脚本在 `scripts/` 目录下
- 结果输出到 `results/` 目录
- 不要提交大文件（>10MB）到 git

## 问题反馈
- 使用 GitHub Issues 报告 bug
- 描述清楚环境、步骤和预期行为
