# RNAseqFlow 中文快速开始

**版本**：v0.1.0 开发预览。已通过 Python 输入校验与自动测试，尚未在此开发环境完成外部 R/Snakemake 与真实 FASTQ 端到端验证。

## 安装

建议 Ubuntu/Linux，并安装 Mamba/Conda：

```bash
mamba env create -f envs/full.yaml
conda activate rnaseqflow
pip install -e .
```

## 运行自带示例

```bash
rnaseqflow validate -c config/config.yaml
rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --mem-mb 32000 --dry-run
rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --mem-mb 32000
```

示例包含 24 个**虚构样本**、200 个**人工基因**与若干演示基因集，不代表真实的GO/KEGG注释或真实研究。

## 三个数据入口

- `counts`：`samples.tsv` + 原始整数基因 counts 矩阵。只做统计分析时最方便。
- `salmon`：已存在每个样本的 `quant.sf` + `tx2gene.tsv`，跳过FASTQ质控、索引与定量。
- `fastq`：原始FASTQ + 转录本参考序列（transcript FASTA）+ `tx2gene.tsv`，自动进行FastQC、fastp、Salmon和下游分析。

`tx2gene.tsv` 可以从带有转录本信息的 GTF/GFF3 自动提取：

```bash
rnaseqflow tx2gene --annotation genes.gff3 --output tx2gene.tsv
```

修改 `config/config.yaml` 的路径、输入模式、分组与对比组。设计模型默认 `~ batch + condition`，如没有批次列或无需控制批次，改为 `~ condition`。不要仅删除样本表中的 batch 列却保留模型变量。

## 分析层级

- `--level basic`：原始数据处理（若提供FASTQ）、定量和样本QC。
- `--level intermediate`：增加DESeq2差异分析和可选GMT通路富集。
- `--level advanced`：增加探索性的WGCNA与性状相关分析。WGCNA默认至少20个样本；不足时记录SKIPPED。

可以运行：

```bash
mkdir -p logs
nohup rnaseqflow run -c config/config.yaml --level advanced \
  --cores 24 --jobs 6 --mem-mb 90000 > logs/rna.log 2>&1 &
```

`--cores 24` 表示最大CPU线程数，`--jobs 6` 表示最多同时执行6个Snakemake任务，两者均实际生效，`--mem-mb` 约束调度声明的全局内存。

结果入口：`results/analysis/report/index.html`，常用原始表格：`normalized_counts.tsv`、`contrasts/<比较名>/deg_all.tsv`、`deg_selected.tsv`、`enrichment/ora.tsv`与`gsea.tsv`。

## 重要注意

- **TPM/FPKM不允许作为DESeq2的原始输入 counts**。
- GO/KEGG需要真实、可靠、与基因ID对应的功能注释集合（GMT）；示例文件仅为人工基因集。
- Padj(FDR)<0.05才可报告显著富集；名义P值需明确说明。
- 要求合理的生物学重复数，避免处理组与批次完全共线。
- 高级模式的共表达模块只是相关证据，不等于基因调控因果机制。
- 正式用于SCI论文之前，请先核对软件版本、参数、样本配对、比对/定量情况，并与独立分析结果对照验证。

## 发布到独立 GitHub 仓库

由于当前 ChatGPT 的 GitHub 接口不能创建仓库，项目包提供了可选的 GitHub CLI 发布脚本。在你自己的电脑安装并登录 `gh` 后，从项目目录运行：

```bash
gh auth login
bash scripts/publish_to_github.sh
```

默认创建 `wangxiaopeng-git/RNAseqFlow` 公共仓库。如需先私有发布，可设 `GH_VISIBILITY=private`；如需另一个仓库名，可设 `GH_REPO=your-new-name`。脚本发现同名远程仓库后会拒绝覆盖。
