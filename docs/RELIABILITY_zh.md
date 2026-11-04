# A 版增强流程

本版延续 Python CLI + Snakemake + R，支持 counts、已有 Salmon quant.sf 和 FASTQ→Salmon 三种输入。STAR/featureCounts 留待 B 版。

## 新增与修复

- 外部 tx2gene 全表检查空 ID、列数、一对多冲突；完全相同的重复映射可接受。R 导入再次检查。
- `analysis.covariate_types` 显式声明 `categorical` / `numeric`。contrast 的 factor 始终是分类变量，所以 0/1 等数字标签也能安全使用；声明为 numeric 的协变量必须全为有限数值。
- 当前仅支持加性模型；含 `*` / `:` 的交互公式会被拒绝，避免把参考水平的主效应误报成总体效应。
- 可选 `analysis.lfc_shrinkage: normal` 使用 DESeq2 自带正态先验。输出增加 `log2FoldChange_shrunk` / `lfcSE_shrunk` 和收缩 MA 图。原始 log2FC、Wald 检验、padj 和 UP/DOWN 筛选不变；默认 `none`。不将正态收缩宣称为适用于一切数据的最优方法。
- WGCNA 调用期间加载其要求的相关函数，退出时恢复包搜索环境；固定随机种子并限制其内部并行。固定 soft power 仍需真实项目诊断，不表示自动选择了最佳网络。
- `provenance.json` 记录输入/输出 SHA-256、字节数、有效配置、Python/工具版本和代码校验和；R sessionInfo 继续记录 R 包版本。FASTQ 输入的完整哈希会增加顺序读取开销。
- 关键结果、分析脚本及其 R 辅助文件纳入 Snakemake 依赖。`rnaseqflow run` 在复用结果前验证所有报告引用的分析文件；缺失或内容改变会强制重跑分析。
- R 每次在全新暂存目录中生成结果，成功后复制本次产物与摘要，避免阈值改变后把旧网络、旧富集图混入新报告；原目录中未引用的旧文件不会被主动删除。
- `rnaseqflow verify --results /path/to/results/analysis` 独立检查报告引用的分析结果；退出码 0 表示检查通过，1 表示缺失或损坏。
- HTML 报告直接显示 SVG/PNG 和上下调数量，并链接模型矩阵、配置和溯源记录。
- 环境清单及 doctor 加入 svglite；测试保证源码和 wheel 内 assets 相同。

## 配置示例

```yaml
analysis:
  level: intermediate
  design: "~ batch + age + condition"
  covariate_types:
    batch: categorical
    age: numeric
    condition: categorical
  lfc_shrinkage: normal
  contrasts:
    - {name: treated_vs_control, factor: condition, numerator: treated, denominator: control}
```

明确类型优先；未声明且非对比因素的其他列保持自动推断，兼容旧配置。不要以整数编码批次却不声明为 categorical。数值连续变量不能作为当前命名分类对比的 factor。

## 验收与边界

```bash
python -m pytest -q
Rscript tests/integration/r_checks.R
RNASEQFLOW_INTEGRATION=1 python -m pytest tests/integration -q
```

默认 Python 测试不要求 R；集成测试明确为 opt-in，并在 CI 完整环境任务中执行。真实轻量级示例及独立数值对照见 [pasilla](../examples/pasilla/README.md)。CI 文件已加入，远端 GitHub Actions 是否通过需在 PR 中观察，本地验证不代表远端执行。

`verify` 校验分析结果文件，不是数据真实性/生物学正确性的证明，也不是数字签名。原始输入和源代码的校验和用于审计；Snakemake 的输入/代码文件依赖用于调度。HTML 索引由 report 规则维护，不对 manifest 自身或 HTML 做循环哈希。完整性检查不能替代 FASTQ、Salmon、基因集覆盖率和研究设计的专业审阅。不要在同一个项目目录并发运行不同配置；更改层级或设计优先创建独立项目，旧的未引用文件不代表当前结果。

本次真实数据验收覆盖 counts/intermediate；WGCNA 使用合成数据作软件回归，不宣称获得真实网络发现。真实 FASTQ/Salmon、真实通路富集和剪接仍需后续数据验证。
