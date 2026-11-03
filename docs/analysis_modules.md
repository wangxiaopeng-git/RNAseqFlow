# Analysis modules and staged roadmap

## Implemented modules (v0.1.0 preview)

| Stage | Inputs | Primary executable | Key outputs | Scientific checks |
|---|---|---|---|---|
| Input validation | YAML + sample sheet | `rnaseqflow validate` | preflight status | unique IDs, integer counts, design contrasts, required files |
| Reference helper | GTF/GFF3 | `rnaseqflow tx2gene` | tx2gene.tsv | transcript-to-gene cardinality; no silent ID stripping |
| Raw QC | FASTQ | FastQC | sample HTML/ZIP | base quality / adapter reads / over-representation |
| Trim | FASTQ | fastp | cleaned FASTQ, JSON, HTML | pair preservation; per-sample adapters/quality |
| Aggregate QC | fastp/FastQC | MultiQC | one HTML | across-sample QC review |
| Transcript quantification | transcriptome + FASTQ | Salmon | quant.sf for every sample | annotation version; reference completeness |
| Import Salmon | quant.sf + tx2gene | tximport | counts with DESeq2 offsets; TPM table | no TPM in DESeq2 count constructor |
| Gene-level analysis | integer counts or tximport | DESeq2 | normalized counts, size factors | counts QC, low-expression filter |
| Sample QC | normalized expression | ggplot2 / pheatmap | PCA, correlations, sample depth | sequencing batches, outliers |
| Differential expression | count model & design | DESeq2 | per-contrast full and selected tables | controls and replicates, FDR, full rank |
| DEG plots | DE results | ggplot2 / DESeq2 / pheatmap | volcano/MA/heatmap SVG | heatmap Z-score is not absolute expression |
| Gene-set ORA | custom GMT + DE list | hypergeometric | ora.tsv, SVG | eligible tested background, BH adjustment |
| Preranked GSEA | custom GMT + Wald stats | fgseaMultilevel | gsea.tsv | consistent IDs and adequate pathway sizes |
| Signed WGCNA | transformed counts + traits | WGCNA | modules, eigengenes, trait FDR | >= 20 samples by default; soft threshold diagnostic needed |
| Cross-contrast comparison | selected DEG lists | R Jaccard | DEG set overlap table/heatmap | >= 2 contrasts; no causal interpretation |
| Provenance | config + environment | R + JSON | sessionInfo.txt, effective_config.yaml, summary.json | module availability and reported versions |
| Reporting | generated outputs | Python | offline index.html | links to original tables and editable images |

## Basic analysis

Use to check whether a dataset is usable. For raw FASTQ this includes QC, trimming, quantification. For all modes it includes gene expression summarization and sample-level QC. It does **not** run DESeq2 significance tests or pathway inference. Useful for exploratory design checks or counts with no replicated treatment.

## Intermediate analysis

Requires an appropriate replicated experiment. Run batch-aware negative-binomial differential expression (DESeq2), export effect estimates, raw p-values, BH FDR, and numerical quality-control metrics. Build custom GO/KEGG and other annotation-based GMT pathways for non-model organisms. A pathway with nominal P < 0.05 but BH FDR >= 0.05 is **not significant**.

## Advanced analysis

Adds multi-factor design support through the R design formula (with conventional named contrasts) and WGCNA module discovery and module–trait association, when sample size and expression matrix are sufficient. WGCNA is exploratory, has potentially unstable modules for small n, and should not be used uncritically with 20 samples. A fixed beta is not sufficient to claim scale-free topology; sensitivity and diagnostics are needed.

## Future modules — not part of v0.1

1. Decoy-aware Salmon indexing; optional STAR/HISAT2 + featureCounts branch, strandedness and mapping QC.
2. Transcript-level differential usage (DRIMSeq / stageR), splicing (rMATS/MAJIQ), and isoform switching.
3. Time-course differential expression; nested repeated-measures designs; robust effect-size threshold tests and shrinkage diagnostics.
4. QC enhancements: gene-body coverage, 5'/3' bias, rRNA fraction, contamination, ribosomal mapping, strandedness, sample fingerprinting.
5. Functional annotation interoperability: GO OBO gene mappings; KEGG KOs, species-specific background generation; eggNOG support.
6. Cell-type deconvolution, multi-study RNA-seq meta-analysis, cross-species ortholog conversion, multiple-condition enrichment comparison.
7. Exploratory differential coexpression, network topology and hub statistics with replication/permutation checks.
8. Optional cluster/HPC profiles (Slurm), container images, pinned package lock files, versioned reference indexes, runtime and storage benchmarks.

No future module is advertised as implemented. Additional domains (small RNA, single-cell, spatial, metatranscriptomics) require their own validated methods rather than blindly reusing bulk RNA-seq thresholds.
