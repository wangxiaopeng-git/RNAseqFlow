# RNAseqFlow

**Modular bulk RNA-seq analysis from raw FASTQ, Salmon `quant.sf`, or a gene-count matrix.**

RNAseqFlow provides a unified command-line interface for RNA-seq processing, exploratory QC, differential expression, pathway analysis, and coexpression networks. It does **not** assume a human or mouse annotation database. Run local reference transcriptomes and custom gene-set files for non-model organisms.

**Status: v0.1.0 development preview.** Python CLI/validation and synthetic-input tests have been exercised. R/Bioconductor execution, Snakemake DAGs, real sequencing data, and Conda environment resolution **must be verified on an RNA-seq-capable Linux installation before publication or production use**. No real scientific result has been validated by this release.

## Three levels

| Module | Basic | Intermediate | Advanced | Implementation |
|---|:---:|:---:|:---:|---|
| FASTQ validation and per-sample fastp | ✓ | ✓ | ✓ | Executable workflow (FASTQ mode) |
| Raw FastQC + MultiQC | ✓ | ✓ | ✓ | Executable workflow (FASTQ mode) |
| Salmon transcript quantification | ✓ | ✓ | ✓ | Executable workflow (FASTQ mode) |
| Import existing `quant.sf` with `tximport` | ✓ | ✓ | ✓ | Implemented |
| Input gene-count matrix | ✓ | ✓ | ✓ | Implemented; raw nonnegative integers only |
| Normalized counts / TPM (when available) | ✓ | ✓ | ✓ | Implemented; TPM is not used as DESeq2 counts |
| Library sizes, PCA, sample correlations | ✓ | ✓ | ✓ | Implemented |
| Full-rank design checking, multi-factor DESeq2 | — | ✓ | ✓ | Implemented: standard designs + named contrasts |
| DEG tables, volcano, MA, top-DEG heatmap | — | ✓ | ✓ | Implemented |
| GMT-based gene-set ORA (BH correction) | — | ✓ | ✓ | Implemented, optional |
| Preranked GSEA via `fgseaMultilevel` | — | ✓ | ✓ | Implemented, optional |
| Signed WGCNA and module-trait associations | — | — | ✓ | Implemented, optional, min 20 samples by default |
| Cross-contrast DEG Jaccard similarity | — | — | ✓ | Implemented for two or more contrasts |
| Sample-level metadata, R session log, HTML index | ✓ | ✓ | ✓ | Implemented |
| STAR / HISAT2 + featureCounts alignment | — | — | — | Planned, not implemented |
| rMATS alternative splicing / isoform switching | — | — | — | Planned, not implemented |
| DRIMSeq / stageR differential transcript usage | — | — | — | Planned, not implemented |
| Fusion, circRNA, ASE, time-series, multi-omics | — | — | — | Planned, not implemented |

The input routes (counts / Salmon / FASTQ) are **alternatives**, not separate analysis levels. Higher levels include the lower-level analytical stages that apply to the chosen route.

## Installation (Linux)

Install [Mamba](https://mamba.readthedocs.io/) or an appropriate Conda-compatible solver, then:

```bash
mamba env create -f envs/full.yaml
conda activate rnaseqflow
python -m pip install -e .
rnaseqflow --help
rnaseqflow doctor
```

The workflow uses command-line `snakemake`, `fastp`, `fastqc`, `multiqc`, `salmon`, and `Rscript` from the active environment. The all-in-one environment specification is in `envs/full.yaml`; Conda solver/package compatibility depends on platform and repository updates. `envs/python.yaml` and `envs/r-analysis.yaml` are split alternatives for advanced/manual deployment, not a complete one-command execution environment.

## First run: small synthetic count data

The repository ships 24 fictional samples with 200 artificial genes and five fake gene sets. They test file and result handling, **not biological correctness or realistic sequencing performance**.

```bash
# Run in the repository root
rnaseqflow validate -c config/config.yaml
rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --mem-mb 32000 --dry-run
rnaseqflow run -c config/config.yaml --cores 8 --jobs 2 --mem-mb 32000
```

Results are written to `results/analysis/`, with the landing page at `results/analysis/report/index.html`. `--jobs` is an actual Snakemake resource limit (`concurrent_jobs`) that applies to **every** rule; `--cores` separately limits total threads. The global `--mem-mb` resource prevents jobs exceeding the declared per-rule memory requests, but memory requirements still need adjustment for each real dataset.

For a separate project with the same synthetic example:

```bash
rnaseqflow init --dest /path/to/my_rnaseq_project
rnaseqflow validate -c /path/to/my_rnaseq_project/config/config.yaml
rnaseqflow run -c /path/to/my_rnaseq_project/config/config.yaml --level basic --cores 8 --jobs 2
```

## Using your own data

Edit the `config/config.yaml` and sample sheet before running:

**1. Raw counts (quickest):** set `project.input_mode: counts`, `project.counts_matrix` to a tab-delimited matrix (`gene_id` as the first column, one column per sample), and `project.samples` to a `samples.tsv`. Counts must be **raw integer gene counts**, not TPM, FPKM, CPM, transformed values, or Excel-rounded normalized counts.

**2. Existing Salmon:** set `project.input_mode: salmon`; include `sample_id`, `condition` (and design columns), and a `quant_dir` containing `quant.sf` for every sample. Provide a `tx2gene.tsv` with `transcript_id` / `gene_id` header. Both transcript IDs and gene IDs must match the same annotation/reference version.

**3. FASTQ:** set `project.input_mode: fastq`, `project.transcript_fasta` to the transcriptome FASTA, `project.tx2gene` to the gene map, and `fastq_1`/`fastq_2` in the sample sheet (leave fastq_2 blank for single-end). One project can use single-end or paired-end reads, but not a mixture. Validated path: raw FASTQ → fastp → Salmon → tximport → analysis, with FastQC/MultiQC. This path currently requires a preassembled **transcriptome FASTA**, not only a genomic FASTA.

For an annotation GTF/GFF3, create a mapping without an online database:

```bash
rnaseqflow tx2gene --annotation reference/annotation.gtf --output reference/tx2gene.tsv
# GFF3 example: --annotation reference/genes.gff3
```

GFF3 extraction expects transcript features (e.g. `mRNA`) with `ID` and `Parent` pointing to the gene. Some annotation flavors require manual normalization. Transcript IDs from mapping and Salmon index must be checked exactly; versions/prefixes are **not silently removed**.

For custom pathways, supply `project.gene_sets_gmt`: each row is `set_name<TAB>description<TAB>gene1<TAB>gene2...`. The demo sets labeled `TOY_` are **not true GO or KEGG**. Enrichments use genes detected and retained in the analysis as the background; IDs must be consistent across counts and gene sets. Obtain species-specific GO/KEGG mappings and convert to GMT before reporting biological pathway names.

## Minimal sample metadata

```tsv
sample_id	condition	batch	fastq_1	fastq_2	quant_dir
S01	control	B1	/path/S01_R1.fastq.gz	/path/S01_R2.fastq.gz	
S02	control	B2	/path/S02_R1.fastq.gz	/path/S02_R2.fastq.gz	
S03	treated	B1	/path/S03_R1.fastq.gz	/path/S03_R2.fastq.gz	
S04	treated	B2	/path/S04_R1.fastq.gz	/path/S04_R2.fastq.gz	
```

For counts input, FASTQ and quant paths can be blank. `analysis.design: "~ batch + condition"` controls the full model, and contrast fields identify the factor, numerator, and denominator. Non-full-rank designs (such as a treatment perfectly confounded with batch) stop with an error rather than returning misleading DEGs. Standard condition contrasts are supported. Advanced interaction designs may require explicit contrast coefficient handling that is not yet exposed by the CLI; do not assume an arbitrary complex R formula can be interpreted as a simple contrast.

## Run on a server

```bash
mkdir -p logs
nohup rnaseqflow run -c config/config.yaml --level intermediate \
  --cores 24 --jobs 6 --mem-mb 90000 > logs/rnaseqflow.log 2>&1 &
```

Change `--cores`, `--jobs`, and `--mem-mb` for your server; all three are consumed by Snakemake. Run from any directory: paths in the YAML are resolved relative to the project root. Avoid multiple simultaneous runs using the same project `results/` directory. Subsequent runs reuse complete Snakemake outputs and rerun modified dependencies. To avoid overwriting results when switching levels or experimental designs, use separate project directories or archive old results.

## Reproducibility and safety

- The CLI checks missing files, sample IDs, duplicate genes, integer counts, required metadata, per-group sample size, comparison names, and gene-set availability; R checks design-matrix rank and degrees of freedom.
- The default example uses 12 samples per condition and two balanced batches. For real differential analysis, biological replication and appropriate controls are essential; two per group is a minimum check, **not a design recommendation**.
- Gene selection: `padj < fdr` (BH) and `abs(log2FC) >= log2fc` is a *reporting filter after testing against zero*, not a formal equivalence/large-effect test.
- GO/KEGG/custom ORA uses a hypergeometric test and BH across all eligible pathways; GSEA ranks by the DESeq2 Wald statistic.
- Network correlations are exploratory. The default WGCNA soft threshold is user-provided and should be diagnosed for real data; module-trait tests use BH across reported tests. Biological causation is not implied.
- Default FASTQ Salmon index uses a transcriptome-only reference, not a decoy-aware index; false transcript assignment can occur. Decoy-aware generation and independent mapping QC are planned upgrades.
- This is bulk RNA-seq analysis, not single-cell RNA-seq, spatial transcriptomics, or a metatranscriptomic taxonomic pipeline.

See [Module design](docs/analysis_modules.md), [Statistical methods & limitations](docs/methods_and_limits.md), [Linux/HPC notes](docs/install_hpc.md), and [中文快速开始](docs/QUICKSTART_zh.md).

## Testing

```bash
python -m pip install -e '.[test]'
python -m pytest -q
rnaseqflow validate -c config/config.yaml
```

The CI workflow checks Python validation only. Before releasing a stable version, add external-tool integration tests and compare all R outputs against independently verified RNA-seq reference datasets.

## License

MIT. The bundled synthetic dataset is generated for demonstration only.

## Publish as a separate repository

If you own the GitHub account and have the [GitHub CLI](https://cli.github.com/) installed and authenticated, from the **extracted** project root:

```bash
gh auth login
bash scripts/publish_to_github.sh
```

This script creates a *new* `wangxiaopeng-git/RNAseqFlow` repo (public by default) with `gh repo create`, commits the extracted project and pushes `main`. It refuses to overwrite an existing remote repository. For a private repository set `GH_VISIBILITY=private`; set `GH_OWNER`/`GH_REPO` to customize. GitHub creation is a user-executed step, not an action that has already occurred.
