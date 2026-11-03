# Real, lightweight pasilla example

This is real gene-count data from the Drosophila Pasilla RNAi study (Brooks et al., Genome Research 2011, PMID 20921232), distributed by Bioconductor pasilla 1.38.0. The count file is 498,373 bytes: 14,599 genes, three treated and four untreated samples. These are processed counts for the genes supplied by that package, not a fresh FASTQ alignment or an annotation-complete genome benchmark.

`prepare.py` downloads from a fixed upstream Git commit and checks hard-coded SHA-256 hashes, including the original sample annotation and DESCRIPTION/license metadata. It refuses to overwrite an existing project. No token is required; HTTPS access to raw.githubusercontent.com is required.

In the full environment, from the repository root:

```bash
python examples/pasilla/prepare.py --dest /tmp/pasilla
rnaseqflow validate -c /tmp/pasilla/config/config.yaml
rnaseqflow run -c /tmp/pasilla/config/config.yaml --cores 2 --jobs 1 --mem-mb 24000
rnaseqflow verify --results /tmp/pasilla/results/analysis
Rscript examples/pasilla/compare_reference.R /tmp/pasilla
```

The model is `~ library_type + condition`, adjusting for single-read versus paired-end libraries. Positive log2FC means treated / untreated. Count filtering requires at least three samples with count >=10. BH FDR is .05; |log2FC| >=1 is an additional reporting filter, not a test of a twofold minimum effect.

The observed development run retained 8,148 genes, with 1,108 at FDR <.05, and 114 up/110 down after the additional fold-change filter. Package versions can affect these numbers: the independent comparison is computed again rather than treating these numbers as universal ground truth. It checks all six original DESeq2 statistics at absolute tolerance 1e-7, using matching sample order and factor references. Both paths use DESeq2; this checks orchestration, not independent biological validity.

Open `/tmp/pasilla/results/analysis/report/index.html` locally for PCA, volcano, MA and heatmap plots. PCA uses VST without removing library type; the heatmap uses genes selected by differential testing and is not independent confirmation. Seven samples are insufficient for the workflow's default WGCNA requirement. Enrichment is disabled because this example does not supply a verified FlyBase-matched GMT. Gene-level differential expression does not validate Pasilla's alternative-splicing effects.

Try `analysis.lfc_shrinkage: normal` to add normal-prior shrinkage columns and `ma_shrunk.svg`; original Wald P values, BH adjustment and selection remain unchanged. Use separate project directories when comparing designs.
