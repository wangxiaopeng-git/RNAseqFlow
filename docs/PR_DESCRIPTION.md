# Harden RNAseqFlow analysis, provenance and real-data regression

RNAseqFlow's original counts workflow runs, but ambiguous transcript mappings can change gene assignments, WGCNA can resolve the wrong correlation function, and missing/old optional results can appear complete. This change preserves the existing Python/Snakemake/R architecture and strengthens the counts, Salmon and FASTQ interfaces without adding a new alignment route.

## Resulting behavior

- Reject conflicting/empty transcript mappings in Python and R. Preserve identifiers and numeric-coded categories; allow explicit finite numeric covariates. Reject unsupported interaction contrasts before computation.
- Make WGCNA's correlation dependency explicit and restore the prior R search environment on exit.
- Offer optional native DESeq2 normal LFC shrinkage in separate columns/MA plot; preserve original Wald/BH statistics and DEG selection.
- Build each R analysis in a fresh staging directory so old optional artifacts cannot enter a new report. Track core outputs and scripts in Snakemake; verify every advertised output before reuse and regenerate changed/missing artifacts.
- Record active-route input/output SHA-256, effective configuration, source hashes and tool/package versions. Add `rnaseqflow verify` and inline report plots/design/contrast counts.
- Include svglite in environment/doctor, bundled-asset consistency tests, Python and opt-in full R CI jobs, Chinese documentation, and a fixed-source 0.5 MB real pasilla example with independent DESeq2 numerical comparison.

## Validation

```bash
python -m pytest -q
Rscript tests/integration/r_checks.R
RNASEQFLOW_INTEGRATION=1 python -m pytest tests/integration -q
python examples/pasilla/prepare.py --dest /tmp/pasilla
rnaseqflow run -c /tmp/pasilla/config/config.yaml --cores 2 --jobs 1 --mem-mb 24000
rnaseqflow verify --results /tmp/pasilla/results/analysis
Rscript examples/pasilla/compare_reference.R /tmp/pasilla
```

Local validation results and exact scope are recorded in `docs/VALIDATION.md`. GitHub Actions is configured but its remote result must be checked after opening this PR.

## Compatibility and limits

Default shrinkage remains `none`; additive old configurations continue to work. Contrast factors are categorical; explicitly declare numeric-coded batches as categorical. Interaction formulas are intentionally rejected until a complete coefficient interface exists. The first reuse of pre-manifest results reruns analysis to establish integrity metadata. Input hashing adds sequential I/O, especially for FASTQ files. Existing unreferenced artifacts are retained on disk but do not represent current results; avoid concurrent runs sharing a project.

Real data validation covers counts/intermediate, not real FASTQ/Salmon, annotation-matched enrichment or biological network validity. WGCNA regression uses synthetic data. No splice analysis or STAR/featureCounts route is introduced. The seven-sample real example deliberately skips WGCNA and enrichment. Fixed manifests improve auditability but are not signed attestations or a fully locked environment.
