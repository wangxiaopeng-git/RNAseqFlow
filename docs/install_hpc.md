# Deployment on shared Linux servers

## Resource controls

`--cores` is the Snakemake thread pool, `--jobs` is passed as a resource named `concurrent_jobs`, used as `concurrent_jobs=1` on **all** rules, and `--mem-mb` is a global resource limiter. This permits explicit controls on actual simultaneous jobs rather than merely setting unused environment variables.

Jobs request fixed memory in the example workflow; for real 1–2 Gb vertebrate/crustacean transcriptomes, retune Salmon indexing and R analysis memory based on input dimensions. A Snakemake resource request is not an OS memory hard-limit.

```bash
mkdir -p logs
nohup rnaseqflow run -c config/config.yaml --cores 32 --jobs 4 --mem-mb 131072 \
  > logs/rnaseq.log 2>&1 &
```

Sample arrays, RAM needs, and job count must be chosen based on the available memory and CPU, with samples not oversubscribing per-rule threads. Never share one output directory across concurrent divergent configurations. To run two experimental designs at once, use two distinct project roots.

## Storage and privacy

Do not commit FASTQ/BAM/CRAM, unpublished genotype/phenotype data, credentials, or large Salmon indexes to GitHub. The `.gitignore` excludes common large inputs and outputs, but treat it as a safeguard rather than a security control. Sensitive project paths and data access policy are managed by the researcher.

## Reproducibility

Keep transcriptome FASTA, tx2gene mapping, GMT version, sample sheet, comparisons, and quantifier versions together. In v0.1 the R `sessionInfo()` and effective YAML are saved, but tool-level checksums and fully pinned software versions are not yet implemented. For public release, provide lock files or containers and record SHA-256 of input references and package versions.

## Workload manager

Snakemake supports HPC executors (e.g. Slurm), but this development preview only implements local-mode scheduling. Launching a single process within a scheduler reservation is appropriate when used within allocated CPU/memory resources. Avoid running heavy jobs directly on HPC login nodes.
