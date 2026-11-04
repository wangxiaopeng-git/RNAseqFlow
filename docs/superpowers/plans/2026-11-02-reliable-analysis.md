# Reliable analysis A Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Deliver a reviewable PR improving RNAseqFlow reliability and scientific transparency.
**Architecture:** Preserve CLI/Snakemake/R. Validate before compute; record file-level provenance after analysis; verify outputs before reuse. Canonical scripts mirrored into wheel assets with equality tests.
**Tech Stack:** Python >=3.10, Snakemake >=8,<10, R/DESeq2/WGCNA.
**Spec:** docs/superpowers/specs/2026-11-02-reliable-analysis-design.md

## Global Constraints
- Default LFC shrinkage none; optional normal outputs separate estimates.
- Preserve raw Wald/BH inference and reporting filters.
- Existing isolated checkout, feature branch; no push/merge.
- No new alignment route or external services.

## Review Focus
- Ambiguous transcript IDs and empty TSV cells must fail before quantification import.
- Numeric-coded categories must remain factors, continuous covariates finite.
- Missing/corrupt optional results must invalidate cached analysis.
- Report paths must remain inside analysis root, with HTML escaped.
- Optional shrinkage must not silently change P values or selection.

### Task 1: Scientific validation and R execution
Files: configuration.py, scripts/run_analysis.R, scripts/network.R, environments, config examples, Python and R tests.
- [x] Add failing regressions for ambiguous mappings, covariate types, interaction rejection and invalid thresholds; run pytest.
- [x] Implement validation; R explicit metadata types, independent map checks, WGCNA lookup fix, normal shrinkage and explicit RNG seeds.
- [x] Run full Python suite and actual R synthetic integration (original WGCNA failure must disappear). Commit.

### Task 2: Integrity, provenance and report
Files: new provenance.py, cli.py, workflow/Snakefile, build_report.py, tests.
Interfaces: write_manifest(summary_path, config_path), verify_manifest(path) -> list[str]. Manifest contains SHA-256 output records, input records, configuration and versions.
- [x] Write failing missing/corrupt output and path safety/report tests; run them.
- [x] Track scripts and core outputs; write manifest after R; force analysis if output verification fails; CLI verify command.
- [x] Embed result plots and analysis metrics; include provenance links. Run full suite, deletion/corruption integration. Commit.

### Task 3: Reproducible example and delivery
Files: examples/pasilla, tests/integration, CI, README, Chinese docs, PR description.
- [x] Add executable fixed-source pasilla preparation and independent DESeq2 comparison, mirror consistency check.
- [x] Run real counts, normal shrinkage, numeric factor/continuous covariate, advanced synthetic and wheel installation checks.
- [x] Document validation evidence and limitations; fresh whole-branch review, fix important findings with regressions. Commit and produce portable patch/bundle.
