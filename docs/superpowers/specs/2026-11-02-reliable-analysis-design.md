# Reliable analysis A

User authorized direct implementation of approach A; defer STAR/featureCounts to B.
Preserve Python + Snakemake + R and existing counts/Salmon/FASTQ routes. Fix WGCNA correlation lookup, reject ambiguous tx2gene, detect missing/corrupt outputs and track scripts. Add explicit covariate types and optional native DESeq2 normal LFC shrinkage while retaining unshrunk inference. Reject interaction contrasts until explicitly supported. Add SHA-256 input/output provenance, software versions and report plots/metrics. Ship a fixed-source lightweight real pasilla example with independent DESeq2 regression and synthetic R integration tests.

Alternatives: extending existing architecture (selected) versus adding alignment routes (deferred B). Native normal shrinkage avoids adding an external R dependency; none remains default. Existing isolated cloud checkout is used on a feature branch; no extra worktree. Top-level scripts remain canonical; bundled copies must match in tests. No new external service or credentials. No fabricated enrichment; seven-sample pasilla skips WGCNA and enrichment.

Acceptance: Python regressions; actual synthetic advanced execution without R profile workarounds; real pasilla counts and numeric comparison; optional shrinkage execution; missing/corrupt output regeneration; safe report paths; wheel assets and entrypoint. Preserve existing source behavior unless documented; no automatic push or merge.
