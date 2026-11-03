# Local validation of approach A

Validated in the prepared Linux cloud environment, using the source and wheel for RNAseqFlow 0.2.0. These are local results; GitHub Actions has not yet run on this branch remotely.

| Check | Observed result |
|---|---|
| Default Python suite | 43 passed; 2 opt-in integration tests skipped |
| Opt-in actual R/Snakemake suite | 2 passed, 0 skipped (97.40 s in final run) |
| Standalone R metadata/mapping/WGCNA checks | Passed; WGCNA executes without an external R profile and restores package search path |
| Fresh real pasilla counts analysis | 3/3 Snakemake steps completed; output verification passed |
| Independent direct DESeq2 comparison | All 8,148 genes and six original statistics matched; max absolute error 4.511094e-10 |
| Built wheel, installed outside source | Actual package import under independent site-packages; basic analysis 3/3 steps completed; output verification passed |
| Source vs bundled assets | Equality test passed |

The integration suite verifies actual numeric-coded sample categories, preservation of gene IDs with leading zeros, normal LFC shrinkage without changing original inference, WGCNA completion, exclusion of prior network artifacts when the next run skips WGCNA, and automatic regeneration after deleting an optional PCA plot. Python regressions also cover ambiguous maps, zero-valued YAML contrast labels, invalid thresholds/types, corrupt output hashes, malformed manifests, report path containment and URL encoding.

An independent whole-branch review identified stale conditional outputs and zero-valued YAML labels. Failing regressions reproduced both; fixes then passed the suite. The report root/encoding findings were also fixed with failing-to-passing regressions. R now builds each analysis in an empty private temporary directory and publishes only that run's inventory. Previous unreferenced files are preserved rather than deleting the user's output directory.

The real example used 7 samples, 14,599 input genes, and retained 8,148 genes. There were 1,108 at BH FDR <.05; additionally requiring |log2FC| >=1 produced 114 up and 110 down. The model adjusts for single-read/paired-end library type. These observations are not a biological gold standard. Both numerical paths use the same DESeq2 engine; they test orchestration and input alignment.

Scope intentionally not validated: real FASTQ→Salmon, real external Salmon imports, species-matched enrichment or real WGCNA biological findings. The seven-sample real example disables enrichment and network analysis. Full integration currently exercises synthetic counts and WGCNA, not every sequencing route.

Some R runs emit `stack imbalance` warnings from native dependency loading. They did not prevent successful execution or numerical agreement, but their package-level cause was not resolved here. Versions and logs are retained with the cloud validation artifacts; repeat in a fully pinned clean environment before production use. The environment specifications retain version ranges rather than representing a universal lockfile.
