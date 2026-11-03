# Signed WGCNA module detection. Input is a samples x genes normalized expression matrix.
run_network <- function(log_expression, samples, settings, outdir) {
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  min_samples <- as.integer(settings$min_samples %||% 20)
  if (ncol(log_expression) < min_samples) {
    writeLines(sprintf("SKIPPED: %d samples, minimum configured is %d", ncol(log_expression), min_samples),
               file.path(outdir,"status.txt"))
    return(invisible(FALSE))
  }
  genes <- order(apply(log_expression, 1, var), decreasing=TRUE)
  ntop <- min(length(genes), as.integer(settings$top_variable_genes %||% 5000))
  mat <- t(log_expression[genes[seq_len(ntop)], , drop=FALSE])
  good <- WGCNA::goodSamplesGenes(mat, verbose=0)
  if (!good$allOK) mat <- mat[good$goodSamples, good$goodGenes, drop=FALSE]
  min_module <- as.integer(settings$min_module_size %||% 20)
  if (ncol(mat) < max(2 * min_module, 30) || nrow(mat) < min_samples) {
    writeLines("SKIPPED: insufficient variable genes/samples after quality filtering", file.path(outdir,"status.txt"))
    return(invisible(FALSE))
  }
  modules <- WGCNA::blockwiseModules(mat, power=as.integer(settings$soft_power %||% 6),
                                    maxBlockSize=5000, networkType="signed", TOMType="signed",
                                    minModuleSize=min_module, mergeCutHeight=0.25,
                                    saveTOMs=FALSE, verbose=0)
  color <- modules$colors
  out <- data.frame(gene_id=colnames(mat), module=color)
  write.table(out,file=file.path(outdir,"gene_modules.tsv"),sep="\t", quote=FALSE,row.names=FALSE)
  eig <- WGCNA::moduleEigengenes(mat, colors=color)$eigengenes
  eig <- WGCNA::orderMEs(eig)
  write.table(data.frame(sample_id=rownames(eig), eig, check.names=FALSE),
              file=file.path(outdir,"module_eigengenes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  traits <- settings$traits
  phenotype <- samples[match(rownames(mat), samples$sample_id), ,drop=FALSE]
  traitdata <- list()
  for (name in traits) {
    x <- phenotype[[name]]
    if (is.numeric(x)) {
      if (anyNA(x)) stop(sprintf("Trait %s contains missing values", name))
      traitdata[[name]] <- x
    } else {
      x <- factor(x)
      if (nlevels(x) == 2) {
        traitdata[[paste0(name,"_",levels(x)[2],"_vs_",levels(x)[1])]] <- as.numeric(x == levels(x)[2])
      } else if (nlevels(x) > 2) {
        for (lev in levels(x)[-1]) traitdata[[paste0(name,"_",lev,"_vs_",levels(x)[1])]] <- as.numeric(x == lev)
      }
    }
  }
  if (length(traitdata)) {
    tm <- as.data.frame(traitdata, check.names=FALSE)
    rownames(tm) <- phenotype$sample_id
    cors <- cor(eig, tm, use="pairwise.complete.obs")
    pval <- WGCNA::corPvalueStudent(cors, nrow(eig))
    tab <- data.frame(module=rep(rownames(cors), times=ncol(cors)),
                      trait=rep(colnames(cors), each=nrow(cors)),
                      correlation=as.vector(cors), pvalue=as.vector(pval))
    tab$padj <- p.adjust(tab$pvalue, "BH")
    write.table(tab, file=file.path(outdir,"module_trait.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    if (nrow(cors)>1 && ncol(cors)>0) {
      svg(file.path(outdir,"module_trait.svg"),width=9,height=max(4, nrow(cors)*0.38+2))
      pheatmap::pheatmap(cors, cluster_rows=FALSE, cluster_cols=FALSE,
                        main="Module-trait Pearson correlations", display_numbers=TRUE,
                        number_format="%.2f")
      dev.off()
    }
  }
  writeLines("COMPLETED: signed WGCNA. Correlations are associations, not evidence of causality.",
             file.path(outdir,"status.txt"))
  invisible(TRUE)
}
