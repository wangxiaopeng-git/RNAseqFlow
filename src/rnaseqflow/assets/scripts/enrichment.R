# GMT-based ORA and preranked GSEA, independent of species-specific online APIs.
run_enrichment <- function(res, universe, gmt_path, outdir, alpha, log2fc, min_size, max_size) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pathways <- fgsea::gmtPathways(gmt_path)
  pathways <- lapply(pathways, function(g) intersect(unique(g), universe))
  pathways <- pathways[vapply(pathways, length, 0L) >= min_size &
                         vapply(pathways, length, 0L) <= max_size]
  if (!length(pathways)) {
    writeLines("No gene sets met the overlap-size criteria.", file.path(outdir, "status.txt"))
    return(invisible(FALSE))
  }
  fdr <- !is.na(res$padj) & res$padj < alpha
  hits <- res$gene_id[fdr & !is.na(res$log2FoldChange) & abs(res$log2FoldChange) >= log2fc]
  hits <- intersect(hits, universe)
  N <- length(universe)
  ora <- do.call(rbind, lapply(names(pathways), function(name) {
    genes <- pathways[[name]]
    overlap <- intersect(genes, hits)
    m <- length(genes); n <- length(hits); k <- length(overlap)
    data.frame(pathway=name, universe_size=N, pathway_size=m, deg_size=n,
               overlap_count=k, pvalue=phyper(k-1, m, N-m, n, lower.tail=FALSE),
               overlap_genes=paste(overlap, collapse=";"))
  }))
  ora$padj <- p.adjust(ora$pvalue, method="BH")
  ora <- ora[order(ora$padj, ora$pvalue), , drop=FALSE]
  write.table(ora, file=file.path(outdir,"ora.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  finite <- !is.na(res$stat) & is.finite(res$stat) & !duplicated(res$gene_id)
  ranks <- setNames(as.numeric(res$stat[finite]), res$gene_id[finite])
  ranks <- sort(ranks, decreasing=TRUE)
  pathways_rank <- lapply(pathways, function(g) intersect(g, names(ranks)))
  pathways_rank <- pathways_rank[vapply(pathways_rank, length, 0L) >= min_size &
                                   vapply(pathways_rank, length, 0L) <= max_size]
  if (length(ranks) >= 15 && length(pathways_rank) > 0 && length(unique(ranks)) >= 5) {
    set.seed(20261102)
    gsea <- suppressWarnings(fgsea::fgseaMultilevel(pathways=pathways_rank, stats=ranks,
                                                    minSize=min_size, maxSize=max_size,
                                                    eps=0, nPermSimple=1000))
    if (nrow(gsea)) {
      gsea$leadingEdge <- vapply(gsea$leadingEdge, function(x) paste(x, collapse=";"), character(1))
      gsea <- as.data.frame(gsea)
      gsea <- gsea[order(gsea$padj, gsea$pval), , drop=FALSE]
      write.table(gsea, file=file.path(outdir,"gsea.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
  sig <- subset(ora, !is.na(padj) & padj < alpha & overlap_count > 0)
  if (nrow(sig)) {
    sig <- head(sig, 15)
    svg(file.path(outdir,"ora_top.svg"), width=9, height=6)
    par(mar=c(5, max(10, min(26, max(nchar(sig$pathway)) / 1.5)), 4, 2))
    barplot(rev(-log10(pmax(sig$padj, .Machine$double.xmin))), names.arg=rev(sig$pathway),
            horiz=TRUE, las=1, xlab="-log10(BH-adjusted P)", col="#237e91", main="ORA significant pathways")
    dev.off()
  }
  invisible(TRUE)
}
