#!/usr/bin/env Rscript
# Bulk RNA-seq: basic (QC), intermediate (DE/enrichment), advanced (WGCNA).
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript run_analysis.R CONFIG.yaml OUTPUT_DIR")
config_path <- normalizePath(args[1], mustWork=TRUE)
out_dir <- args[2]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
`%||%` <- function(a,b) if (is.null(a)) b else a
cfg <- yaml::read_yaml(config_path)
level <- cfg$analysis$level
mode <- cfg$project$input_mode
root <- normalizePath(getwd())
resolve_path <- function(x) {
  if (is.null(x) || !nzchar(x)) stop("Missing path")
  if (grepl("^/", x)) x else file.path(root, x)
}
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
if (!length(script_arg)) stop("Cannot determine the R script path")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE))
source(file.path(script_dir, "enrichment.R"))
source(file.path(script_dir, "network.R"))

samples <- read.delim(resolve_path(cfg$project$samples), check.names=FALSE,
                      stringsAsFactors=FALSE)
if (anyDuplicated(samples$sample_id)) stop("Duplicate sample IDs")
rownames(samples) <- samples$sample_id
samples$sample_id <- as.character(samples$sample_id)
for (name in setdiff(names(samples), c("sample_id","fastq_1","fastq_2","quant_dir"))) {
  if (is.character(samples[[name]])) samples[[name]] <- factor(samples[[name]])
}

if (mode == "counts") {
  tab <- read.delim(resolve_path(cfg$project$counts_matrix), check.names=FALSE,
                    stringsAsFactors=FALSE)
  if (anyDuplicated(tab$gene_id) || anyNA(tab$gene_id)) stop("Invalid gene IDs")
  count_data <- as.matrix(tab[, -1, drop=FALSE])
  rownames(count_data) <- tab$gene_id
  if (!setequal(colnames(count_data), samples$sample_id)) stop("Count/sample ID mismatch")
  count_data <- count_data[, samples$sample_id, drop=FALSE]
  storage.mode(count_data) <- "integer"
  if (anyNA(count_data) || any(count_data < 0)) stop("Counts must be non-negative integers")
  txi <- NULL
} else {
  txmap <- read.delim(resolve_path(cfg$project$tx2gene), check.names=FALSE, stringsAsFactors=FALSE)
  txmap <- unique(txmap[,c("transcript_id","gene_id")])
  files <- if (mode == "fastq") file.path(root, "results", "quant", samples$sample_id, "quant.sf")
           else vapply(samples$quant_dir, function(p) file.path(resolve_path(p), "quant.sf"), character(1))
  names(files) <- samples$sample_id
  if (!all(file.exists(files))) stop("One or more Salmon quant.sf files are missing")
  txi <- tximport::tximport(files, type="salmon", tx2gene=txmap, countsFromAbundance="no")
  count_data <- txi$counts
  write.table(data.frame(gene_id=rownames(txi$abundance), txi$abundance, check.names=FALSE),
              file=file.path(out_dir,"gene_tpm.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}

design <- if (level == "basic") ~ 1 else as.formula(cfg$analysis$design)
if (level != "basic") {
  allcols <- all.vars(design)
  if (anyNA(samples[,allcols,drop=FALSE])) stop("Missing design metadata")
  x <- model.matrix(design, data=samples)
  if (qr(x)$rank < ncol(x)) stop("Design matrix not full rank: check confounded factors/batches")
  if (nrow(x) <= ncol(x)) stop("Not enough residual degrees of freedom for design")
}
if (is.null(txi)) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=count_data, colData=samples, design=design)
} else {
  dds <- DESeq2::DESeqDataSetFromTximport(txi=txi, colData=samples, design=design)
}
min_count <- as.integer(cfg$analysis$min_count %||% 10)
min_samples <- as.integer(cfg$analysis$min_samples %||% 3)
keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
dds <- dds[keep,]
if (nrow(dds) < 2) stop("Too few genes after count filtering")
dds <- DESeq2::estimateSizeFactors(dds)
normal <- DESeq2::counts(dds, normalized=TRUE)
write.table(data.frame(gene_id=rownames(normal),normal,check.names=FALSE),
            file=file.path(out_dir,"normalized_counts.tsv"),sep="\t", quote=FALSE,row.names=FALSE)
sf <- DESeq2::sizeFactors(dds)
if (is.null(sf)) {
  # tximport's gene-specific length offsets can lead to normalization factors rather
  # than one scalar per sample; report that explicitly instead of inventing factors.
  sf_tab <- data.frame(sample_id=colnames(dds), size_factor=NA_real_,
                       normalization="gene_specific_normalization_factors")
} else {
  sf_tab <- data.frame(sample_id=colnames(dds), size_factor=as.numeric(sf),
                       normalization="sample_size_factor")
}
write.table(sf_tab,file=file.path(out_dir,"size_factors.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(samples, file=file.path(out_dir,"sample_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Basic QC uses log2(normalized counts + 1), not a claim of variance stabilization.
log_expr <- log2(normal + 1)
if (level != "basic") {
  dds <- DESeq2::DESeq(dds, quiet=TRUE)
  log_expr <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind=FALSE))
}
write.table(data.frame(gene_id=rownames(log_expr), log_expr,check.names=FALSE),
            file=file.path(out_dir,"transformed_expression.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

qcdir <- file.path(out_dir,"qc")
dir.create(qcdir, showWarnings=FALSE)
libsize <- colSums(DESeq2::counts(dds))
library_df <- data.frame(sample_id=names(libsize), mapped_gene_count=as.numeric(libsize))
write.table(library_df,file=file.path(qcdir,"library_sizes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

p <- ggplot2::ggplot(library_df, ggplot2::aes(x=sample_id, y=mapped_gene_count)) +
  ggplot2::geom_col(fill="#277e97") + ggplot2::coord_flip() +
  ggplot2::theme_bw() + ggplot2::labs(x="Sample", y="Gene-assigned fragment counts")
ggplot2::ggsave(file.path(qcdir,"library_sizes.svg"),p,width=8,height=max(4,nrow(library_df)*0.26))

if (ncol(log_expr) >= 3 && nrow(log_expr) >= 2) {
  pc <- prcomp(t(log_expr), center=TRUE, scale.=FALSE)
  pc_df <- data.frame(sample_id=rownames(pc$x), PC1=pc$x[,1], PC2=pc$x[,2], check.names=FALSE)
  pc_df <- merge(pc_df,samples,by="sample_id",sort=FALSE)
  write.table(pc_df,file=file.path(qcdir,"pca.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  colorvar <- if ("condition" %in% colnames(pc_df)) "condition" else NULL
  p <- ggplot2::ggplot(pc_df, ggplot2::aes(x=PC1,y=PC2)) + ggplot2::geom_point(size=3) +
    ggplot2::theme_bw() + ggplot2::labs(
       x=sprintf("PC1 (%.1f%%)",100*summary(pc)$importance[2,1]),
       y=sprintf("PC2 (%.1f%%)",100*summary(pc)$importance[2,2]),title="Expression PCA")
  if (!is.null(colorvar)) p <- p + ggplot2::aes(color=.data[[colorvar]])
  ggplot2::ggsave(file.path(qcdir,"pca.svg"),p,width=8,height=6)
}

cor_mat <- cor(log_expr, method="pearson")
write.table(data.frame(sample_id=rownames(cor_mat),cor_mat,check.names=FALSE),
            file=file.path(qcdir,"sample_correlation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if (ncol(cor_mat) >= 2) {
  grDevices::svg(file.path(qcdir,"sample_correlation.svg"),width=8,height=7)
  pheatmap::pheatmap(cor_mat, main="Pearson sample correlations", border_color=NA)
  grDevices::dev.off()
}

qc_files <- paste0("qc/", list.files(qcdir, full.names=FALSE))
sections <- list(list(title="Input and expression QC", detail="Check group balance, sequencing or quantification quality, PCA and library size outliers before interpreting DEGs.",
                      files=c("normalized_counts.tsv","transformed_expression.tsv","size_factors.tsv",
                              "sample_metadata.tsv",qc_files)))

if (level != "basic") {
  deg_files <- character()
  contrast_results <- list()
  for (ct in cfg$analysis$contrasts) {
    factor <- ct$factor
    if (!factor %in% names(samples) || !is.factor(samples[[factor]])) stop("Contrast factor must be categorical")
    if (!all(c(ct$numerator,ct$denominator) %in% levels(samples[[factor]]))) stop("Contrast levels not found")
    cdir <- file.path(out_dir,"contrasts",ct$name)
    dir.create(cdir,recursive=TRUE,showWarnings=FALSE)
    result <- DESeq2::results(dds, contrast=c(factor,ct$numerator,ct$denominator),
                             alpha=as.numeric(cfg$analysis$fdr))
    res <- as.data.frame(result)
    res$gene_id <- rownames(res)
    res <- res[,c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    alpha <- as.numeric(cfg$analysis$fdr)
    log2fc <- as.numeric(cfg$analysis$log2fc %||% 1)
    res$status <- "NS"
    res$status[!is.na(res$padj) & res$padj < alpha & res$log2FoldChange >= log2fc] <- "UP"
    res$status[!is.na(res$padj) & res$padj < alpha & res$log2FoldChange <= -log2fc] <- "DOWN"
    res <- res[order(res$padj, na.last=TRUE),]
    contrast_results[[ct$name]] <- res
    write.table(res, file=file.path(cdir,"deg_all.tsv"),sep="\t", quote=FALSE,row.names=FALSE)
    write.table(subset(res, status != "NS"), file=file.path(cdir,"deg_selected.tsv"),
                sep="\t",quote=FALSE,row.names=FALSE)
    g <- ggplot2::ggplot(res, ggplot2::aes(x=log2FoldChange,
                               y=-log10(pmax(padj, .Machine$double.xmin)),color=status)) +
      ggplot2::geom_point(alpha=0.65,size=1.3,na.rm=TRUE) +
      ggplot2::scale_color_manual(values=c("DOWN"="#3274b2","NS"="#999999","UP"="#bf4a47")) +
      ggplot2::theme_bw() + ggplot2::labs(y="-log10(BH-adjusted P)",title=ct$name)
    ggplot2::ggsave(file.path(cdir,"volcano.svg"),g,width=8,height=6)
    grDevices::svg(file.path(cdir,"ma.svg"),width=8,height=6)
    DESeq2::plotMA(result, ylim=c(-6,6), main=ct$name)
    grDevices::dev.off()
    top <- head(res$gene_id[res$status != "NS"], as.integer(cfg$analysis$top_heatmap_genes %||% 30))
    if (length(top) >= 2) {
      grDevices::svg(file.path(cdir,"top_deg_heatmap.svg"),width=8,height=max(5,length(top)*0.3))
      pheatmap::pheatmap(log_expr[top,,drop=FALSE], scale="row", main="Top selected DEGs (row Z-score)")
      grDevices::dev.off()
    }
    files <- paste0("contrasts/",ct$name,"/",c("deg_all.tsv","deg_selected.tsv","volcano.svg","ma.svg"))
    if (length(top) >= 2) files <- c(files,paste0("contrasts/",ct$name,"/top_deg_heatmap.svg"))
    if (isTRUE(cfg$analysis$enrichment$enabled)) {
      ep <- file.path(cdir,"enrichment")
      run_enrichment(res, rownames(dds), resolve_path(cfg$project$gene_sets_gmt),ep,
                     alpha,log2fc,
                     as.integer(cfg$analysis$enrichment$min_set_size %||% 5),
                     as.integer(cfg$analysis$enrichment$max_set_size %||% 500))
      en_files <- list.files(ep,full.names=FALSE)
      files <- c(files,paste0("contrasts/",ct$name,"/enrichment/",en_files))
    }
    sections[[length(sections)+1]] <- list(
      title=paste0("Differential expression: ",ct$name),
      detail=sprintf("DESeq2 Wald tests, BH FDR < %.3g; |log2FC| >= %.3g is a post-test reporting filter, not an effect-size hypothesis test.",alpha,log2fc),
      files=files)
  }
  if (level == "advanced" && length(contrast_results) >= 2) {
    names_ct <- names(contrast_results)
    comp_dir <- file.path(out_dir,"comparative")
    dir.create(comp_dir, showWarnings=FALSE, recursive=TRUE)
    jaccard <- matrix(0, nrow=length(names_ct), ncol=length(names_ct),
                      dimnames=list(names_ct,names_ct))
    for (i in seq_along(names_ct)) {
      hits_i <- contrast_results[[i]]$gene_id[contrast_results[[i]]$status != "NS"]
      for (j in seq_along(names_ct)) {
        hits_j <- contrast_results[[j]]$gene_id[contrast_results[[j]]$status != "NS"]
        denom <- length(union(hits_i,hits_j))
        jaccard[i,j] <- if (denom) length(intersect(hits_i,hits_j))/denom else NA_real_
      }
    }
    write.table(data.frame(contrast=rownames(jaccard), jaccard,check.names=FALSE),
                file=file.path(comp_dir,"deg_jaccard.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    if (all(is.finite(jaccard))) {
      svg(file.path(comp_dir,"deg_jaccard.svg"),width=7,height=6)
      pheatmap::pheatmap(jaccard,cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=TRUE,
                        main="DEG set Jaccard similarity")
      dev.off()
    }
    sections[[length(sections)+1]] <- list(title="Advanced cross-contrast DEG comparison",
      detail="Jaccard similarity of the selected DEG sets across multiple contrasts; an exploratory summary.",
      files=paste0("comparative/",list.files(comp_dir)))
  }
}

if (level == "advanced" && isTRUE(cfg$analysis$wgcna$enabled)) {
  net_dir <- file.path(out_dir,"network")
  status <- tryCatch(run_network(log_expr, samples, cfg$analysis$wgcna,net_dir),
                     error=function(e) {writeLines(paste("FAILED:",conditionMessage(e)),
                                                   file.path(net_dir,"status.txt")); stop(e)})
  network_files <- list.files(net_dir,full.names=FALSE)
  sections[[length(sections)+1]] <- list(title="Coexpression network (exploratory WGCNA)",
      detail="Signed coexpression modules and multiple-testing-corrected module-trait associations; independent biological validation is required.",
      files=paste0("network/",network_files))
}

# Put generated session and configuration records into the manifest.
sink(file.path(out_dir,"sessionInfo.txt"))
print(sessionInfo())
sink()
writeLines(yaml::as.yaml(cfg),file.path(out_dir,"effective_config.yaml"))
sections[[length(sections)+1]] <- list(title="Reproducibility", detail="Exact run settings and R/Bioconductor session versions.",
                                      files=c("sessionInfo.txt","effective_config.yaml"))
summary_obj <- list(project=cfg$project$name %||% "RNAseqFlow", level=level,
                    samples=nrow(samples), genes_after_filter=nrow(dds), sections=sections)
jsonlite::write_json(summary_obj, path=file.path(out_dir,"summary.json"), auto_unbox=TRUE, pretty=TRUE)
message(sprintf("RNAseqFlow completed: %d samples; %d genes retained; level=%s",nrow(samples),nrow(dds),level))
