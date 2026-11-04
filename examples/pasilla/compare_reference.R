# Independent DESeq2 reconstruction from the downloaded official metadata.
args <- commandArgs(TRUE)
if(length(args)!=1) stop('Usage: Rscript examples/pasilla/compare_reference.R PROJECT')
root <- normalizePath(args[1],mustWork=TRUE)
suppressPackageStartupMessages(library(DESeq2))
tab <- read.delim(file.path(root,'data/counts.tsv'),row.names=1,check.names=FALSE)
meta <- read.csv(file.path(root,'data/sample_annotation.csv'),check.names=FALSE)
rownames(meta) <- sub('fb$','',meta$file)
meta$condition <- factor(meta$condition)
meta$type <- factor(meta$type)
dds <- DESeqDataSetFromMatrix(as.matrix(tab[,rownames(meta)]),meta,~type+condition)
dds <- dds[rowSums(counts(dds)>=10)>=3,]
dds <- DESeq(dds,quiet=TRUE)
ref <- as.data.frame(results(dds,contrast=c('condition','treated','untreated'),alpha=.05))
actual <- read.delim(file.path(root,'results/analysis/contrasts/treated_vs_untreated/deg_all.tsv'),row.names='gene_id')
stopifnot(setequal(rownames(ref),rownames(actual)))
ref <- ref[rownames(actual),]
columns <- c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
errors <- sapply(columns,function(k) {
  stopifnot(identical(is.na(ref[[k]]),is.na(actual[[k]])))
  max(abs(ref[[k]]-actual[[k]]),na.rm=TRUE)
})
print(errors)
stopifnot(all(errors<1e-7))
write.table(data.frame(column=names(errors),max_absolute_error=errors),
            file.path(root,'reference-comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: all genes and six original DESeq2 statistics agree within 1e-7\n')
