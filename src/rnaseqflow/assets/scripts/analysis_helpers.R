# Explicit metadata and mapping contracts shared by CLI analysis and R tests.
read_metadata <- function(path, types=list(), contrasts=list()) {
  samples <- read.delim(path, check.names=FALSE, colClasses='character', na.strings='')
  if (anyNA(samples$sample_id) || anyDuplicated(samples$sample_id)) stop('Invalid sample IDs')
  rownames(samples) <- samples$sample_id
  for (ct in contrasts) types[[ct$factor]] <- 'categorical'
  for (name in setdiff(names(samples), c('sample_id','fastq_1','fastq_2','quant_dir'))) {
    kind <- types[[name]]
    if (identical(kind, 'categorical')) {
      if (anyNA(samples[[name]])) stop(paste('Missing categorical covariate',name))
      samples[[name]] <- factor(samples[[name]])
    } else if (identical(kind, 'numeric')) {
      x <- suppressWarnings(as.numeric(samples[[name]]))
      if (any(!is.finite(x))) stop(paste('Invalid numeric covariate',name))
      samples[[name]] <- x
    } else {
      x <- type.convert(samples[[name]], as.is=TRUE)
      samples[[name]] <- if (is.character(x)) factor(x) else x
    }
  }
  samples
}

validate_txmap <- function(txmap) {
  txmap <- txmap[,c('transcript_id','gene_id'),drop=FALSE]
  if (!nrow(txmap) || anyNA(txmap) || any(!nzchar(trimws(as.matrix(txmap))))) stop('Empty transcript/gene IDs')
  if (any(as.matrix(txmap) != trimws(as.matrix(txmap)))) stop('Whitespace in transcript/gene IDs')
  txmap <- unique(txmap)
  if (anyDuplicated(txmap$transcript_id)) stop('Transcript maps to conflicting genes')
  txmap
}
