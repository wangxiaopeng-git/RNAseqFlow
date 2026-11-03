# Execute with Rscript from repository root; real package calls, no R profile workaround.
source('scripts/analysis_helpers.R')
source('scripts/network.R')
`%||%` <- function(a,b) if(is.null(a)) b else a
fails <- function(expr) inherits(tryCatch(force(expr),error=identity),'error')
tmp <- tempfile(); dir.create(tmp)
sheet <- file.path(tmp,'samples.tsv')
writeLines(c('sample_id\tcondition\tage','001\t0\t2.5','002\t1\t3.5'),sheet)
s <- read_metadata(sheet,list(age='numeric'),list(list(factor='condition')))
stopifnot(identical(s$sample_id,c('001','002')),is.factor(s$condition),is.numeric(s$age))
stopifnot(fails(validate_txmap(data.frame(transcript_id=c('T1','T1'),gene_id=c('G1','G2')))))
stopifnot(fails(validate_txmap(data.frame(transcript_id='T1',gene_id=''))))
set.seed(123)
expr <- matrix(rnorm(200*24),200,24,dimnames=list(paste0('G',1:200),paste0('S',1:24)))
# Add coherent expression blocks to exercise actual module detection.
for(i in 1:4) expr[((i-1)*50+1):(i*50),] <- expr[((i-1)*50+1):(i*50),] + rep(rnorm(24)*2,each=50)
samples <- data.frame(sample_id=colnames(expr),condition=factor(rep(c('A','B'),12)))
rownames(samples) <- samples$sample_id
original_search <- search()
stopifnot(run_network(expr,samples,list(min_samples=20,min_module_size=20,traits='condition'),file.path(tmp,'network')))
stopifnot(identical(search(),original_search))
stopifnot(file.exists(file.path(tmp,'network/gene_modules.tsv')))
stopifnot(grepl('COMPLETED',readLines(file.path(tmp,'network/status.txt'))))
cat('R metadata, mapping and actual WGCNA checks passed\n')
