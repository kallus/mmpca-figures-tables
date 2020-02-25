# utility methods
row_apply <- function(x, f) {
  res <- apply(x, 1, f)
  if (is.matrix(res)) {
    return(t(res))
  }
  return(res)
}
col_apply <- function(x, f) apply(x, 2, f)

top_k_mad_columns <- function(x, k) {
  x[, rank(-col_apply(x, function(x) mad(x, na.rm=TRUE)),
    ties.method='first') <= k]
}

with_columns <- function(x, names) {
  y <- matrix(NA, nrow(x), length(names), dimnames=list(rownames(x), names))
  is <- intersect(names, colnames(x))
  y[, is] <- x[, is]
  return(y)
}

with_rows <- function(x, names) {
  y <- matrix(NA, length(names), ncol(x), dimnames=list(names, colnames(x)))
  is <- intersect(names, rownames(x))
  y[is, ] <- x[is, ]
  return(y)
}

with_rows_columns <- function(x, rnames, cnames) {
  with_rows(with_columns(x, cnames), rnames)
}

# methods for normalization
quantnormalize <- function(x.na) {
  x <- x.na[!is.na(x.na)]
  x <- rank(x, ties.method='average')
  x <- qnorm(x/(length(x)+1))
  x <- (x-mean(x))/sd(x)
  x.na[!is.na(x.na)] <- x
  return(x.na)
}

center_scale_columns <- function(x) scale(x, center=TRUE, scale=TRUE)
center_columns <- function(x) scale(x, center=TRUE, scale=FALSE)
center_rows <- function(x) t(scale(t(x), center=TRUE, scale=FALSE))
quant_normalize_rows <- function(x) row_apply(x, quantnormalize)

if (!file.exists('biodata.rds.gz')) {
# use LDGC
wd <- getwd()
setwd('~/sync/phd/ldgc/ldgc-main')
source('main.R')
setwd(wd)

TCGA <- loadTCGA(dataPath, c('expression_rnaseq', 'methylation'))
TCGA <- list(expr=TCGA$expression_rnaseq, meth=TCGA$methylation)

TCGA$expr <- quant_normalize_rows(TCGA$expr)
TCGA$meth <- center_rows(TCGA$meth)

# remove 0 MAD columns
TCGA$expr <- TCGA$expr[,
  which(col_apply(TCGA$expr, function(x) mad(x, na.rm=TRUE)) > 0)]
TCGA$meth <- TCGA$meth[,
  which(col_apply(TCGA$meth, function(x) mad(x, na.rm=TRUE)) > 0)]

# remove the half of genes with lowest variance
qv <- quantile(col_apply(TCGA$expr, function(x) var(x, na.rm=TRUE)), 0.5)
TCGA$expr <- TCGA$expr[,
  which(col_apply(TCGA$expr, function(x) var(x, na.rm=TRUE)) > qv)]

# remove the half of probes with lowest variance
qv <- quantile(col_apply(TCGA$meth, function(x) var(x, na.rm=TRUE)), 0.5)
TCGA$meth <- TCGA$meth[,
  which(col_apply(TCGA$meth, function(x) var(x, na.rm=TRUE)) > qv)]

# normalize
TCGA$expr <- center_scale_columns(TCGA$expr)
TCGA$meth <- center_columns(TCGA$meth)

cohorts <- list(
  TCGA_expr=setdiff(rownames(TCGA$expr), rownames(TCGA$meth)),
  TCGA_both=intersect(rownames(TCGA$meth), rownames(TCGA$expr)),
  TCGA_meth=setdiff(rownames(TCGA$meth), rownames(TCGA$expr)))

patcov <- list(
  cov(t(TCGA$expr[cohorts$TCGA_expr, ]), t(TCGA$expr[cohorts$TCGA_both, ]),
    use='na.or.complete'),
  cov(t(TCGA$meth[cohorts$TCGA_both, ]), t(TCGA$meth[cohorts$TCGA_meth, ]),
    use='na.or.complete')
  )

# gene sets
geneset <- list() # ordered by length:
geneset$mtor <- read.csv('~/data/KEGG_pathways/map04150_genelist',
  stringsAsFactors=FALSE)[, 1]
geneset$ras <- read.csv('~/data/KEGG_pathways/map04014_genelist',
  stringsAsFactors=FALSE)[, 1]
geneset$mapk <- read.csv('~/data/KEGG_pathways/map04010_genelist',
  stringsAsFactors=FALSE)[, 1]
geneset$pi3kakt <- read.csv('~/data/KEGG_pathways/map04151_genelist',
  stringsAsFactors=FALSE)[, 1]

pw <- unique(do.call(c, geneset[1:4]))

data_genes <- unique(colnames(TCGA$expr))
pw <- intersect(pw, data_genes)

data_probes <- unique(colnames(TCGA$meth))

assoc <- readr::read_tsv('~/data/ldgcData/OTHER/methylation-gene.tsv')
pw_probes <- unique(subset(assoc, UCSC_RefGene_Name %in% pw)$Name)
pw_probes <- intersect(pw_probes, data_probes)

types <- list(pw_genes=pw, pw_probes=pw_probes)


pw_assoc <- loadMethylationGeneAssociations(types$pw_probes, types$pw_genes)
pw_assoc <- as.matrix(Reduce('+', pw_assoc))
pw_assoc[pw_assoc > 0] <- -1
pw_assoc <- center_columns(center_rows(pw_assoc))

# construct input to MVPCA
x <- list(
  with_rows_columns(TCGA$expr, cohorts$TCGA_expr, types$pw_genes),
  with_rows_columns(TCGA$expr, cohorts$TCGA_both, types$pw_genes),
  with_rows_columns(TCGA$meth, cohorts$TCGA_both, types$pw_probes),
  with_rows_columns(TCGA$meth, cohorts$TCGA_meth, types$pw_probes),
  patcov[[1]], patcov[[2]])
x <- lapply(x, function(data) {
  # scale according to first PC and number of observations
  data0 <- data
  data0[is.na(data0)] <- 0
  data * sqrt(nrow(data)) / svd(data0)$d[1]
})
inds <- matrix(NA, 6, 2)
inds[1, ] <- c(1, 4)
inds[2, ] <- c(2, 4)
inds[3, ] <- c(2, 5)
inds[4, ] <- c(3, 5)
inds[5, ] <- c(1, 2)
inds[6, ] <- c(2, 3)

#pdf('pw.pdf', width=19, height=20)
#NAmat <- function(rmat, cmat) matrix(NA, nrow(rmat), ncol(cmat))
#alldata <- cbind(
#  rbind(x[[1]], x[[2]], NAmat(x[[4]], x[[1]])),
#  rbind(NAmat(x[[1]], x[[3]]), x[[3]], x[[4]]))
#imm(alldata, asp=TRUE)
#abline(h=cumsum(sapply(x[1:2], nrow)),
#  v=cumsum(sapply(x[c(2, 4, 5)], ncol)))
#dev.off()
data <- list(x=x, inds=inds)
saveRDS(data, 'biodata.rds.gz')
}

data <- readRDS('biodata.rds.gz')

set.seed(1)
nparam <- 3
lambda <- matrix(rep(exp(seq(-8, 0, length=10)), nparam), ncol=nparam)
results <- mmpca::mmpca(data$x, data$inds, k=40, trace=1, lambda=lambda, 
  init_theta=TRUE, cachepath='cache', enable_sparsity=TRUE,
  enable_variable_selection=FALSE)
fh <- gzfile('gbm-mvpca-analysis-pw-patcov.rds.gz')
saveRDS(list(x=x, inds=inds, mvpca_res=results), fh)
close(fh)
