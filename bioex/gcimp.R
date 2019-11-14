# find list of probes that are associated with GCIMP genes Noushmehr2010
gcimp <- c('G0S2', 'RBP1', 'FABP5', 'CA3', 'RARRES2', 'OCIAD2', 'CBR1', 'PDPN',
  'LGALS3', 'CTHRC1', 'CCNA1', 'ARMC3', 'CHST6', 'C11orf63', 'GJB2', 'KIAA0746',
  'MOSC2', 'CHI3L1', 'RARRES1', 'AQP5', 'SPON2', 'RAB36', 'CHRDL2', 'TOM1L1',
  'BIRC3', 'LDHA', 'SEMA3E', 'FMOD', 'C10orf107', 'FLNC', 'TMEM22', 'TCTEX1D1',
  'DKFZP586H2123', 'TRIP4', 'SLC39A12', 'FLJ21963', 'CRYGD', 'LECT1', 'EPHX2',
  'LGALS8', 'C7orf46', 'F3', 'TTC12', 'ITGBL1', 'B3GNT5', 'NMNAT3', 'FZD6',
  'FKBP5', 'SLC25A20', 'MMP9')
assoc <- readr::read_tsv('~/data/ldgcData/OTHER/methylation-gene.tsv')

# use LDGC
wd <- getwd()
setwd('~/sync/phd/ldgc/ldgc-main')
source('main.R')
setwd(wd)

TCGA <- loadTCGA(dataPath, c('expression_rnaseq', 'methylation'))
TCGA <- list(expr=TCGA$expression_rnaseq, meth=TCGA$methylation)

gcimp_probes <- unique(subset(assoc,
  Name %in% colnames(TCGA$meth) & UCSC_RefGene_Name %in% gcimp)$Name)

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
pwgenes <- do.call(c, geneset)
geneset$gcimp <- gcimp

pw_probes <- unique(subset(assoc,
  Name %in% colnames(TCGA$meth) & UCSC_RefGene_Name %in% pwgenes)$Name)

examine_overlap <- function(set) {
  overlap <- matrix(NA, length(set), length(set))
  colnames(overlap) <- rownames(overlap) <- names(set)
  for (i in 1:length(set)) {
    for (j in 1:length(set)) {
      overlap[i, j] <- length(intersect(set[[i]], set[[j]]))
    }
  }
  print(overlap)
}
examine_overlap(geneset)

