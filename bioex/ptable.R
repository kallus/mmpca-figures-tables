

# -----
require("R.utils")

sourceDirectory('src/bioex/src')

source('src/bioex/utils.R')
 
res <- readRDS('src/bioex/gbm-mvpca-analysis-pw-patcov.rds.gz')
Vorder <- mvpca_order(res$inds, res$mvpca_res)


# use LDGC
wd <- getwd()
#ldgcPath <- '/home/pj/code/ldgc/ldgc-main' #
#ldgcPath <- 'D:/Code/ldgc/ldgc-main'
ldgcPath <- '~/sync/phd/ldgc/ldgc-main'
setwd(ldgcPath)
source('main.R')
TCGA <- loadTCGA(dataPath, c('clinical', 'mutation','expression_rnaseq'))
setwd(wd)


sol <- res$mvpca_res$solution

V <- sol$V
D <- sol$D
U <- lapply(1:5, function(i) V[[i]] %*% diag(D[, i]))
patients <- c(rownames(res$x[[1]]),rownames(res$x[[2]]), rownames(res$x[[4]]))

patientCohorts <- c(array(1,nrow(res$x[[1]])),array(2,nrow(res$x[[2]])),array(3,nrow(res$x[[4]])))

Upat <- do.call(rbind, U[1:3])
Uvar <- do.call(rbind, U[4:5])

rownames(Upat) <- patients
rownames(Uvar) <- c(colnames(res$x[[1]]), colnames(res$x[[3]]))

clinical <- TCGA$clinical[patients, ]
expression <- TCGA$expression_rnaseq


stemness <- readr::read_tsv('./data/stemnessCombined.txt')
stemness$TrimmedID <- gsub('.{1}$', '', str_replace_all(stemness$TrimmedID,'-','.'))

rownames(clinical)%in%stemness$TrimmedID
stopifnot(all(stemness$TrimmedID %in% rownames(clinical)))

clinical$mRNAstemness <- array(NA,nrow(clinical))
clinical$mRNAstemness[match(stemness$TrimmedID,rownames(clinical))] <- stemness$mRNAsi

clinical$methStemness <- array(NA,nrow(clinical))
clinical$methStemness[match(stemness$TrimmedID,rownames(clinical))] <- stemness$mDNAsi


# ---------

#ldgcPath <- 'D:/Code/ldgc/'
#ldgcPath <- '/home/pj/code/ldgc/'
ldgcPath <- '~/data/'
methAssoc <- readr::read_tsv(paste(ldgcPath, 'ldgcData/OTHER/methylation-gene.tsv',sep=""))

# ----------


# gene sets
geneset <- list() # ordered by length:
geneset$mtor <- read.csv('../../data/KEGG_pathways/map04150_genelist',stringsAsFactors=FALSE)[, 1]
geneset$ras <- read.csv('../../data/KEGG_pathways/map04014_genelist',stringsAsFactors=FALSE)[, 1]
geneset$mapk <- read.csv('../../data/KEGG_pathways/map04010_genelist',stringsAsFactors=FALSE)[, 1]
geneset$pi3kakt <- read.csv('../../data/KEGG_pathways/map04151_genelist',stringsAsFactors=FALSE)[, 1]


# ---------

require("rms")

nComp <- 25

#Survival analysis
pSurvival <- array(NA,nComp)
pSurvivalAdjustedAge <- array(NA,c(2,nComp))
pSurvivalAdjustedAgeSex <- array(NA,c(3,nComp))

pSubtype <- array(NA,nComp)
pGCIMP <- array(NA,nComp)

pAge <- array(NA,nComp)
pSex <- array(NA,nComp)
pCohort <- array(NA,nComp)
pTissue <- array(NA,nComp)
pRnaStemness <- array(NA,nComp)
pMethStemness <- array(NA,nComp)

for(j in 1:nComp){
  
  y <- Upat[, j]
  z <- Uvar[, j]
  
  if(sd(y)==0){ next }
  
  
  layout(matrix(c(1,1,1,1,2,2,3,3,3,4,4,4,
                  5,5,5,6,6,6,7,7,7,8,8,8,
                  9,9,9,10,10,10,11,11,11,12,12,12,
                  13,13,13,14,14,14,15,15,15,16,16,16), 
                nrow = 4, ncol = 12, byrow = TRUE))
  
  #Survival
  f <- plot_survival(y,clinical,z,j,methAssoc)
  pSurvival[j] <- f[[1]]
  pSurvivalAdjustedAge[,j] <- f[[2]]
  pSurvivalAdjustedAgeSex[,j] <- f[[3]]
  #Age:
  pAge[j] <- cor.test(clinical$age,y)$p.value
  plot_scatter(clinical$age,y,xlab="Age",ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Empty
  plot.new()
  
  #
  geneNames <- meth_genes(names(z),methAssoc)
  isMeth <- names(z)%in%methAssoc$Name
  
  lapply(geneset,function(x){
    ind <- sapply(geneNames,function(n){
      any(str_split(n,';')%in%x)
    })
    
    zPlot <- z[ind]
    isMethPlot <- isMeth[ind]
    plotOrder <- order(zPlot)
    
    plot(zPlot[plotOrder],ylim=range(z), pch=16, col = if(isMethPlot[plotOrder]) "purple" else "#444444")
    
  })
  
  
  #Subtype:
  pSubtype[j] <- kruskal.test(y~as.factor(clinical$subtype))$p.value
  
  xFactor <- factor(clinical$subtype,levels=c("CL","MS","NL","PN"))
  plot_boxpoints(xFactor,y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$gcimp)
  
  #G-CIMP
  xFactor <- factor(clinical$gcimp,levels=c("NON G-CIMP","G-CIMP"))
  pGCIMP[j] <- t.test(y~xFactor)$p.value
  plot_boxpoints(xFactor,y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Sex
  xFactor <- factor(clinical$sex,levels=c("MALE","FEMALE"))
  pSex[j] <- t.test(y~xFactor)$p.value
  
  plot_boxpoints(xFactor,y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  #Cohorts
  pCohort[j] <- plot_boxpoints(as.factor(patientCohorts),y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  #Sample type
  pTissue[j] <- plot_boxpoints(as.factor(clinical$sampletype),y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  pRnaStemness[j] <- plot_scatter(clinical$mRNAstemness,y,xlab="mRNA stemness",ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  pMethStemness[j] <- plot_scatter(clinical$methStemness,y,xlab="meth stemness",ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  
  
}

tbl <- data.frame(survival=pSurvival, surv_adj=pSurvivalAdjustedAge[1, ],
  age=pAge, subtype=pSubtype, gcimp=pGCIMP, sex=pSex, cohort=pCohort,
  tissue=pTissue, mRNA_stemness=pRnaStemness, meth_stemness=pMethStemness)

load('data/gsea/pathMat_1.RData')

tbl <- cbind(tbl[1:6, ], pathMat[1:6, c(39, 2, 50, 46, 11, 40)])

dot <- function(c) {
  if (c == '.') {
    return('$\\cdot$')
  }
  if (c == ' ') {
    return('')
  }
  return(c)
}
for (i in 1:ncol(tbl)) {
  cat(names(tbl)[i])
  cat(' & ')
  cat(dot(gtools::stars.pval(tbl[1, i])))
  for (j in 2:6) {
    cat(' & ')
    cat(dot(gtools::stars.pval(tbl[j, i])))
  }
  cat(' \\\\\n')
}

