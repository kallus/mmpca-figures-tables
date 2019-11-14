

# -----
require("R.utils")

sourceDirectory('src/bioex/src')

source('src/bioex/utils.R')
 
res <- readRDS('src/bioex/gbm-mvpca-analysis-pw-patcov.rds.gz')
Vorder <- mvpca_order(res$inds, res$mvpca_res)


# use LDGC
wd <- getwd()
#ldgcPath <- '/home/pj/code/ldgc/ldgc-main' #
ldgcPath <- 'D:/Code/ldgc/ldgc-main'
ldgcPath <- '~/phd/ldgc/ldgc-main'
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

ldgcPath <- 'D:/Code/ldgc/'
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

nComp <- 6

#Survival analysis
pSurvival <- array(NA,nComp)
pSurvivalAdjustedAge <- array(NA,c(2,nComp))
pSurvivalAdjustedAgeSex <- array(NA,c(3,nComp))

pSubtype <- array(NA,nComp)
pGCIMP <- array(NA,nComp)

pAge <- array(NA,nComp)
pSex <- array(NA,nComp)

source('src/bioex/src/plot_survival.R')
for(j in 1:nComp){
pdf(paste("fig/compSurvival", j, ".pdf", sep=''),height=9,width=9)
par(mar=c(4.0, 4.5, 2, 1))
  
  y <- Upat[, j]
  z <- Uvar[, j]
  
  if(sd(y)==0){ next }
  
  
  layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,
                  5,5,5,6,6,6,7,7,7,8,8,8,8,
                  9,9,9,9,10,10,10,10), 
                nrow = 3, ncol = 12, byrow = TRUE))
  
  #Survival
  f <- plot_survival(y,clinical,z,j,methAssoc)
  pSurvival[j] <- f[[1]]
  pSurvivalAdjustedAge[,j] <- f[[2]]
  pSurvivalAdjustedAgeSex[,j] <- f[[3]]
  
  #
#  geneNames <- meth_genes(names(z),methAssoc)
#  isMeth <- names(z)%in%methAssoc$Name
#  
#  lapply(geneset,function(x){
#    ind <- sapply(geneNames,function(n){
#      any(str_split(n,';')%in%x)
#    })
#    
#    zPlot <- z[ind]
#    isMethPlot <- isMeth[ind]
#    plotOrder <- order(zPlot)
#    
#    plot(zPlot[plotOrder],ylim=range(z), pch=16, col = if(isMethPlot[plotOrder]) "purple" else "#444444")
#    
#  })
  
  
  #Subtype:
  pSubtype[j] <- kruskal.test(y~as.factor(clinical$subtype))$p.value
  
  xFactor <- factor(clinical$subtype,levels=c("CL","MS","NL","PN"))
  plot_boxpoints(xFactor,y,ylab=paste("MM-PC",j,sep=""),mtext.cex =
    0.8,pointColors=clinical$gcimp, xlab='Subtype')
  
  #Sample type
  tissuetype <- clinical$sampletype
  levels(tissuetype) <- c('NA', 'Primary', 'Recur.', 'Normal')
  plot_boxpoints(as.factor(tissuetype),y,ylab=paste("MM-PC",j,sep=""),mtext.cex
    = 0.8,pointColors=clinical$subtype, xlab='Tissue type')
  
  #G-CIMP
  xFactor <- factor(clinical$gcimp,levels=c("NON G-CIMP","G-CIMP"))
  pGCIMP[j] <- t.test(y~xFactor)$p.value
  plot_boxpoints(xFactor,y,ylab=paste("MM-PC",j,sep=""),mtext.cex =
    0.8,pointColors=clinical$subtype, xlab='G-CIMP status')
  
  #Sex
  xFactor <- factor(clinical$sex,levels=c("MALE","FEMALE"))
  pSex[j] <- t.test(y~xFactor)$p.value
  
  plot_boxpoints(xFactor,y,ylab=paste("MM-PC",j,sep=""),mtext.cex =
    0.8,pointColors=clinical$subtype, xlab='Sex')
  
  
  #Cohorts
  patFact <- as.factor(patientCohorts)
  levels(patFact) <- c('GE', 'Both', 'Meth.')
  plot_boxpoints(patFact,y,ylab=paste("MM-PC",j,sep=""),mtext.cex
  = 0.8,pointColors=clinical$subtype, xlab='Cohort')
  
  
  f <- plot_bar_signature(y,clinical,z,j,methAssoc, xlab=paste("MM-PC",j,sep=""))
 
  #Age:
  pAge[j] <- cor.test(clinical$age,y)$p.value
  plot_scatter(clinical$age,y,xlab="Age",ylab=paste("MM-PC",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  plot_scatter(clinical$mRNAstemness,y,xlab="mRNA stemness",ylab=paste("MM-PC",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  plot_scatter(clinical$methStemness,y,xlab="Meth. stemness",ylab=paste("MM-PC",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  
dev.off()
  
}

#plotP <- rbind(outP1,outP2[1,,drop=F])
#colnames(plotP) <- 1:ncol(plotP)

#plotP[is.na(plotP)] <- 1
#barplot(-log10(plotP), ylab="-log10(p)",xlab="Component",beside=TRUE)
#abline(h=-log10(0.05),lwd=2,lty=2,col="darkolivegreen")


# -------
# Comparing survival with clustering:

clusters <- readRDS('data/cluster.rds.gz')

table(clusters[[1]]$r4_cl_id,clinical$gcimp)
table(clusters[[1]]$r4_cl_id,clinical$subtype)


source('src/bioex/src/plot_survival_clusters.R')
pdf('fig/Survival_clusters.pdf', width=9, height=5)
par(mar=c(4.0, 4.5, 2, 1))
par(mfrow=c(1,2))
plot_survival_clusters(clusters[[1]]$r3_cl_id,clinical,clusters$r3cols)
plot_survival_clusters(clusters[[1]]$r4_cl_id,clinical,clusters$r4cols)
dev.off()






