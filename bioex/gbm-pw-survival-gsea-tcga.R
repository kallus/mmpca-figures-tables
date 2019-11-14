

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

pdf("fig/compSurvival.pdf",height=16,width=16)
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
  plot_boxpoints(as.factor(patientCohorts),y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  #Sample type
  plot_boxpoints(as.factor(clinical$sampletype),y,ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  plot_scatter(clinical$mRNAstemness,y,xlab="mRNA stemness",ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  #Stemness
  plot_scatter(clinical$methStemness,y,xlab="meth stemness",ylab=paste("MVPC_",j,sep=""),mtext.cex = 0.8,pointColors=clinical$subtype)
  
  
  
  
}

#plotP <- rbind(outP1,outP2[1,,drop=F])
#colnames(plotP) <- 1:ncol(plotP)

#plotP[is.na(plotP)] <- 1
#barplot(-log10(plotP), ylab="-log10(p)",xlab="Component",beside=TRUE)
#abline(h=-log10(0.05),lwd=2,lty=2,col="darkolivegreen")

dev.off()

# --------


pdf('fig/Subtype_bars.pdf', width=0.8, height=2.9)
par(mar=c(0.25, 1.75, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
y <- -log10(pSubtype[nComp:1])
y[is.na(y)] <- 0
names <- rep('', 25)
barplot(y, horiz=TRUE,
        axes=FALSE, yaxs='i',xlim=c(0,40))
abline(v=-log10(0.05),col="tomato",lwd=1)
title(expression('E) Subtype'), adj=0, line=6)
title(ylab='Rank', line=0.75)
axis(2, c(30, seq(25.15, 1.25, length=5))-0.5, c(1, seq(5, 25, 5)), line=-1, tick=FALSE)
axis(3, c(0, 15, 30), c(0, 15, 30), line=0.25)
mtext(expression('-log'[10]*'(p)'), side=3, line=2.5, cex=0.5)
dev.off()


pdf('fig/Gcimp_bars.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.5, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
y <- -log10(pGCIMP[nComp:1])
y[is.na(y)] <- 0
barplot(y, horiz=TRUE,
        axes=FALSE, yaxs='i',xlim=c(0,25))
abline(v=-log10(0.05),col="tomato",lwd=1)
title(expression('F) G-CIMP'), adj=0, line=6)
axis(3, c(0, 10, 20), c(0, 10, 20), line=0.25)
mtext(expression('-log'[10]*'(p)'), side=3, line=2.5, cex=0.5)
dev.off()



pdf('fig/Sex_bars.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.5, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
y <- -log10(pSex[nComp:1])
y[is.na(y)] <- 0
barplot(y, horiz=TRUE,
        axes=FALSE, yaxs='i',xlim=c(0,80))
abline(v=-log10(0.05),col="tomato",lwd=1)
title(expression('G) Sex'), adj=0, line=6)
axis(3, c(0, 35, 70), c(0, 35, 70), line=0.25)
mtext(expression('-log'[10]*'(p)'), side=3, line=2.5, cex=0.5)
dev.off()


pdf('fig/Age_bars.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.5, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
y <- -log10(pAge[nComp:1])
y[is.na(y)] <- 0
barplot(y, horiz=TRUE,
        axes=FALSE, yaxs='i',xlim=c(0,15))
abline(v=-log10(0.05),col="tomato",lwd=1)
title(expression('H) Age'), adj=0, line=6)
axis(3, c(0, 5, 10), c(0, 5, 10), line=0.25)
mtext(expression('-log'[10]*'(p)'), side=3, line=2.5, cex=0.5)
dev.off()


pdf('fig/Survival_bars_adj.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.5, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
y <- -log10(pSurvival[nComp:1])
y[is.na(y)] <- 0
y2 <- -log10(pSurvivalAdjustedAge[1,nComp:1])
y2[is.na(y2)] <- 0
y3 <- -log10(pSurvivalAdjustedAgeSex[1,nComp:1])
y3[is.na(y3)] <- 0

barplot(rbind(y,y3), beside = T, horiz=TRUE,
        axes=FALSE, yaxs='i',xlim=c(0,4.5))
abline(v=-log10(0.05),col="tomato",lwd=1)
title(expression('I) Survival'), adj=0, line=6)
axis(3, c
     (0, 1.5, 3), c(0, 1.5, 3), line=0.25)
mtext(expression('-log'[10]*'(p)'), side=3, line=2.5, cex=0.5)
dev.off()



# -------
# Comparing survival with clustering:

#clusters <- readRDS('data/cluster.rds.gz')
#
#table(clusters[[1]]$r4_cl_id,clinical$gcimp)
#table(clusters[[1]]$r4_cl_id,clinical$subtype)
#
#
#pdf('fig/Survival_clusters.pdf', width=8, height=16)
#
#par(mfrow=c(2,1))
#plot_survival_clusters(clusters[[1]]$r3_cl_id,clinical,clusters$r3cols)
#plot_survival_clusters(clusters[[1]]$r4_cl_id,clinical,clusters$r4cols)
#
#dev.off()
#





# ------



pdf('fig/cluster_subtypes_bar2.pdf', width=1.5, height=1.5)
subtype_gcimp <- clinical$subtype
levels(subtype_gcimp) <- c('', 'CL', 'MS', 'NL', 'PN', 'PN-G')
subtype_gcimp[clinical$subtype  == 'PN' & clinical$gcimp == 'G-CIMP'] <- 'PN-G'
clusters <- readRDS('data/cluster.rds.gz')
par(mar=c(2, 3.5, 2, 0), cex=0.5)
cl <- clusters[[1]]$r3_cl_id
clusters[[1]]$r3_cl_id[cl == 2] <- 3
clusters[[1]]$r3_cl_id[cl == 3] <- 2
v <- table(clusters[[1]]$r3_cl_id,clinical$subtype)
barplot(prop.table(t(v),2),col=unlist(sapply(colnames(v),color_function)))
title(expression('E) Subtype/G-CIMP by cluster'), adj=0, line=1.5)
title(ylab='Subtype/G-CIMP proportion', line=2.25)
title(xlab='Clusters (rank 2)', line=1)
cols <- unlist(sapply(colnames(v),color_function))
v <- table(clusters[[1]]$r3_cl_id,subtype_gcimp)
barplot(prop.table(t(v),2),col=c(cols, 1), add=T, density=40)
dev.off()

clusters <- readRDS('data/cluster.rds.gz')
pdf('fig/cluster_subtypes_bar3.pdf', width=1.3, height=1.5)
par(mar=c(2, 1, 2, 0), cex=0.5)
cl <- clusters[[1]]$r4_cl_id
clusters[[1]]$r4_cl_id[cl == 2] <- 6
clusters[[1]]$r4_cl_id[cl == 3] <- 2
clusters[[1]]$r4_cl_id[cl == 4] <- 3
clusters[[1]]$r4_cl_id[cl == 5] <- 7
clusters[[1]]$r4_cl_id[cl == 6] <- 8
clusters[[1]]$r4_cl_id[cl == 7] <- 4
clusters[[1]]$r4_cl_id[cl == 8] <- 5
space <- rep(0.2, 8)
space[c(3, 5, 7)] <- 0.6
v <- table(clusters[[1]]$r4_cl_id,clinical$subtype)
barplot(prop.table(t(v),2),col=unlist(sapply(colnames(v),color_function)),
  space=space)
cols <- unlist(sapply(colnames(v),color_function))
v <- table(clusters[[1]]$r4_cl_id,subtype_gcimp)
barplot(prop.table(t(v),2),col=c(cols, 1), add=T, density=40, space=space)
title(xlab='Clusters (rank 3)', line=1)
dev.off()




# ------
enrichFun <- function(x,y){
  
  uX <- setdiff(unique(x),NA)
  uY <- setdiff(unique(y),NA)
  
  nUX <- length(uX)
  nUY <- length(uY)
   
  outSignif <- array(NA,c(nUX,nUY))
  rownames(outSignif) <- uX
  colnames(outSignif) <- uY
  
  for(i in 1:nUX){
    for(j in 1:nUY){
      outSignif[i,j] <- fisher.test(table(x==uX[i],y==uY[j]))$p.value
    }
  }
  
  outSignif
}

table1 <- enrichFun(clusters[[1]]$r3_cl_id,clinical$subtype)
table2 <- enrichFun(clusters[[1]]$r4_cl_id,clinical$subtype)


pdf('fig/cluster_subtypes.pdf', width=1.5, height=1.5)
par(mar=c(1, 1, 8, 0))
par(mfrow=c(1, 1), cex=0.5)
imm(table1)
axis(2, 1:nrow(table1), rownames(table1),line=-4, tick=FALSE, las=2)
axis(3, 1:ncol(table1), colnames(table1),line=-0.5, tick=FALSE, las=2)
title('Cluster overlaps', adj=0, line=6)
dev.off()

pdf('fig/cluster_subtypes2.pdf', width=1.5, height=1.5)
par(mar=c(1, 1, 8, 0))
par(mfrow=c(1, 1), cex=0.5)
imm(table2)
axis(2, 1:nrow(table2), rownames(table2),line=-5.5, tick=FALSE, las=2)
axis(3, 1:ncol(table2), colnames(table2),line=-0.5, tick=FALSE, las=2)
title('Cluster overlaps', adj=0, line=6)
dev.off()


# ------- 

#GSEA analysis
outputPathGSEA <- "data/gsea"
gseaMat <- Upat[,1:nComp]
colnames(gseaMat) <- 1:nComp

for(j in 1:ncol(gseaMat)){
  
  y <- gseaMat[, j]
  if(sd(y)==0){ next }
  
  isNA <- is.na(y)
  
  yNA <- is.na(y)	
  ySub <- y[!yNA]
  
  iPat <- intersect(names(ySub),rownames(expression))
  
  ySub <- ySub[match(iPat,names(ySub))]
  expressionSub <- expression[match(iPat,rownames(expression)),]
  
  myDF <- cbind(Genes = colnames(expressionSub), t(expressionSub))
  write.table(myDF,paste(outputPathGSEA,'/data/GSEA_',j,'.txt',sep=""),sep="\t",quote=F,row.names=F) # Needs one more column name
  
  dataSubSig <- as.character(ySub)
  
  printer = file(paste(outputPathGSEA,'/data/GSEA_',j,'.cls',sep=""),"w") 
  writeLines("#numeric",con=printer)
  writeLines(paste("#cluster_",j,sep=""),con=printer) 
  writeLines(paste(dataSubSig,sep="",collapse="\t"),con=printer) 
  close(printer) 
  
}

runV <- '3'
dir.create(paste(outputPathGSEA,'/output_',runV,sep=""), showWarnings = FALSE, recursive =T)
printer = file(paste(outputPathGSEA,'/GSEA_',runV,'.bat',sep=""),"w") 
wd <- ''
gSet <- 'c2.cp.reactome.v6.2.symbols.gmt'  #  'c2.cp.biocarta.v6.2.symbols.gmt' 'h.all.v6.2.symbols.gmt' #  #'msigdb.v6.2.symbols.gmt' #'c6.all.v6.2.symbols.gmt' #h.all.v6.2.symbols.gmt
for(j in 1:ncol(gseaMat)){
  y <- gseaMat[, j]
  if(sd(y)==0){ next }
  
  writeLines(paste('java -Xmx8g -cp ', wd, 'gsea-3.0.jar xtools.gsea.Gsea -res ', wd, 'data\\GSEA_',j,'.txt -cls "', wd, 'data\\GSEA_',j,'.cls#cluster_',j, ' -gmx .\\MSigDB\\', gSet,' -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label Clusters -metric Pearson -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ', wd, '.\\output_',runV,'\\output_',j, ' -gui false',sep=""),con=printer)
}
close(printer) 

exportGSEA("2",outputPathGSEA,gseaMat)
