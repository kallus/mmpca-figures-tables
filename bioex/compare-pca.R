source('src/bioex/utils.R')
res <- readRDS('src/bioex/gbm-mvpca-analysis-pw-patcov.rds.gz')
V <- res$mvpca_res$solution$V
x <- res$x
rownames(V[[1]]) <- rownames(x[[1]])
rownames(V[[2]]) <- rownames(x[[2]])
rownames(V[[3]]) <- rownames(x[[4]])
rownames(V[[4]]) <- colnames(x[[1]])
rownames(V[[5]]) <- colnames(x[[3]])
for (i in 1:5) {
  colnames(V[[i]]) <- paste('MVPC', 1:ncol(V[[i]]), sep='')
}
res$mvpca_res$solution$V <- V

# pdf('data.pdf', height=11.7, width=16.5)
# par(mar=c(0, 0, 0, 0))
# mvpca_data_plot(res$x, res$inds, res$mvpca_res, Vorder)
# dev.off()
#
# pdf('model.pdf', height=11.7, width=16.5)
# par(mar=c(0, 0, 0, 0))
# mvpca_model_plot(res$inds, res$mvpca_res, Vorder, 1:11)
# dev.off()
#
# pdf('model_all.pdf', height=7.5, width=16.5)
# par(mar=c(0, 0, 0, 0))
# mvpca_model_plot(res$inds, res$mvpca_res, Vorder, c(1, 2, 5, 11))
# dev.off()
#
# pdf('model_ge.pdf', height=7.5, width=16.5)
# par(mar=c(0, 0, 0, 0))
# mvpca_model_plot(res$inds, res$mvpca_res, Vorder, c(3, 7, 8))
# dev.off()
#
# pdf('model_meth.pdf', height=7.5, width=16.5)
# par(mar=c(0, 0, 0, 0))
# mvpca_model_plot(res$inds, res$mvpca_res, Vorder, c(4, 6, 9, 10))
# dev.off()

mvcl <- mvpcaclust(res$x, res$inds, res$mvpca_res, Vorder)

for (i in 1:5) {
  fname <- paste('tree', i, '.pdf', sep='')
#  pdf(fname, height=2, width=11.7*length(Vorder[[i]])/384)
#  par(mar=c(1, 0, 3, 0))
#  plot(mvcl[[i]], main='', axes=F, ylab='', xlab='', sub='', hang=-1, cex=0.2,
#    labels=FALSE)
#  dev.off()
}

#pdf('fig/D.pdf', width=0.75, height=2.9)
D <- res$mvpca_res$solution$D
D[D == 0] <- NA
par(mar=c(0.25, 2, 7, 0.25))
par(mfrow=c(1, 1), cex=0.5)
imm(D[1:25, ])
axis(2, c(1, seq(5, 25, 5)), line=-1.35, tick=FALSE)
axis(3, 1:5, c('Cohort 1', 'Cohort 2', 'Cohort 3', 'Genes', 'Meth. sites'),
  line=-0.75, tick=FALSE, las=2)
title(expression('A) View'), adj=0, line=6)
title(ylab='Rank', line=0.5)
#dev.off()

# for each block: which factors are important
R2 <- res$mvpca_res$solution$R2_blockwise
for (i in 1:6) R2[, i] <- R2[, i] / sum(R2[, i])
R2[R2 == 0] <- NA
#pdf('fig/R2.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.25, 7, 0.25))
par(mfrow=c(1, 1), cex=0.5)
imm(R2[1:25, ])
axis(3, 1:6, c(
    expression('GE in'~C[1]),
    expression('GE in'~C[2]),
    expression('Meth. in'~C[2]),
    expression('Meth. in'~C[3]),
    expression('Sim.'~C[1]*C[2]),
    expression('Sim.'~C[2]*C[3])),
  line=-0.75, tick=FALSE, las=2)
title(expression('B) Block'), adj=0, line=6)
#dev.off()

row_apply <- function(x, f) {
  res <- apply(x, 1, f)
  if (is.matrix(res)) {
    return(t(res))
  }
  return(res)
}

# color signature for structure
Dbin <- ifelse(is.na(D), 0, 1)
Dstruct <- as.numeric(as.factor(row_apply(Dbin, digest::digest)))[1:24]
Dstruct <- c(Dstruct, NA)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#pdf('fig/Dstruct.pdf', width=0.15, height=2.9)
par(mar=c(0.25, 0.25, 7, 0.25))
par(mfrow=c(1, 1), cex=0.5)
imm(matrix(Dstruct, length(Dstruct), 1),
  colors=gg_color_hue(max(Dstruct, na.rm=TRUE)))
axis(3, 1, 'Structure', line=-0.75, tick=FALSE, las=2)
title(expression('C)'), adj=0, line=6)
#dev.off()

#pdf('fig/Total.pdf', width=0.65, height=2.9)
par(mar=c(0.25, 0.5, 7, 0))
par(mfrow=c(1, 1), cex=0.5)
R2_total <- res$mvpca_res$solution$R2_total
barplot(R2_total[25:1]/sum(R2_total), horiz=TRUE, xlim=c(0, 0.27),
  axes=FALSE, yaxs='i')
title(expression('D) Total'), adj=0, line=6)
axis(3, c(0, 0.1, 0.2), c(0, 10, 20), line=0.25)
mtext(expression(R^2~'(% of tot.)'), side=3, line=2.5, cex=0.5)
#dev.off()


# use LDGC
wd <- getwd()
setwd('~/sync/phd/ldgc/ldgc-main') # setwd('~/Code/ldgc-main')
source('main.R')
TCGA <- loadTCGA(dataPath, c('clinical', 'mutation'))
setwd(wd)


sol <- res$mvpca_res$solution

with_rows <- function(x, names) {
  y <- matrix(NA, length(names), ncol(x), dimnames=list(names, colnames(x)))
  is <- intersect(names, rownames(x))
  y[is, ] <- x[is, ]
  return(y)
}

agecl <- list(cutree(mvcl[[2]], k=10) == 7, cutree(mvcl[[3]], k=16) == 7)

V <- sol$V
D <- sol$D
U <- lapply(1:5, function(i) V[[i]] %*% diag(D[, i]))
Upat <- do.call(rbind, U[2:3])
patients <- c(rownames(res$x[[2]]), rownames(res$x[[4]]))
clinical <- TCGA$clinical[patients, ]
clinical$cl <- do.call(c, agecl)
clinical$MVPC2 <- Upat[, 2]
clinical$IDH <- ifelse(with_rows(TCGA$mutation, patients)[, 'IDH1'], 'Mutated',
  'Wild type')
library(ggplot2)
#pdf('fig/age.pdf', width=4, height=4)
plotdf <- subset(clinical, subtype == 'PN')# | subtype == 'NL')
fit <- lm(age~MVPC2, plotdf)
paste("Adj R2 =",signif(summary(fit)$adj.r.squared, 2))
ggplot(plotdf, aes(x=MVPC2, y=age)) +
  geom_point() + stat_smooth(method = "lm", col = "red") + theme_bw() +
  ylab('Patient age (years)') +  ggtitle('B) Patient subset: proneural subtype')
#dev.off()

patients <- c(rownames(res$x[[1]]), rownames(res$x[[2]]), rownames(res$x[[4]]))
clinical <- TCGA$clinical[patients, ]
clinical$IDH <- ifelse(with_rows(TCGA$mutation, patients)[, 'IDH1'], 'Mutated',
  'Wild type')
clinical$IDH[is.na(clinical$IDH)] <- 'Unknown'
levels(clinical$subtype)[1] <- 'Unknown'
Upat <- do.call(rbind, U[1:3])
for (i in 1:24) {
  pc <- paste('MVPC', i, sep='')
  clinical[[pc]] <- Upat[, i]
}

#pdf('fig/subtype.pdf', width=5, height=4)
ggplot(clinical, aes(x=MVPC1, y=MVPC2, color=subtype, shape=IDH)) +
  geom_point() + theme_bw() + ggtitle('A) Patients') +
  labs(color='Subtype', shape='IDH status')
#dev.off()

sge <- svd(rbind(x[[1]], x[[2]]))
smeth <- svd(rbind(ifelse(is.na(x[[3]]), 0, x[[3]]),
  ifelse(is.na(x[[4]]), 0, x[[4]])))

jivesol <- r.jive::jive(
  list(t(res$x[[2]]), t(ifelse(is.na(res$x[[3]]), 0, res$x[[3]]))),
  method='perm', center=FALSE, scale=FALSE, orthIndiv=FALSE)
jiveu <- svd(do.call(rbind, jivesol$joint))$v

spge <- nsprcomp::nsprcomp(t(rbind(x[[1]], x[[2]])), ncomp=6, k=120)
spmeth <- nsprcomp::nsprcomp(t(rbind(ifelse(is.na(x[[3]]), 0, x[[3]]),
  ifelse(is.na(x[[4]]), 0, x[[4]]))), ncomp=6, k=230)

for (i in 1:6) {
  pc <- paste('GEPC', i, sep='')
  clinical[[pc]] <- c(sge$u[i, ], rep(NA, nrow(x[[4]])))
  pc <- paste('MethPC', i, sep='')
  clinical[[pc]] <- c(rep(NA, nrow(x[[1]])), smeth$u[i, ])
  pc <- paste('SpGEPC', i, sep='')
  clinical[[pc]] <- c(spge$rotation[, i], rep(NA, nrow(x[[4]])))
  pc <- paste('SpMethPC', i, sep='')
  clinical[[pc]] <- c(rep(NA, nrow(x[[1]])), spmeth$rotation[, i])
  pc <- paste('JIVEPC', i, sep='')
  clinical[[pc]] <- c(rep(NA, nrow(x[[1]])), jiveu[i, ], rep(NA, nrow(x[[4]])))
}

gcimp <- as.numeric(clinical$gcimp)
pdf('fig/compare-pca.pdf', width=15, height=15)
pairs(clinical[, c('MVPC1', 'MVPC2', 'MVPC3', 'MVPC4', 'MVPC5', 'MVPC6')],
  col=clinical$subtype, pch=gcimp)
pairs(clinical[, c('SpMethPC1', 'SpMethPC2', 'SpMethPC3', 'SpMethPC4', 'SpMethPC5', 'SpMethPC6')],
  col=clinical$subtype, pch=gcimp)
pairs(clinical[, c('MethPC1', 'MethPC2', 'MethPC3', 'MethPC4', 'MethPC5', 'MethPC6')],
  col=clinical$subtype, pch=gcimp)
pairs(clinical[, c('SpGEPC1', 'SpGEPC2', 'SpGEPC3', 'SpGEPC4', 'SpGEPC5', 'SpGEPC6')],
  col=clinical$subtype, pch=gcimp)
pairs(clinical[, c('GEPC1', 'GEPC2', 'GEPC3', 'GEPC4', 'GEPC5', 'GEPC6')],
  col=clinical$subtype, pch=gcimp)
pairs(clinical[, c('JIVEPC1', 'JIVEPC2', 'JIVEPC3', 'JIVEPC4', 'JIVEPC5', 'JIVEPC6')],
  col=clinical$subtype, pch=gcimp)
dev.off()

D <- res$mvpca_res$solution$D
r2dir <- matrix(NA, 6, 6)
for (i in 1:6) {
  for (j in 1:6) {
    active <- 0 != row_apply(D[, res$inds[i, ]], prod)
    r2dir[i, j] <- sum(row_apply(D[active, res$inds[j, ]], prod)^2) /
      sum(res$x[[j]]^2, na.rm=TRUE)
  }
}

#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length = n + 1)
#  hcl(h = hues, l = 65, c = 100, alpha=0.625)[1:n]
#}
#pdf('fig/graph.pdf', width=3, height=3)
#par(mar=rep(0, 4))
#wmat <- r2dir[1:4, 1:4]
#wmat[wmat < 0.35] <- 0
#g <- igraph::graph.adjacency(wmat^4, mode='directed', weighted=TRUE, diag=FALSE)
#blocknames <- c(
#    expression('GE in'~C[1]),
#    expression('GE in'~C[2]),
#    expression('Meth. in'~C[2]),
#    expression('Meth. in'~C[3]))
#igraph::V(g)$name <- blocknames
#igraph::V(g)$shape <- c('circle', 'circle', 'square', 'square')
#igraph::V(g)$color <- gg_color_hue(3)[c(1, 2, 2, 3)]
#l <- matrix(NA, 4, 2)
#l[1, ] <- c(-1, 1)
#l[2, ] <- c(-1, -1/3)
#l[3, ] <- c(1, 1/3)
#l[4, ] <- c(1, -1)
#l[, 1] <- 0.75*l[, 1]
#plot(g, edge.width=100*igraph::E(g)$weight, edge.curved=TRUE,
#  edge.arrow.size=1.00, vertex.label.color='black', vertex.label.family='sans',
#  vertex.size=50, vertex.label.cex=0.5, layout=l, rescale=FALSE,
#  xlim=c(-1.5, 1.5))
#dev.off()

stop('stop')

pdf('fig/graphmatrix.pdf', width=1.5, height=1.5)
par(mar=c(1, 1, 8, 0))
par(mfrow=c(1, 1), cex=0.5)
imm(r2dir[1:4, 1:4])
axis(2, 1:4, c(
    expression('GE in'~C[1]),
    expression('GE in'~C[2]),
    expression('Meth. in'~C[2]),
    expression('Meth. in'~C[3])),
  line=-4.5, tick=FALSE, las=2)
axis(3, 1:4, c(
    expression('GE in'~C[1]),
    expression('GE in'~C[2]),
    expression('Meth. in'~C[2]),
    expression('Meth. in'~C[3])),
  line=-0.5, tick=FALSE, las=2)
title(expression('B) Directed '*R^2*' matrix'), adj=0, line=6)
dev.off()

pdf('fig/ramp.pdf', width=0.75*3/5, height=2.5*3/5)
par(mar=c(2.5, 0, 8, 3.75), cex=0.5)
minval <- min(r2dir[1:4, 1:4])
maxval <- max(r2dir[1:4, 1:4])
colfunc <- viridisLite::viridis
legend_image <- as.raster(matrix(rev(colfunc(256)), ncol=1))
plot(c(0, 1), c(minval, maxval), type='n', xaxs='i', yaxs='i', axes=FALSE,
  xlab='')
axis(4, c(0.2, 0.3, 0.4))
rasterImage(legend_image, 0, minval, 1, maxval)
dev.off()


# patient and variable clusters
obsord <- mvpca_order_subset(res$inds, res$mvpca_res, 1:3, 1:4)
varord <- mvpca_order_subset(res$inds, res$mvpca_res, 4:5, 1:4)

x <- res$x
obstype <- c(rep(1, nrow(x[[1]])), rep(2, nrow(x[[2]])), rep(3, nrow(x[[4]])))
obstype <- obstype[obsord]
vartype <- c(rep(1, ncol(x[[1]])), rep(2, ncol(x[[3]])))
vartype <- vartype[varord]
#imm(t(t(obstype)), colors=c('#000000', '#BBBBBB', '#FFFFFF'))
#imm(t(vartype), colors=c('#000000', '#FFFFFF'))

xhat <- function(V, D, component_subset, row, col) {
  x1 <- V[[row]][, component_subset, drop=FALSE] %*%
    diag(D[component_subset, row] * D[component_subset, col],
      length(component_subset)) %*%
    t(V[[col]][, component_subset, drop=FALSE])
  probit(x1)
}
V <- res$mvpca_res$solution$V
D <- res$mvpca_res$solution$D
bigx <- matrix(NA, length(obsord), length(varord))
colnames(bigx) <- c(colnames(x[[1]]), colnames(x[[3]]))
rownames(bigx) <- c(rownames(x[[1]]), rownames(x[[3]]), rownames(x[[4]]))
bigx[1:nrow(x[[1]]), 1:ncol(x[[1]])] <- xhat(V, D, 1:11, 1, 4)
bigx[1:nrow(x[[2]]) + nrow(x[[1]]), 1:ncol(x[[2]])] <- xhat(V, D, 1:11, 2, 4)
bigx[1:nrow(x[[2]]) + nrow(x[[1]]), 1:ncol(x[[3]]) + ncol(x[[1]])] <-
  xhat(V, D, 1:11, 2, 5)
bigx[1:nrow(x[[4]]) + nrow(x[[1]]) + nrow(x[[2]]),
  1:ncol(x[[3]]) + ncol(x[[1]])] <- xhat(V, D, 1:11, 3, 5)
#pdf('modelintegrated.pdf', height=11.7, width=16.5)
#par(mar=c(0, 1, 0, 1))
#imm(bigx[obsord, varord], na.hl=FALSE)
#dev.off()

databigx <- bigx
databigx[1:nrow(x[[1]]), 1:ncol(x[[1]])] <- probit(x[[1]])
databigx[1:nrow(x[[2]]) + nrow(x[[1]]), 1:ncol(x[[2]])] <- probit(x[[2]])
databigx[1:nrow(x[[2]]) + nrow(x[[1]]), 1:ncol(x[[3]]) + ncol(x[[1]])] <-
probit(x[[3]])
databigx[1:nrow(x[[4]]) + nrow(x[[1]]) + nrow(x[[2]]),
  1:ncol(x[[3]]) + ncol(x[[1]])] <- probit(x[[4]])
#pdf('dataintegrated.pdf', height=11.7, width=16.5)
#par(mar=c(0, 1, 0, 1))
#imm(databigx[obsord, varord], na.hl=FALSE)
#dev.off()

bigx[1:nrow(x[[1]]), 1:ncol(x[[3]]) + ncol(x[[1]])] <- xhat(V, D, 1:11, 1, 5)
bigx[1:nrow(x[[4]]) + nrow(x[[1]]) + nrow(x[[2]]), 1:ncol(x[[1]])] <-
  xhat(V, D, 1:11, 3, 4)
imputed <- matrix(FALSE, nrow(bigx), ncol(bigx))
imputed[is.na(databigx)] <- TRUE


require("viridis")
require("colorspace")

source("src/bioex/utils2.R")
source("src/bioex/heatmap.4f.R")

Xplot <- bigx[obsord, varord]
Ximp <- imputed[obsord, varord]



obscl <- mvpcaclust_subset(res$x, res$inds, res$mvpca_res, obsord, 1:3, 1:4)
obsclid <- cutree(obscl, k=6) #list(cutree(obscl[[2]], k=10) == 7, cutree(mvcl[[3]], k=16) == 7)
obsclid2 <- cutree(obscl, k=10) #list(cutree(obscl[[2]], k=10) == 7, cutree(mvcl[[3]], k=16) == 7)

varcl <- mvpcaclust_subset(res$x, res$inds, res$mvpca_res, varord, 4:5, 1:4)
varclid <- cutree(obscl, k=6) #list(cutree(obscl[[2]], k=10) == 7, cutree(mvcl[[3]], k=16) == 7)

patients <- c(rownames(res$x[[1]]), rownames(res$x[[2]]), rownames(res$x[[4]]))
clinical <- TCGA$clinical[patients, ]
clinical$IDH <- ifelse(with_rows(TCGA$mutation, patients)[, 'IDH1'], 3,2)
clinical$IDH[is.na(clinical$IDH)] <- NA
levels(clinical$subtype)[1] <- NA
levels(clinical$gcimp)[1] <- NA


subtype <- clinical$subtype
IDHstat <- clinical$IDH
gcimp <- as.numeric(clinical$gcimp)
gcimp[gcimp == 1] <- 3

cohort <- c(rep(1, nrow(x[[1]])), rep(2, nrow(x[[2]])), rep(3, nrow(x[[4]])))

maxVal <- max(Xplot,na.rm=T)
minVal <- min(Xplot,na.rm=T)

steps <- c(minVal,seq(minVal, maxVal,0.01), maxVal)
totalScale <-  viridis(length(steps)-1) #c(colorRampPalette(c('#386CB0',"white"), space = "rgb")(sum(steps<meanVal)-1),'#ffffff',colorRampPalette(c('white',"#e41a1c"), space = "rgb")(sum(steps>meanVal)-1))

hexColorRamp <- function(data) {
  fun <- colorRamp(rev(viridisLite::inferno(256)))
  apply(data, 1:2, function(x) {
    if (is.na(x)) return(NA)
    rgbmat <- fun(x) / 256
    rgb(rgbmat[1], rgbmat[2], rgbmat[3])
  })
}

obscl1 <- cutree(obscl, k=8)[obsord] #s=100
obscl2 <- cutree(obscl, k=17)[obsord] #s=65

cmap <- unique(obscl1)
dtc1 <- sapply(obscl1, function(i) 1 + (which(cmap == i) %% 9))
cmap <- unique(obscl2)
dtc2 <- sapply(obscl2, function(i) 1 + (which(cmap == i) %% 9))
rCclust <- apply(rbind(dtc1, dtc2), 1:2, color_function)
rC <- apply(rbind(cohort, subtype, gcimp), 1:2, color_function)
rCcont <- rbind((clinical$age-min(clinical$age, na.rm=T))/max(clinical$age, na.rm=T),
  (clinical$survival-min(clinical$survival, na.rm=T))/max(clinical$survival, na.rm=T))
rC <- rbind(rCclust, rbind(rC, hexColorRamp(rCcont))[, obsord])
rownames(rC) <- c('Clusters (rank 2)', 'Clusters (rank 3)', 'Cohort', 'Subtype', 'GCIMP', 'Age', 'Survival')

varcl1 <- cutree(varcl, k=9)[varord] #s=100
varcl2 <- cutree(varcl, k=22)[varord] #s=65

cmap <- unique(varcl1)
dtc1 <- sapply(varcl1, function(i) color_function(1 + (which(cmap == i) %% 9)))
cmap <- unique(varcl2)
dtc2 <- sapply(varcl2, function(i) color_function(1 + (which(cmap == i) %% 9)))

blue <- color_function(2)
cC <- matrix(c(rep(blue, ncol(x[[1]])), rep('#FFFFFF', ncol(x[[3]])))[varord], ncol=1)
cC <- cbind(cC, dtc2, dtc1)
colnames(cC) <- rev(c('Clusters (rank 2)', 'Clusters (rank 3)', 'Gene / Meth.'))

#assoc <- readr::read_tsv('~/data/ldgcData/OTHER/methylation-gene.tsv')
#replist <- list()
#for (gene in colnames(x[[1]])) {
#  probes <- intersect(colnames(x[[3]]), unique(subset(assoc, UCSC_RefGene_Name == gene)$Name))
#  if (length(probes) > 1) {
#    replist[[length(replist)+1]] <- c(gene, probes)
#  }
#}

source("src/bioex/heatmap.4f.R")
pdf("fig/heatmap_v2.pdf",height=20,width=20)
par(mar=c(4,10,4,4),xpd=F)
heatmap.4f(Xplot,makeGray=Ximp,RowSideColors=rC,ColSideColors=cC,NumColSideColors=0.5,Rowv=F,Colv=F,dendrogram="none",scale="none",margins=c(30,10), symbreaks=FALSE, key=F, symkey=FALSE, density.info="none", denscol="black", trace="none",main="", cexRow=0.4,cexCol=1,col=totalScale, breaks=steps, KeyValueName="", labRow=F, labCol=F)
dev.off()

# heatmap_with_imputation(bigx[obsord, varord], imputed[obsord, varord])


#pdf('modelintegratedimputed.pdf', height=11.7, width=16.5)
#par(mar=c(0, 1, 0, 1))
#imm(bigx[obsord, varord], na.hl=FALSE)
#dev.off()

databigx[is.na(databigx)] <- bigx[is.na(databigx)]
#pdf('dataintegratedimputed.pdf', height=11.7, width=16.5)
#par(mar=c(0, 1, 0, 1))
#imm(databigx[obsord, varord], na.hl=FALSE)
#dev.off()

obscl <- mvpcaclust_subset(res$x, res$inds, res$mvpca_res, obsord, 1:3, 1:4)
#pdf('obstree.pdf', height=2, width=11.7*length(obsord)/384)
#par(mar=c(1, 0, 3, 0))
#plot(obscl, main='', axes=F, ylab='', xlab='', sub='', hang=-1, cex=0.2,
#  labels=FALSE)
#dev.off()

varcl <- mvpcaclust_subset(res$x, res$inds, res$mvpca_res, varord, 4:5, 1:4)
#pdf('vartree.pdf', height=2, width=11.7*length(varord)/384)
#par(mar=c(1, 0, 3, 0))
#plot(varcl, main='', axes=F, ylab='', xlab='', sub='', hang=-1, cex=0.2,
#  labels=FALSE)
#dev.off()

# variables <- c(colnames(res$x[[1]]), colnames(res$x[[3]]))
# Uvar <- do.call(rbind, U[4:5])
# df <- list()
# for (i in 1:24) {
#   pc <- paste('PC', i, sep='')
#   df[[pc]] <- Uvar[, i]
# }
# df <- as.data.frame(df)
# rownames(df) <- variables
# geneset <- list() # ordered by length:
# geneset$mtor <- read.csv('~/data/KEGG_pathways/map04150_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$ras <- read.csv('~/data/KEGG_pathways/map04014_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$mapk <- read.csv('~/data/KEGG_pathways/map04010_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$pi3kakt <- read.csv('~/data/KEGG_pathways/map04151_genelist',
#   stringsAsFactors=FALSE)[, 1]
# for (i in 1:length(geneset)) {
#   df[names(geneset)[i]] <- variables %in% geneset[[i]]
# }
# df$type <- c(rep('gene', ncol(res$x[[1]])), rep('site', ncol(res$x[[3]])))
# 
# ggplot(df, aes(x=PC2, y=PC4, color=type)) + geom_point() + theme_bw()
# 
# pnum <- matrix(NA, 384, 9)
# pnum[, 1] <- clinical$age
# pnum[clinical$sex == 'MALE', 2] <- 1
# pnum[clinical$sex == 'FEMALE', 2] <- 2
# pnum[, 3] <- clinical$vitalstat
# pnum[, 4] <- clinical$survival
# pnum[, 5] <- clinical$IDH
# pnum[, 6:9] <- 0
# pnum[clinical$subtype == 'CL', 6] <- 1
# pnum[clinical$subtype == 'MS', 7] <- 1
# pnum[clinical$subtype == 'NL', 8] <- 1
# pnum[clinical$subtype == 'PN', 9] <- 1
# normimm <- function(x) {
#   x <- x - min(x, na.rm=TRUE)
#   x <- 2*x/max(x, na.rm=TRUE)
#   return(x-1)
# }
# col_apply <- function(x, f) apply(x, 2, f)
# pnum <- col_apply(pnum, normimm)
# 
# ix <- c(Vorder[[1]], Vorder[[2]]+nrow(res$x[[1]]),
#   Vorder[[3]]+nrow(res$x[[1]])+nrow(res$x[[2]]))
# imm(pnum[ix, ])
# 
# gnum <- matrix(0, 4, ncol(res$x[[1]]))
# 
# geneset <- list() # ordered by length:
# geneset$mtor <- read.csv('~/data/KEGG_pathways/map04150_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$ras <- read.csv('~/data/KEGG_pathways/map04014_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$mapk <- read.csv('~/data/KEGG_pathways/map04010_genelist',
#   stringsAsFactors=FALSE)[, 1]
# geneset$pi3kakt <- read.csv('~/data/KEGG_pathways/map04151_genelist',
#   stringsAsFactors=FALSE)[, 1]
# 
# pw <- colnames(res$x[[1]])
# gnum[1, pw %in% geneset[[1]]] <- 1
# gnum[2, pw %in% geneset[[2]]] <- 1
# gnum[3, pw %in% geneset[[3]]] <- 1
# gnum[4, pw %in% geneset[[4]]] <- 1
# 
# ix <- Vorder[[4]]
# imm(gnum[, ix])
# 
# x <- res$x
# inds <- res$inds
# result <- res$mvpca_res
# reorder <- TRUE
# 
#   zero2na <- function(x) {
#     x[x == 0] <- NA
#     return(x)
#   }
# 
#   imm <- function(x, xlab='column index', ylab='row index', asp=TRUE,
#       legend.outside=NA, ...) {
#     if (is.na(legend.outside)) {
#       legend.outside <- !asp
#     }
#     if (asp) {
#       asp <- 1
#     } else {
#       asp <- NA
#     }
#     d <- dim(x)+1
#     colors <- viridisLite::viridis(256)
#     image(1:d[2]-0.5, 1:d[1]-0.5, t(x), xlim=c(1, d[2])-0.5,
#       ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
#       col=colors, ...)
#     axis(1, pretty(1:d[2]), pos=d[1]-0.5)
#     axis(2, pretty(1:d[1]), pos=0.5)
#     # highlight NA, NaN, Inf
#     ix <- which(is.na(t(x)), arr.ind=TRUE)
#     rect(ix[, 1]-0.5, ix[, 2]-0.5, ix[, 1]+0.5, ix[, 2]+0.5, border='red')
#     ix <- which(is.infinite(t(x)), arr.ind=TRUE)
#     rect(ix[, 1]-0.5, ix[, 2]-0.5, ix[, 1]+0.5, ix[, 2]+0.5, border='red')
#     maxv <- max(x[!is.infinite(x)], na.rm=TRUE)
#     minv <- min(x[!is.infinite(x)], na.rm=TRUE)
#     delta <- (maxv-minv)/4
#     legs <- c(maxv, maxv-delta, maxv-2*delta, maxv-3*delta, minv)
#     legcols <- colors[c(256, 192, 128, 64, 1)]
#     if (legend.outside) {
#       par(mar=c(5.1, 4.1, 4.1, 5.1))
#       legend('topright', legend=signif(legs, 2), col=legcols, pt.bg=legcols,
#         pch=22, inset=c(-0.15, 0), xpd=TRUE)
#     } else {
#       legend('topright', legend=signif(legs, 2), col=legcols, pt.bg=legcols,
#         pch=22)
#     }
#   }
# 
#   probit <- function(x) pnorm(x, mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
#   probitxy <- function(x, y) {
#     if (all(is.na(y))) {
#       pnorm(x, mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
#     } else {
#       pnorm(x, mean(y, na.rm=TRUE), sd(y, na.rm=TRUE))
#     }
#   }
#   last <- function(x) x[length(x)]
# 
#   n <- max(inds)
#   V <- result$solution$V
#   D <- result$solution$D
#   xhat <- result$solution$xhat
# 
#   # find order for rows of each V
#   Vorder <- lapply(1:n, function(j) {
#    v <- V[[j]]
#    count <- sum(D[, j] > 0)
#    ix <- order(-D[, j])[1:count]
#    do.call(order, c(lapply(ix, function(i) {
#      vv <- v[, i]
#      vv[v[, i] < 0] <- 1
#      vv[v[, i] == 0] <- 0
#      vv[v[, i] > 0] <- -1
#      return(vv)
#    }), list(v[, ix[1]])))
#   })
# 
#   if (!reorder) {
#     Vorder <- lapply(1:n, function(j) 1:nrow(V[[j]]))
#   }
# 
#   # find position of each block
#   ixrow <- list()
#   rowcount <- 0
#   for (i in unique(inds[, 1])) {
#     ixrow[[i]] <- 1:nrow(V[[i]])+rowcount
#     rowcount <- rowcount + nrow(V[[i]])
#   }
#   ixcol <- list()
#   colcount <- 0
#   for (i in unique(inds[, 2])) {
#     ixcol[[i]] <- 1:nrow(V[[i]])+colcount
#     colcount <- colcount + nrow(V[[i]])
#   }
#   rowborders <- sapply(-sort(-unique(inds[, 1]))[-1],
#     function(i) last(ixrow[[i]]))
#   if (length(rowborders) == 0) {
#     rowborders <- -1
#   }
#   colborders <- sapply(-sort(-unique(inds[, 2]))[-1],
#     function(i) last(ixcol[[i]]))
#   if (length(colborders) == 0) {
#     colborders <- -1
#   }
# 
#   # construct matrices of all blocks for data, xhat and rank 1 xhat
#   blocks <- matrix(NA, rowcount, colcount)
#   for (i in 1:length(x)) {
#     blocks[ixrow[[inds[i, 1]]], ixcol[[inds[i, 2]]]] <-
#       probit(x[[i]][Vorder[[inds[i, 1]]], Vorder[[inds[i, 2]]]])
#   }
#   blockshat <- matrix(NA, rowcount, colcount)
#   for (i in 1:length(x)) {
#     row <- inds[i, 1]
#     col <- inds[i, 2]
#     xhat <- V[[row]] %*% diag(D[, row] * D[, col]) %*% t(V[[col]])
#     blockshat[ixrow[[row]], ixcol[[col]]] <-
#       probitxy(xhat[Vorder[[row]], Vorder[[col]]],
#         x[[i]][Vorder[[row]], Vorder[[col]]])
#   }
# 
#   blocks1 <- matrix(NA, rowcount, colcount)
#   for (i in 1:length(x)) {
#     row <- inds[i, 1]
#     col <- inds[i, 2]
#     x1 <- V[[row]][, component_subset, drop=FALSE] %*%
#       diag(D[component_subset, row] * D[component_subset, col],
#         length(component_subset)) %*%
#       t(V[[col]][, component_subset, drop=FALSE])
#     blocks1[ixrow[[row]], ixcol[[col]]] <-
#       probitxy(x1[Vorder[[row]], Vorder[[col]]],
#         x[[i]][Vorder[[inds[i, 1]]], Vorder[[inds[i, 2]]]])
#   }
# 
#   # prepare Vs for plotting
#   vv <- lapply(1:n, function(i) {
#     u <- V[[i]] %*% diag(D[, i])
#     u[u == 0] <- NA
#     return(u)
#   })
# 
#   for (i in 1:length(x)) {
#     row <- inds[i, 1]
#     col <- inds[i, 2]
#     if (!is.null(rownames(x[[i]]))) {
#       rownames(vv[[row]]) <- rownames(x[[i]])
#     }
#     if (!is.null(colnames(x[[i]]))) {
#       rownames(vv[[col]]) <- colnames(x[[i]])
#     }
#   }
# 
#   imm(zero2na(abs(D)), 'View', 'Rank')
#   imm(zero2na(abs(result$solution$R2_blockwise)), 'Block', 'Rank')
#   for (i in 1:n) {
#     imm(vv[[i]][Vorder[[i]], , drop=FALSE], 'Rank')
#   }
#   imm(blocks)
#   abline(h=0.5+rowborders, v=0.5+colborders)
#   imm(blockshat)
#   abline(h=0.5+rowborders, v=0.5+colborders)
#   imm(blocks1)
#   abline(h=0.5+rowborders, v=0.5+colborders)
#   for (i in 1:n) {
#     vvv <- rowSums(vv[[i]][, component_subset, drop=FALSE])
#     if (all(is.na(vvv))) next
#     ix <- 1:nrow(vv[[i]])
#     if (nrow(vv[[i]]) > 100) {
#       tmp <- -abs(vvv[rev(Vorder[[i]])])
#       tmp[is.na(tmp)] <- 0
#       val <- sort(tmp)[100]
#       ix <- which(tmp < val)
#     }
#     par(mai=c(1,2,1,1))
#     mp <- barplot(vvv[rev(Vorder[[i]])][ix], horiz=TRUE,
#       names.arg=rownames(vvv[rev(Vorder[[i]])])[ix], las=1,
#       cex.names=0.45)
#     abline(h=mp, lty=3, col='gray')
#   }
