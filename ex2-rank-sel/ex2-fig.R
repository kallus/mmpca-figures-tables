row_apply <- function(x, f) t(apply(x, 1, f))
col_apply <- function(x, f) apply(x, 2, f)
normalize <- function(x) x/sqrt(sum(x^2))
norm <- function(x) sqrt(sum(x^2))
angle <- function(x, y) {
  a <- acos(sum(x*y)/sqrt(sum(x^2)*sum(y^2)))
  min(a, pi-a)
}
angle_to <- function(v) {
  function(x) angle(x, v)
}
joint_dir <- function(V, D, x) {
  Dtmp <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  if (all(row_apply(Dtmp, min) == 0) && all(Dtmp[4, ] == 0)) {
    # return pc1 of residual
    residual <- lapply(1:3, function(i) x[[i]] - V[[i]] %*%
      diag(D[, i] * D[, 4]) %*% t(V[[4]]))
    return(svd(do.call(rbind, residual))$v[, 1])
  }
  if (all(row_apply(Dtmp, min) == 0)) {
    # look at second smallest
    return(V[[4]][, which.max(row_apply(Dtmp, function(x) sort(x)[2]))])
  }
  V[[4]][, which.max(row_apply(Dtmp, min))]
}
joint_dir_jive <- function(jivesol, x) {
  if (jivesol$rankJ == 0 || sum(abs(do.call(rbind, jivesol$joint))) == 0) {
    # return pc1 of residual
    residual <- lapply(1:3, function(i) x[[i]] - jivesol$individual[[i]])
    return(svd(do.call(rbind, residual))$v[, 1])
  }
  svd(do.call(rbind, jivesol$joint))$v[, 1]
}
correct <- function(v, dirs) {
  if (4 == which.min(row_apply(dirs, angle_to(v)))) {
    return(1)
  } else {
    return(0)
  }
}

individual_comp_count <- function(D) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM != 0] <- 1
  sum(rowSums(DM) == 1)
}

individual_and_pair_comp_count <- function(D) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM != 0] <- 1
  sum(rowSums(DM) == 1) + sum(rowSums(DM) == 2)
}

joint_comp_count <- function(D) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM != 0] <- 1
  sum(rowSums(DM) == 3)
}

matrix_rank <- function(D, i) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM != 0] <- 1
  sum(DM[, i])
}

individual_comp_count_thr <- function(D, thr) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM <= thr] <- 0
  DM[DM > thr] <- 1
  sum(rowSums(DM) == 1)
}

individual_and_pair_comp_count_thr <- function(D, thr) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM <= thr] <- 0
  DM[DM > thr] <- 1
  sum(rowSums(DM) == 1) + sum(rowSums(DM) == 2)
}

joint_comp_count_thr <- function(D, thr) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM <= thr] <- 0
  DM[DM > thr] <- 1
  sum(rowSums(DM) == 3)
}

matrix_rank_thr <- function(D, i, thr) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM <= thr] <- 0
  DM[DM > thr] <- 1
  sum(DM[, i])
}

total_rank_thr <- function(D, thr) {
  DM <- D[, 1:3] * D[, 4] %*% t(rep(1, 3))
  DM[DM <= thr] <- 0
  DM[DM > thr] <- 1
  sum(rowSums(DM) > 0)
}

n <- 100
d <- 25
options(stringsAsFactors = F)
tbl <- data.frame(Seed=numeric(0), Jointprop=numeric(0), Method=character(0),
  Correct=numeric(0), Individual=numeric(0), IandP=numeric(0), Joint=numeric(0),
  Matrix1=numeric(0), Matrix2=numeric(0), Matrix3=numeric(0), Rank=numeric(0))
for (seed in 1:500) {
  for (jointprop in seq(0, 1, 0.0005)) {
    filename <- paste('cache/ex2', jointprop, '-', seed, '-', n, '-', d, '.rds.gz', sep='')
    if (file.exists(filename)) {
      res <- readRDS(filename)
      # add noise dir
      dirs <- rbind(res$indata$dirs, svd(do.call(rbind, res$indata$noisemats))$v[, 1])
      # MM-PCA
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='MM-PCA',
        Correct=correct(joint_dir(res$MMPCA$V, abs(res$MMPCA$D), res$indata$x), dirs),
        Individual=individual_comp_count(res$MMPCA$D),
        IandP=individual_and_pair_comp_count(res$MMPCA$D),
        Joint=joint_comp_count(res$MMPCA$D),
        Matrix1=matrix_rank(res$MMPCA$D, 1),
        Matrix2=matrix_rank(res$MMPCA$D, 2),
        Matrix3=matrix_rank(res$MMPCA$D, 3),
        Rank=sum(rowSums(abs(res$MMPCA$D)) > 0))
      # CMF
      cmfV <- lapply(res$CMF$U, function(U) col_apply(U, normalize))
      cmfD <- do.call(cbind, lapply(res$CMF$U, function(U) col_apply(U, norm)))
      thr <- 0.1
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='CMF',
        Correct=correct(joint_dir(cmfV, cmfD, res$indata$x), dirs),
        Individual=individual_comp_count_thr(cmfD, thr),
        IandP=individual_and_pair_comp_count_thr(cmfD, thr),
        Joint=joint_comp_count_thr(cmfD, thr),
        Matrix1=matrix_rank_thr(cmfD, 1, thr),
        Matrix2=matrix_rank_thr(cmfD, 2, thr),
        Matrix3=matrix_rank_thr(cmfD, 3, thr),
        Rank=total_rank_thr(cmfD, thr))
      # PCA
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='PCA',
        Correct=correct(svd(do.call(rbind, res$indata$x))$v[, 1], dirs),
        Individual=NA,
        IandP=NA,
        Joint=NA,
        Matrix1=NA,
        Matrix2=NA,
        Matrix3=NA,
        Rank=NA)
      # JIVE
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='JIVE',
        Correct=correct(joint_dir_jive(res$JIVE, res$indata$x), dirs),
        Individual=sum(res$JIVE$rankA),
        IandP=0,
        Joint=res$JIVE$rankJ,
        Matrix1=res$JIVE$rankA[1] + res$JIVE$rankJ,
        Matrix2=res$JIVE$rankA[2] + res$JIVE$rankJ,
        Matrix3=res$JIVE$rankA[3] + res$JIVE$rankJ,
        Rank=res$JIVE$rankJ + sum(res$JIVE$rankA))
      # JIVE-OC
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='JIVE-OC',
        Correct=correct(joint_dir_jive(res$JIVEO, res$indata$x), dirs),
        Individual=sum(res$JIVEO$rankA),
        IandP=0,
        Joint=res$JIVEO$rankJ,
        Matrix1=res$JIVEO$rankA[1] + res$JIVEO$rankJ,
        Matrix2=res$JIVEO$rankA[2] + res$JIVEO$rankJ,
        Matrix3=res$JIVEO$rankA[3] + res$JIVEO$rankJ,
        Rank=res$JIVEO$rankJ + sum(res$JIVEO$rankA))
    }
  }
}

library(dplyr)
library(ggplot2)
pdf('fig/ex2a.pdf', width=3, height=2.5)
plotdf <- subset(tbl, Method != 'Random') %>%
  group_by(Method, Jointprop) %>%
  summarize(Accuracy=mean(Correct), se=2*sqrt(Accuracy*(1-Accuracy)/n()))
print(ggplot(plotdf, aes(x=Jointprop, y=100*Accuracy, ymin=100*(Accuracy-se),
    ymax=100*(Accuracy+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  ggtitle('A)') + xlab('Joint proportion') + ylab('Accuracy (%)') +
  expand_limits(y=c(0, 100), x=c(1/40, 20/40)) +
  theme_bw() + theme(legend.position='none'))
dev.off()
pdf('fig/ex2c.pdf', width=3, height=2.5)
plotdf <- group_by(tbl, Method, Jointprop) %>%
  summarize(Mean=mean(Joint), se=sd(Joint)/sqrt(n()))
print(ggplot(plotdf, aes(x=Jointprop, y=Mean, ymin=(Mean-se),
    ymax=(Mean+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  geom_hline(yintercept=1, linetype=2) +
  ggtitle('C)') + xlab('Joint proportion') + ylab('Global components') +
  theme_bw() + theme(legend.position='none'))
dev.off()
pdf('fig/ex2e.pdf', width=3, height=2.5)
plotdf <- group_by(tbl, Method, Jointprop) %>%
  summarize(Mean=mean(Matrix1), se=sd(Matrix1)/sqrt(n()))
print(ggplot(plotdf, aes(x=Jointprop, y=Mean, ymin=(Mean-se),
    ymax=(Mean+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  geom_hline(yintercept=2, linetype=2) +
  ggtitle('E)') + xlab('Joint proportion') + ylab('Matrix components') +
  theme_bw() + theme(legend.position='none'))
dev.off()

n <- 10
d <- 40
options(stringsAsFactors = F)
tbl <- data.frame(Seed=numeric(0), Jointprop=numeric(0), Method=character(0),
  Correct=numeric(0), Individual=numeric(0), IandP=numeric(0), Joint=numeric(0),
  Matrix1=numeric(0), Matrix2=numeric(0), Matrix3=numeric(0), Rank=numeric(0))
for (seed in 1:500) {
  for (jointprop in seq(0, 1, 0.0005)) {
    filename <- paste('cache/ex2', jointprop, '-', seed, '-', n, '-', d, '.rds.gz', sep='')
    if (file.exists(filename)) {
      res <- readRDS(filename)
      # add noise dir
      dirs <- rbind(res$indata$dirs, svd(do.call(rbind, res$indata$noisemats))$v[, 1])
      # MM-PCA
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='MM-PCA',
        Correct=correct(joint_dir(res$MMPCA$V, abs(res$MMPCA$D), res$indata$x), dirs),
        Individual=individual_comp_count(res$MMPCA$D),
        IandP=individual_and_pair_comp_count(res$MMPCA$D),
        Joint=joint_comp_count(res$MMPCA$D),
        Matrix1=matrix_rank(res$MMPCA$D, 1),
        Matrix2=matrix_rank(res$MMPCA$D, 2),
        Matrix3=matrix_rank(res$MMPCA$D, 3),
        Rank=sum(rowSums(abs(res$MMPCA$D)) > 0))
      # CMF
      cmfV <- lapply(res$CMF$U, function(U) col_apply(U, normalize))
      cmfD <- do.call(cbind, lapply(res$CMF$U, function(U) col_apply(U, norm)))
      thr <- 0.1
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='CMF',
        Correct=correct(joint_dir(cmfV, cmfD, res$indata$x), dirs),
        Individual=individual_comp_count_thr(cmfD, thr),
        IandP=individual_and_pair_comp_count_thr(cmfD, thr),
        Joint=joint_comp_count_thr(cmfD, thr),
        Matrix1=matrix_rank_thr(cmfD, 1, thr),
        Matrix2=matrix_rank_thr(cmfD, 2, thr),
        Matrix3=matrix_rank_thr(cmfD, 3, thr),
        Rank=total_rank_thr(cmfD, thr))
      # PCA
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='PCA',
        Correct=correct(svd(do.call(rbind, res$indata$x))$v[, 1], dirs),
        Individual=NA,
        IandP=NA,
        Joint=NA,
        Matrix1=NA,
        Matrix2=NA,
        Matrix3=NA,
        Rank=NA)
      # JIVE
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='JIVE',
        Correct=correct(joint_dir_jive(res$JIVE, res$indata$x), dirs),
        Individual=sum(res$JIVE$rankA),
        IandP=0,
        Joint=res$JIVE$rankJ,
        Matrix1=res$JIVE$rankA[1] + res$JIVE$rankJ,
        Matrix2=res$JIVE$rankA[2] + res$JIVE$rankJ,
        Matrix3=res$JIVE$rankA[3] + res$JIVE$rankJ,
        Rank=res$JIVE$rankJ + sum(res$JIVE$rankA))
      # JIVE-OC
      tbl[nrow(tbl) + 1, ] <- list(seed, jointprop,
        Method='JIVE-OC',
        Correct=correct(joint_dir_jive(res$JIVEO, res$indata$x), dirs),
        Individual=sum(res$JIVEO$rankA),
        IandP=0,
        Joint=res$JIVEO$rankJ,
        Matrix1=res$JIVEO$rankA[1] + res$JIVEO$rankJ,
        Matrix2=res$JIVEO$rankA[2] + res$JIVEO$rankJ,
        Matrix3=res$JIVEO$rankA[3] + res$JIVEO$rankJ,
        Rank=res$JIVEO$rankJ + sum(res$JIVEO$rankA))
    }
  }
}

library(dplyr)
library(ggplot2)
pdf('fig/ex2b.pdf', width=4.27, height=2.5)
plotdf <- subset(tbl, Method != 'Random') %>%
  group_by(Method, Jointprop) %>%
  summarize(Accuracy=mean(Correct), se=2*sqrt(Accuracy*(1-Accuracy)/n()))
print(ggplot(plotdf, aes(x=Jointprop, y=100*Accuracy, ymin=100*(Accuracy-se),
    ymax=100*(Accuracy+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  ggtitle('B)') +
  xlab('Joint proportion') + ylab('Accuracy (%)') +
  expand_limits(y=c(0, 100), x=c(1/40, 20/40)) +
  theme_bw())
dev.off()
pdf('fig/ex2d.pdf', width=3, height=2.5)
plotdf <- group_by(tbl, Method, Jointprop) %>%
  summarize(Mean=mean(Joint), se=sd(Joint)/sqrt(n()))
print(ggplot(plotdf, aes(x=Jointprop, y=Mean, ymin=(Mean-se),
    ymax=(Mean+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  geom_hline(yintercept=1, linetype=2) +
  ggtitle('D)') +
  xlab('Joint proportion') + ylab('Global components') +
  theme_bw() + theme(legend.position='none'))
dev.off()
pdf('fig/ex2f.pdf', width=3, height=2.5)
plotdf <- group_by(tbl, Method, Jointprop) %>%
  summarize(Mean=mean(Matrix1), se=sd(Matrix1)/sqrt(n()))
print(ggplot(plotdf, aes(x=Jointprop, y=Mean, ymin=(Mean-se),
    ymax=(Mean+se), color=Method, shape=Method)) +
  geom_ribbon(alpha=0.2, color=NA) +
  geom_line() +
  geom_point(size=1.25) +
  geom_hline(yintercept=2, linetype=2) +
  ggtitle('F)') +
  xlab('Joint proportion') + ylab('Matrix components') +
  theme_bw() + theme(legend.position='none'))
dev.off()
