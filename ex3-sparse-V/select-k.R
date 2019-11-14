CMF <- function(data, views, K, iter.max) {
  D <- rep(NA, max(views))
  for (i in 1:nrow(views)) {
    D[views[i, 1]] <- nrow(data[[i]])
    D[views[i, 2]] <- ncol(data[[i]])
  }
  data <- lapply(data, CMF::matrix_to_triplets)
  opts <- CMF::getCMFopts()
  opts$verbose <- 0
  opts$iter.max <- iter.max
  opts$useBias <- F
  CMF::CMF(data, views, K, rep('gaussian', length(data)), D, opts=opts)
}

frob <- function(x) {
  sqrt(sum(x^2))
}

mcc <- function(v, vhat) {
  vhat[abs(vhat) < 0.01] <- 0
  vhat[vhat != 0] <- 1
  v <- c(v[, 1], v[, 2])
  vh <- c(vhat[, 1], vhat[, 2])
  m1 <- mccr::mccr(v, vh)
  vh <- c(vhat[, 2], vhat[, 1])
  m2 <- mccr::mccr(v, vh)
  return(max(m1, m2))
}

col_apply <- function(x, f) apply(x, 2, f)
norm <- function(x) sqrt(sum(x^2))

if (!file.exists('sparse-v-rank.rds.gz')) {
  set.seed(1)
  reps <- 100
  snrs <- c(1, 2, 4, 10, 20)
  results <- matrix(NA, reps*length(snrs)*2, 4)
  method <- rep(NA, nrow(results))
  i <- 0
  for (snr in snrs) {
    print(snr)
    for (i_rep in 1:reps) {
      u <- matrix(0, 30, 2)
      v <- matrix(0, 30, 2)
      u[sample(30, 3), 1] <- 1
      v[sample(30, 3), 1] <- 1
      u[sample(30, 3), 2] <- 1
      v[sample(30, 3), 2] <- 1
      noise <- matrix(rnorm(30*30), 30, 30)
      x <- u %*% t(v) 
      x <- x + frob(x) / frob(noise) * noise / snr
      tc <- system.time(cres <- CMF(list(x), matrix(c(1, 2), 1, 2), 3, 20000))
      cmfD <- do.call(cbind, lapply(cres$U, function(U) col_apply(U, norm)))
      tm <- system.time(mres <- mmpca::mmpca(list(x), matrix(c(1, 2), 1, 2), 3, init_theta=F))
      i <- i + 1
      method[i] <- 'CMF'
      results[i, 1] <- mcc(rbind(u, v), rbind(cres$U[[1]], cres$U[[2]]))
      results[i, 2] <- snr
      results[i, 3] <- tc[3]
      results[i, 4] <- sum(rowSums(abs(cmfD) > 0.01) > 1)
      i <- i + 1
      method[i] <- 'MMPCA'
      results[i, 1] <- mcc(rbind(u, v), rbind(mres$solution$V[[1]], mres$solution$V[[2]]))
      results[i, 2] <- snr
      results[i, 3] <- tm[3]
      results[i, 4] <- sum(rowSums(mres$solution$D != 0) > 1)
    }
  }
  tbl <- data.frame(Method=method, SNR=results[, 2], MCC=results[, 1],
    Time=results[, 3], Rank=results[, 4])
  saveRDS(tbl, 'sparse-v-rank.rds.gz')
} else {
  tbl <- readRDS('sparse-v-rank.rds.gz')
}

plot_labeller <- function(variable, value) {
  return(paste(variable, value, sep=' = '))
}

tbl <- tbl[tbl$SNR != 40, ]
tbl <- tbl[tbl$SNR != 100, ]

levels(tbl$Method)[which(levels(tbl$Method) == 'MMPCA')] <- 'MM-PCA'
tbl <- tbl[tbl$Rank != 0, ]

library(ggplot2)
tbl$SNR <- as.factor(tbl$SNR)
pdf('rank.pdf', width=4.3, height=2.75)
print(ggplot(tbl, aes(x=Method, y=Rank)) + #geom_point(position='jitter') +
  geom_boxplot(outlier.alpha=0, color='darkgray') +
  facet_grid(.~SNR, labeller=plot_labeller) +
  geom_hline(yintercept=2, linetype=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -60, hjust = 0)) +
  coord_cartesian(ylim = c(1, 3)) +
  geom_jitter(width=0.2, height=0.02, size=0.1))
  #geom_point() +
  #geom_smooth(method='lm'))
dev.off()
