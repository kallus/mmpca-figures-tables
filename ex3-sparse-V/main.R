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

if (!file.exists('sparse-v.rds.gz')) {
  set.seed(1)
  reps <- 100
  snrs <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1)
  results <- matrix(NA, reps*length(snrs)*2, 3)
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
      x <- x + snr * frob(x) / frob(noise) * noise
      tc <- system.time(cres <- CMF(list(x), matrix(c(1, 2), 1, 2), 2, 20000))
      tm <- system.time(mres <- mmpca::mmpca(list(x), matrix(c(1, 2), 1, 2), 2, init_theta=F))
      i <- i + 1
      method[i] <- 'CMF'
      results[i, 1] <- mcc(rbind(u, v), rbind(cres$U[[1]], cres$U[[2]]))
      results[i, 2] <- snr
      results[i, 3] <- tc[3]
      i <- i + 1
      rank <- sum(rowSums(mres$solution$D != 0) > 1)
      if (rank > 0) {
        method[i] <- 'MMPCA'
        results[i, 1] <- mcc(rbind(u, v), rbind(mres$solution$V[[1]], mres$solution$V[[2]]))
        results[i, 2] <- snr
        results[i, 3] <- tm[3]
      }
    }
  }
  tbl <- data.frame(Method=method, SNR=results[, 2], MCC=results[, 1],
    Time=results[, 3])
  saveRDS(tbl, 'sparse-v.rds.gz')
} else {
  tbl <- readRDS('sparse-v.rds.gz')
}

library(ggplot2)
#print(ggplot(tbl, aes(SNR, MCC, color=Method)) + geom_point(position='jitter') +
  #geom_smooth(method='lm'))

plot_labeller <- function(variable, value) {
  return(paste(variable, value, sep=' = '))
}

#tbldiff <- tbl[tbl$Method == 'CMF', 2:4]
#tbldiff$MCC <- tbl[tbl$Method == 'MMPCA', 3] - tbldiff$MCC
#tbldiff <- tbldiff[tbldiff$SNR != 0.75, ]
#tbldiff <- tbldiff[tbldiff$SNR != 0.075, ]
#tbldiff <- tbldiff[tbldiff$SNR != 0.025, ]
#tbldiff <- tbldiff[tbldiff$SNR != 0.01, ]
#tbldiff$SNR <- as.factor(1 / tbldiff$SNR)
tbl <- tbl[tbl$SNR != 0.75, ]
tbl <- tbl[tbl$SNR != 0.075, ]
tbl <- tbl[tbl$SNR != 0.025, ]
tbl <- tbl[tbl$SNR != 0.01, ]
tbl <- tbl[!is.na(tbl$SNR), ]
levels(tbl$Method)[which(levels(tbl$Method) == 'MMPCA')] <- 'MM-PCA'
tbl$SNR <- 1 / tbl$SNR
pdf('mccdiff.pdf', width=4.3, height=2.12)
print(ggplot(tbl, aes(Method, MCC)) + 
  geom_boxplot(outlier.alpha=0, color='darkgray') +
  facet_grid(.~SNR, labeller=plot_labeller) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_jitter(width=0.2, height=0, size=0.1) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  xlab('') +
  ylab('Accuracy (MCC)'))
dev.off()
