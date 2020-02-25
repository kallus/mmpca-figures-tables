row_apply <- function(x, f) t(apply(x, 1, f))
col_apply <- function(x, f) apply(x, 2, f)
normalize <- function(x) x/sqrt(sum(x^2))
norm <- function(x) sqrt(sum(x^2))
noise <- function(n, d) matrix(rnorm(n*d), n, d)
data <- function(n, dir) rnorm(n) %*% t(dir)
angle <- function(x, y) {
  a <- acos(sum(x*y)/sqrt(sum(x^2)*sum(y^2)))
  min(a, pi-a)
}
angle_to <- function(v) {
  function(x) angle(x, v)
}

err <- function(D) {
  D <- D[, 1:4] * (D[, 5] %*% t(rep(1, 4)))
  alt1 <- c(D[1, 1:2], D[2, 3:4])
  alt2 <- c(D[2, 1:2], D[1, 3:4])
  sqrt(min(mean(alt1^2), mean(alt2^2)))
}

# currently requires f to return a list of named objects
df_apply <- function(df, f) {
  res <- list()
  for (i in 1:nrow(df)) {
    cat(round(i/nrow(df)*100)); cat('%\r')
    res[[i]] <- do.call(f, lapply(1:ncol(df), function(j) df[i, j]))
  }
  names <- names(res[[1]])
  res <- lapply(1:length(res[[1]]), function(i) sapply(1:nrow(df),
    function(j) res[[j]][[i]]))
  names(res) <- names
  as.data.frame(res)
}

df_apply2 <- function(df, f) {
  res <- list()
  for (i in 1:nrow(df)) {
    cat(round(i/nrow(df)*100)); cat('%\r')
    res[[i]] <- do.call(f, lapply(1:ncol(df), function(j) df[i, j]))
    res[[i]] <- cbind(df[rep(i, nrow(res[[i]])), ], res[[i]])
  }
  do.call(rbind, res)
}

CMF <- function(data, views, K, trace) {
  D <- rep(NA, max(views))
  for (i in 1:nrow(views)) {
    D[views[i, 1]] <- nrow(data[[i]])
    D[views[i, 2]] <- ncol(data[[i]])
  }
  data <- lapply(data, CMF::matrix_to_triplets)
  opts <- CMF::getCMFopts()
  opts$verbose <- trace - 1
  opts$iter.max <- 20000
  CMF::CMF(data, views, K, rep('gaussian', length(data)), D, opts=opts)
}

simulate <- function(n, d, seed, snr) {
  set.seed(seed)

  # pick 2 directions not too similar
  while (TRUE) {
    dirs <- matrix(rnorm(2*d), 2, d)
    dirs <- row_apply(dirs, normalize)
    if (angle(dirs[1, ], dirs[2, ]) > pi/4) {
      break
    }
  }

  # generate data
  block2dir <- c(1, 1, 2, 2)
  xsignal <- lapply(1:4, function(i) {
    data(n, dirs[block2dir[i], ]) / sqrt(n)
  })
  xnoise <- lapply(1:4, function(i) {
    noise(n, d) / sqrt(n)
  })
  snrhat <- sqrt(mean(do.call(rbind, xsignal)^2)) /
    sqrt(mean(do.call(rbind, xnoise)^2))
  x <- lapply(1:4, function(i) {
    xsignal[[i]] + snrhat / snr * xnoise[[i]]
  })
  inds <- matrix(5, 4, 2)
  inds[, 1] <- 1:4
  mvpca <- mmpca::mmpca(x, inds, k=2)
  cmf <- CMF(x, inds, 2, 0)

  D <- list(
    MVPCA=mvpca$solution$D,
    CMF=do.call(cbind, lapply(cmf$U,
      function(U) col_apply(U, norm)))
  )

  xhat <- list(mvpca$solution$xhat,
    lapply(1:4, function(i) cmf$U[[inds[i, 1]]] %*% t(cmf$U[[inds[i, 2]]])))
  xhat <- lapply(xhat, function(x) do.call(rbind, x))

  xsignal <- do.call(rbind, xsignal)
  rmse <- sapply(1:2, function(i) {
    sqrt(mean((xhat[[i]]-xsignal)^2))
  })
  return(data.frame(Method=c('MVPCA', 'CMF'), Error=sapply(D, err), RMSE=rmse))
  #return(c(lapply(D, err)))
}

if (file.exists('cache/ex1-simulations.rds.gz')) {
  simulations <- readRDS('cache/ex1-simulations.rds.gz')
} else {
  simulations <- expand.grid(
    observations=10,#c(5, 10, 25),#, 50),
    variables=10,
    seed=1:100,
    SNR=c(1, 2, 3, 4, 5, 7, 10))
    #SNR=c(10, 5, 2, 1, 0.5, 0.2, 0.1))
    #noise=c(0.05, 0.1, 0.25))#, 0.5))
  #simulations <- cbind(simulations, df_apply(simulations, simulate))
  simulations <- df_apply2(simulations, simulate)
  #simulations <- tidyr::gather(simulations, 'Method', 'Error', MVPCA, CMF)
  fh <- gzfile('cache/ex1-simulations.rds.gz')
  saveRDS(simulations, fh)
  close(fh)
}

plot_labeller <- function(variable, value) {
  return(paste(variable, value, sep=' = '))
}

mysqrt_trans <- function() {
   scales::trans_new("mysqrt",
             transform = base::sqrt,
             inverse = function(x) ifelse(x<0, 0, x^2),
             domain = c(0, Inf))
}

levels(simulations$Method) <- c('CMF', 'MM-PCA')

library(ggplot2)
pdf('fig/ex1.pdf', width=6, height=2.0)
print(ggplot(simulations, aes(x=Method, y=RMSE)) +
  geom_boxplot(outlier.alpha=0, color='darkgray') +
  geom_jitter(width=0.2, height=0, size=0.1) + 
  facet_grid(.~SNR, labeller=plot_labeller) +
  theme_bw() +
  #scale_y_continuous(trans='mysqrt') +
  expand_limits(y=c(0, NA)) +
  xlab('') +
  #ggtitle('A) Simulation 1') +
  ylab('Loss (RMSE)') +
  theme(axis.text.x = element_blank()) +
  theme(plot.margin=unit(c(5.5, 10, -5.5, 5.5), 'points')))
dev.off()
pdf('fig/ex1b.pdf', width=6, height=2.75)
print(ggplot(simulations, aes(x=Method, y=Error)) +
  geom_boxplot(outlier.alpha=0, color='darkgray') +
  geom_jitter(width=0.2, height=0, size=0.1) + 
  facet_grid(.~SNR, labeller=plot_labeller) +
  theme_bw() +
  #scale_y_continuous(trans='mysqrt') +
  expand_limits(y=c(0, NA)) +
  ylab('Structure error (RMSE)') +
  #ylab(expression(Structure~error*' ('*l[2]*'-dist.)')) +
  theme(axis.text.x = element_text(angle = -60, hjust = 0)) +
  theme(plot.margin=unit(c(5.5, 10, 10, 10.5), 'points')))
  #theme(plot.margin=unit(c(5.5, 10, 10, 5.5), 'points')))
dev.off()
#pdf('fig/ex1-bars.pdf', width=1.5, height=2.25)
#library(dplyr)
#plotdf <- simulations %>%
#  group_by(Method, noise, observations) %>%
#  summarize(Exact=sum(Error==0)/n())
#print(plotdf)
#print(ggplot(plotdf, aes(x=Method, y=Exact)) +
#  geom_bar(stat='identity') +
#  expand_limits(y=c(0, 1)) +
#  theme_bw() +
#  ggtitle('B) Structure') + ylab('Accuracy') +
#  theme(axis.text.x = element_text(angle = -45, hjust = 0)))
#dev.off()
