row_apply <- function(x, f) t(apply(x, 1, f))
col_apply <- function(x, f) apply(x, 2, f)
normalize <- function(x) x/sqrt(sum(x^2))
norm <- function(x) sqrt(sum(x^2))
noise <- function(n, d, eps) matrix(rnorm(n*d, 0, eps), n, d)
data <- function(n, dir) rnorm(n) %*% t(dir)
angle <- function(x, y) {
  a <- acos(sum(x*y)/sqrt(sum(x^2)*sum(y^2)))
  min(a, pi-a)
}
angle_to <- function(v) {
  function(x) angle(x, v)
}

# currently requires f to return a list of named objects
df_apply <- function(df, f) {
  args <- lapply(1:nrow(df),
    function(i) lapply(1:ncol(df), function (j) df[i, j]))
  res <- parallel::mclapply(args, function(a) do.call(f, a),
    mc.preschedule=FALSE)
  names <- names(res[[1]])
  res <- lapply(1:length(res[[1]]), function(i) sapply(1:nrow(df),
    function(j) res[[j]][[i]]))
  names(res) <- names
  as.data.frame(res)
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

#simulate <- function(n, d, seed, eps, jointprop) {
n <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
d <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
i <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
seed <- i
eps <- 0.05
jointprop <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
ifilename <- paste('cache/ex2', jointprop, '-', i, '-', n, '-', d, '.rds.gz', sep='')
print(paste(n, d, seed, eps, jointprop))
if (file.exists(ifilename)) {
  stop('Already run')
}

set.seed(seed)

dirs <- matrix(0, 4, d)

# pick 3 directions pi/3 apart pairwise
dirs[1, 1:3] <- c(1, 1, 0)
dirs[2, 1:3] <- c(0, 1, 1)
dirs[3, 1:3] <- c(1, 0, 1)
dirs <- row_apply(dirs, normalize)

# pick one direction at least pi/3 away from the others
while (TRUE) {
  dirs[4, ] <- normalize(matrix(rnorm(d), 1, d))
  if (min(row_apply(dirs[1:3, ], angle_to(dirs[4, ]))) > pi/3) {
    break
  }
}

# generate data
noisemats <- lapply(1:3, function(i) noise(n, d, eps))
x <- lapply(1:3, function(i) {
  (rbind(data(round((1-jointprop)*n), dirs[i, ]),
      data(n-round((1-jointprop)*n), dirs[4, ])) +
    noisemats[[i]]) / sqrt(n)
})
inds <- matrix(4, 3, 2)
inds[, 1] <- 1:3

joint_dir <- function(V, D) {
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
joint_dir_jive <- function(jivesol) {
  if (jivesol$rankJ == 0 || sum(abs(do.call(rbind, jivesol$joint))) == 0) {
    # return pc1 of residual
    residual <- lapply(1:3, function(i) x[[i]] - jivesol$individual[[i]])
    return(svd(do.call(rbind, residual))$v[, 1])
  }
  svd(do.call(rbind, jivesol$joint))$v[, 1]
}
lambda <- matrix(0, 20, 2)
lambda[, 1] <- exp(seq(-6, 0, length=20))
#mmpcacv <- mmpca::mmpca(x, inds, k=10, enable_sparsity=FALSE, lambda=lambda,
#  parallel=F)$solution
mmpcacv <- mmpca::mmpca(x, inds, k=10, enable_sparsity=FALSE, lambda=lambda,
  parallel=F)
for (i in 1:length(mmpcacv$training)) {
  print(mmpcacv$training[[i]]$iterations)
  print(mmpcacv$training[[i]]$status)
  print(mmpcacv$training[[i]]$message)
}
mmpcacv <- mmpcacv$solution
print('solution all data')
print(mmpcacv$iterations)
print(mmpcacv$status)
print(mmpcacv$message)
cmf <- CMF(x, inds, 10, 0)
cmf4 <- CMF(x, inds, 4, 0)
jive <- r.jive::jive(x, method='perm', showProgress=FALSE, center=FALSE,
  scale=FALSE, orthIndiv=FALSE)
jiveo <- r.jive::jive(x, method='perm', showProgress=FALSE, center=FALSE,
  scale=FALSE, orthIndiv=TRUE)

fh <- gzfile(ifilename)
saveRDS(list(MMPCA=mmpcacv, CMF=cmf, CMF4=cmf4, JIVE=jive, JIVEO=jiveo,
#saveRDS(list(CMF=cmf,
  indata=list(dirs=dirs, noisemats=noisemats, x=x, inds=inds)), fh)
close(fh)
