# find order for rows of each V
mvpca_order <- function(inds, result) {
  n <- max(inds)
  V <- result$solution$V
  R2 <- result$solution$R2_blockwise
  Vorder <- lapply(1:n, function(j) {
    v <- V[[j]]
    # ix = comp importance order for this view
    blocks <- unique(c(which(inds[, 1] == j), which(inds[, 2] == j)))
    R2tmp <- rowSums(R2[, blocks])
    count <- sum(R2tmp > 0)
    ix <- order(-R2tmp)[1:count]
    do.call(order, c(lapply(ix, function(i) {
      vv <- v[, i]
      vv[v[, i] < 0] <- 1
      vv[v[, i] == 0] <- 0
      vv[v[, i] > 0] <- -1
      return(vv)
    }), list(v[, ix[1]])))
  })
}

mvpca_order_subset <- function(inds, result, subset, blocks) {
  n <- max(inds)
  V <- result$solution$V
  R2 <- result$solution$R2_blockwise
  v <- do.call(rbind, V[subset])
  # ix = comp importance order for this view
  #blocks <- unique(c(which(inds[, 1] %in% subset),
    #which(inds[, 2] %in% subset)))
  R2tmp <- rowSums(R2[, blocks])
  count <- sum(R2tmp > 0)
  ix <- order(-R2tmp)[1:count]
  print(ix)
  do.call(order, c(lapply(ix, function(i) {
    vv <- v[, i]
    vv[v[, i] < 0] <- 1
    vv[v[, i] == 0] <- 0
    vv[v[, i] > 0] <- -1
    return(vv)
  }), list(v[, ix[1]])))
}

mvpcaclust <- function(x, inds, result, Vorder) {
  n <- max(inds)
  V <- result$solution$V
  R2 <- result$solution$R2_blockwise
  lapply(1:n, function(j) {
    v <- V[[j]]
    # ix = comp importance order for this view
    blocks <- c(which(inds[, 1] == j), which(inds[, 2] == j))
    R2tmp <- rowSums(R2[, blocks])
    R2tmp <- R2tmp/sum(R2tmp)
    count <- sum(R2tmp > 0)
    ix <- order(-R2tmp)[1:count]
    order <- sapply(ix, function(i) {
      vv <- v[, i]
      vv[v[, i] < 0] <- 1
      vv[v[, i] == 0] <- 0
      vv[v[, i] > 0] <- -1
      return(vv)
    })
    dists <- matrix(NA, nrow(v), nrow(v))
    for (i in 1:nrow(v)) {
      for (k in 1:nrow(v)) {
        d <- sum(R2tmp[ix])
        for (l in 1:count) {
          if (order[i, l] == order[k, l]) {
            d <- sum(R2tmp[ix][-(1:l)])
          } else {
            break
          }
        }
        dists[i, k] <- dists[k, i] <- d
      }
    }
    hc <- hclust(as.dist(dists))
    hc$order <- Vorder[[j]]
    if (j %in% inds[, 1]) {
      lbls <- rownames(x[[which(inds[, 1] == j)[1]]])
    } else {
      lbls <- colnames(x[[which(inds[, 2] == j)[1]]])
    }
    hc$labels <- lbls
    return(hc)
  })
}

mvpcaclust_subset <- function(x, inds, result, Vorder, subset, blocks) {
  n <- max(inds)
  V <- result$solution$V
  R2 <- result$solution$R2_blockwise
  v <- do.call(rbind, V[subset])
  # ix = comp importance order for this view
  #blocks <- unique(c(which(inds[, 1] %in% subset),
    #which(inds[, 2] %in% subset)))
  R2tmp <- rowSums(R2[, blocks])
  R2tmp <- R2tmp/sum(R2tmp)
  count <- sum(R2tmp > 0)
  ix <- order(-R2tmp)[1:count]
  order <- sapply(ix, function(i) {
    vv <- v[, i]
    vv[v[, i] < 0] <- 1
    vv[v[, i] == 0] <- 0
    vv[v[, i] > 0] <- -1
    return(vv)
  })
  dists <- matrix(NA, nrow(v), nrow(v))
  for (i in 1:nrow(v)) {
    for (k in 1:nrow(v)) {
      d <- sum(R2tmp[ix])
      for (l in 1:count) {
        if (order[i, l] == order[k, l]) {
          d <- sum(R2tmp[ix][-(1:l)])
        } else {
          break
        }
      }
      dists[i, k] <- dists[k, i] <- d
    }
  }
  hc <- hclust(as.dist(dists))
  hc$order <- Vorder
  hc$labels <- rownames(v)
  return(hc)
}

mvpca_model_plot <- function(inds, result, Vorder, component_subset=NULL) {
  n <- max(inds)
  V <- result$solution$V
  D <- result$solution$D

  if (is.null(component_subset)) {
    component_subset <- 1:ncol(V[[1]])
  }

  # find position of each block
  ixrow <- list()
  rowcount <- 0
  nrows <- length(unique(inds[, 1]))
  for (i in unique(inds[, 1])) {
    ixrow[[i]] <- 1:nrow(V[[i]])+rowcount
    rowcount <- rowcount + nrow(V[[i]])
  }
  ixcol <- list()
  colcount <- 0
  ncols <- length(unique(inds[, 2]))
  for (i in unique(inds[, 2])) {
    ixcol[[i]] <- 1:nrow(V[[i]])+colcount
    colcount <- colcount + nrow(V[[i]])
  }
  rowborders <- sapply(-sort(-unique(inds[, 1]))[-1],
    function(i) last(ixrow[[i]]))
  if (length(rowborders) == 0) {
    rowborders <- -1
  }
  colborders <- sapply(-sort(-unique(inds[, 2]))[-1],
    function(i) last(ixcol[[i]]))
  if (length(colborders) == 0) {
    colborders <- -1
  }

  blocks1 <- matrix(Inf, rowcount + 3*nrows, colcount + 3*ncols)
  for (i in 1:nrow(inds)) {
    row <- inds[i, 1]
    col <- inds[i, 2]
    x1 <- V[[row]][, component_subset, drop=FALSE] %*%
      diag(D[component_subset, row] * D[component_subset, col],
        length(component_subset)) %*%
      t(V[[col]][, component_subset, drop=FALSE])
    blocks1[ixrow[[row]]+3*row, ixcol[[col]]+3*col] <-
      probit(x1[Vorder[[row]], Vorder[[col]]])
  }

  imm(blocks1)
  #abline(h=0.5+rowborders, v=0.5+colborders)
}

mvpca_data_plot <- function(x, inds, result, Vorder) {
  n <- max(inds)
  V <- result$solution$V
  D <- result$solution$D

  # find position of each block
  ixrow <- list()
  rowcount <- 0
  nrows <- length(unique(inds[, 1]))
  for (i in unique(inds[, 1])) {
    ixrow[[i]] <- 1:nrow(V[[i]])+rowcount
    rowcount <- rowcount + nrow(V[[i]])
  }
  ixcol <- list()
  colcount <- 0
  ncols <- length(unique(inds[, 2]))
  for (i in unique(inds[, 2])) {
    ixcol[[i]] <- 1:nrow(V[[i]])+colcount
    colcount <- colcount + nrow(V[[i]])
  }
  rowborders <- sapply(-sort(-unique(inds[, 1]))[-1],
    function(i) last(ixrow[[i]]))
  if (length(rowborders) == 0) {
    rowborders <- -1
  }
  colborders <- sapply(-sort(-unique(inds[, 2]))[-1],
    function(i) last(ixcol[[i]]))
  if (length(colborders) == 0) {
    colborders <- -1
  }

  blocks <- matrix(Inf, rowcount + 3*nrows, colcount + 3*ncols)
  for (i in 1:nrow(inds)) {
    row <- inds[i, 1]
    col <- inds[i, 2]
    blocks[ixrow[[inds[i, 1]]]+3*row, ixcol[[inds[i, 2]]]+3*col] <-
      probit(x[[i]][Vorder[[inds[i, 1]]], Vorder[[inds[i, 2]]]])
  }

  imm(blocks)
  #abline(h=0.5+rowborders, v=0.5+colborders)
}

imm <- function(x, xlab='', ylab='', asp=TRUE, legend.outside=NA, colors=NULL,
    na.hl=TRUE, ...) {
  if (is.na(legend.outside)) {
    legend.outside <- !asp
  }
  if (asp) {
    asp <- 1
  } else {
    asp <- NA
  }
  d <- dim(x)+1
  if (is.null(colors)) {
    colors <- viridisLite::viridis(256)
  }
  image(1:d[2]-0.5, 1:d[1]-0.5, t(x), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
    col=colors, ...)
  #axis(1, pretty(1:d[2]), pos=d[1]-0.5)
  #axis(2, pretty(1:d[1]), pos=0.5)
  # highlight NA, NaN, Inf
  if (na.hl) {
    ix <- which(is.na(t(x)), arr.ind=TRUE)
    rect(ix[, 1]-0.5, ix[, 2]-0.5, ix[, 1]+0.5, ix[, 2]+0.5, border='red')
  }
  ix <- which(is.infinite(t(x)), arr.ind=TRUE)
  #rect(ix[, 1]-0.5, ix[, 2]-0.5, ix[, 1]+0.5, ix[, 2]+0.5, border='red')
  maxv <- max(x[!is.infinite(x)], na.rm=TRUE)
  minv <- min(x[!is.infinite(x)], na.rm=TRUE)
  delta <- (maxv-minv)/4
  legs <- c(maxv, maxv-delta, maxv-2*delta, maxv-3*delta, minv)
  legcols <- colors[c(256, 192, 128, 64, 1)]
  #if (legend.outside) {
  #  par(mar=c(5.1, 4.1, 4.1, 5.1))
  #  legend('topright', legend=signif(legs, 2), col=legcols, pt.bg=legcols,
  #    pch=22, inset=c(-0.15, 0), xpd=TRUE)
  #} else {
  #  legend('topright', legend=signif(legs, 2), col=legcols, pt.bg=legcols,
  #    pch=22)
  #}
}

probit <- function(x) {
  if (sum(abs(x), na.rm=TRUE) == 0) {
    x[] <- 0.5
    return(x)
  }
  pnorm(x, mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}
probitxy <- function(x, y) {
  if (all(is.na(y))) {
    pnorm(x, mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
  } else {
    pnorm(x, mean(y, na.rm=TRUE), sd(y, na.rm=TRUE))
  }
}
last <- function(x) x[length(x)]

with_rows <- function(x, names) {
  y <- matrix(NA, length(names), ncol(x), dimnames=list(names, colnames(x)))
  is <- intersect(names, rownames(x))
  y[is, ] <- x[is, ]
  return(y)
}

# s: max allowed group size (can be vector for many clusterings with diff group
# sizes)
cutree <- function(tree, s=NULL, ...) {
  if (!is.null(s)) {
    tree <- as.hclust(with_member_height(as.dendrogram(tree)))
    return(cutree(tree, h=s+0.5))
  }
  return(stats::cutree(tree, ...))
}

with_member_height <- function(dg) {
  if (is.leaf(dg)) {
    attr(dg, 'height') <- attr(dg, 'members')
    return(dg)
  }
  dg[[1]] <- with_member_height(dg[[1]])
  dg[[2]] <- with_member_height(dg[[2]])
  attr(dg, 'height') <-
    attr(dg[[1]], 'height') +
    attr(dg[[2]], 'height')
  return(dg)
}

