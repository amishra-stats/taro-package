#' Convert a tax table to a phylo object
#'
#' This is mostly the function \code{\link[ape]{as.phylo.formula}} from the ape package
#' but modified so that interior nodes get labeled.
#' Note: assumes that there are NOT single quotes around the labels.
#' @inheritParams ape::as.phylo.formula
#' @author Jacob Bien \email{jbien@usc.edu}
#' @export
tax_table_to_phylo <- function (x, data = parent.frame(), collapse = TRUE, ...) {
  err <- "Formula must be of the kind ~A1/A2/.../An."
  if (any(lapply(data, class) != "factor"))
    stop("Every column of data must be a factor.")
  if (length(x) != 2)
    stop(err)
  if (x[[1]] != "~")
    stop(err)
  f <- x[[2]]
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/")
      stop(err)
    f3.txt <- deparse(f[[3]])
    if (!is.factor(data[[f3.txt]]))
      stop(paste("Variable", f3.txt, "must be a factor"))
    taxo[[f3.txt]] <- data[[f3.txt]]
    if (length(f) > 1)
      f <- f[[2]]
  }
  f.txt <- deparse(f)
  if (!is.factor(data[[f.txt]]))
    stop(paste("Variable", f.txt, "must be a factor."))
  taxo[[f.txt]] <- data[[f.txt]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo))
        warning("leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u-1)])
      t[l] <- paste0("(", paste(x, collapse = ","), ")", "'", levels[l], "'")
    }
    t
  }
  string <- paste0(f.rec(taxo.data), ";")
  #string <- paste0("(", paste(f.rec(taxo.data), collapse = ","), ");")
  phy <- ape::read.tree(text = string)
  if (collapse)
    phy <- ape::collapse.singles(phy)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  phy
}

#' Convert from phylo to the A matrix used in trac.  Note this is similar to
#' the A used in rare, but with the column of all ones (for the root) removed
#' @param phy \code{phylo} object from the \code{ape} package
#' @author Jacob Bien \email{jbien@usc.edu}
#' @export
phylo_to_A <- function(phy) {
  nleaves <- length(phy$tip.label)
  g <- igraph::graph_from_edgelist(phy$edge,
                                   directed = TRUE)
  igraph::V(g)$name <- c(phy$tip.label, phy$node.label)
  num_nodes <- length(igraph::V(g))
  A <- cbind(Matrix::Diagonal(nleaves),
             Matrix::Matrix(0, nleaves, num_nodes - nleaves))
  for (i in seq(nleaves + 1, num_nodes)) {
    # for all nonleaf nodes
    leaves_of_i <- intersect(igraph::subcomponent(g, mode = "out", v = i),
                             1:nleaves)
    A[leaves_of_i, i] <- 1
  }
  A[, seq(nleaves + 1, num_nodes)] <- A[, rev(seq(nleaves + 1, num_nodes))]
  rownames(A) <- phy$tip.label
  colnames(A) <- c(phy$tip.label, rev(phy$node.label))
  A[, -ncol(A)] # remove the root, which is the rightmost column
}



#' Simulate multivariate response matrix using covariates, adjacency matrix,
#' svd decomposition of the coefficient matrix (U,D,V), sample size,
#' intercept, signal to noise ratio and covariance matrix parameter of multivariate
#' error.
#'
#' @param X log transformed microbial abundance data
#' @param A phylogeny inspired adjacency matrix
#' @param U left singular vector of the svd decomposition of the coefficient matrix
#' @param V right singular vector of the svd decomposition of the coefficient matrix
#' @param D singular values of the svd decomposition of the coefficient matrix
#' @param n sample size
#' @param C0 intercept value
#' @param snr signal-to-noise ratio
#' @param rho error matrix correlation
#' @author Aditya Mishra \email{akmishra@mdanderson.org}
simulate_data_mbiome  = function (X,A, U, D, V, n, C0, snr, rho = 0.5) {
  q <- nrow(V)
  svdrr <- eigen(rho^abs(outer(1:q, 1:q, FUN = "-")))
  svdrrinv <- svdrr$vectors %*% diag(svdrr$values^0.5, nrow = q) %*%
    t(svdrr$vectors)
  UU <- matrix(nrow = n, ncol = q, stats::rnorm(n * q, 0, 1)) %*%
    svdrrinv
  nrank <- ncol(U)
  C <- U %*% D %*% t(V)
  Y <- cbind(1,X %*% A) %*% rbind(C0,C)
  Y3 <- X %*% A %*% U[, nrank,drop=F] %*% t(V[, nrank,drop=F]) * D[nrank, nrank]
  sigma <- sqrt(sum(Y3^2)/sum(UU^2))/snr
  UU <- UU * sigma
  Y <- Y + UU
  return(list(Y = Y, X = X, sigma =sigma))
}


