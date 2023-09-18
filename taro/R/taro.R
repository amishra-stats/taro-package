utils::globalVariables(".")
utils::globalVariables("sample_taxa_table")
# devtools::check(remote = TRUE, manual = TRUE)
# devtools::document()


#'Simulation model
#'
#' genertate simulated data for a tree aggregated factor regression model
#'
#' @param n sample size of the simulated data
#' @param q number of continuous multivariate response
#' @param nrank number of simulated latent factor (rank of the coefficient matrix)
#' @param intercept intercept of the true multivariate model
#' @param snr signal to noise ratio
#' @param rho parameter defining correlated error
#' @param taxa_table input taxa table for the simulated microbiome data
#' @param setting choice of the effect size of the coefficient along the nodes of the tree.
#' @param sparsity desired percentage of the non-zero entries
#' Setting 1: primarily selects rare genera to have significant effect;
#' Setting 2: primarily selects most abundant genera to have significant effect
#' Setting 3: primarily selects leaf nodes to have significant effect
#' Setting 4: primarily selects nodes in the phylogentic tree to have significant effect
#' @return return(list(X=X, Y =Y, A = A, U = U, D=D, V=V, C0 = C0))
#'   \item{Y}{Generated response matrix}
#'   \item{X}{Generated predictor matrix}
#'   \item{A}{Taxonomy inspired adjacency matrix}
#'   \item{tree}{Pass reference phylogenetic tree for visualization}
#'   \item{U}{Left singular vectors of the low-rank coefficient matrix}
#'   \item{V}{Right singular vectors of the low-rank coefficient matrix}
#'   \item{D}{Singular values of the low-rank coefficient matrix}
#'   \item{C0}{Simulated intercept of the coefficient matrix}
#' @export
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom stats rnorm sd
#' @importFrom stats sd
#' @import BiocManager
#' @examples
#' #require(taro)
#'
#' # Simulate data from a sparse factor regression model
#' snr <- .25; xrho <- 0.5; nrank <- 3; q <- 50; n = 300; intercept = 0.5
#' rho = 0.5; sparsity = 0.05; setting = 2;  # error correlation
#' \donttest{
#' sim.sample <- taro_sim(n, q, nrank, rho, snr, intercept,
#'                         taxa_table = NULL, setting, sparsity)
#' Y <- sim.sample$Y
#' X <- sim.sample$X
#' }
#' @references
#' Mishra et al. (2023) \emph{Tree aggregated factor regression model}
taro_sim  = function (n = 300, q = 50, nrank = 3, rho = 0.5,
                      snr = 0.25, intercept = 0.5, taxa_table = NULL,
                      setting = 1, sparsity = 0.05){

  ## test for multiple settings
  # library(taro)
  # snr <- .25; xrho <- 0.5; nrank <- 3; q <- 50; n = 300; intercept = 0.5
  # rho = 0.5;taxa_table = NULL # error correlation
  X2 = X3 = X1 = NULL;
  if (is.null(taxa_table)) {
    # load('taro/data/sample_taxa_table.rda') #
    # utils::data('sample_taxa_table')
    taxa_table <- sample_taxa_table
    nf <- nrow(taxa_table)
  } else {nf <- nrow(taxa_table)}

  # Process the data
  Stool_simulation <- SparseDOSSA2::SparseDOSSA2(template = "Stool",
                                                 n_sample = n,
                                                 n_feature = nf,
                                                 verbose = F)
  otuTab <- Stool_simulation$simulated_data
  rownames(otuTab) <- rownames(taxa_table)
  physeq <- phyloseq::phyloseq( phyloseq::otu_table(otuTab, taxa_are_rows = T),
                                phyloseq::tax_table(as.matrix(taxa_table)) )
  physeq <- phyloseq::subset_taxa(physeq,
                                  rowSums(phyloseq::otu_table(physeq) > 0) > 2 )
  df <- phyloseq::tax_table(physeq) %>% data.frame() %>%
    lapply(factor) %>% as.data.frame()
  formula_str <- paste(setdiff(colnames(df),'unique'),collapse = '/') %>%
    paste('~',.,sep = '')
  phy <- tax_table_to_phylo(eval(parse(text = formula_str)),
                            data = df, collapse = F)
  phy$node.label <- gsub("'",'',phy$node.label)

  A <- phylo_to_A(phy) %>% as.matrix()
  # all(phy$tip.label == rownames(A))
  # colnames(A) <- gsub("'",'',colnames(A))
  ## Find columns to drop
  drop_col <- c()
  for (i in 3:(ncol(df)-1)) {
    temp <- table(df[,c(i-1,i)])
    temp <- rowSums(temp>0) # rowSums(temp)
    drop_col %<>% c(names(temp[temp<2]))
  }
  # dim(A)
  temp <- crossprod(as.matrix(A))
  A <- A[,!(colnames(A) %in% drop_col)]
  X <- phyloseq::otu_table(physeq) %>% t()
  X <- log(X + (X==0))
  tind <- colSums(X) !=0
  X <- X[,tind]; A <- A[tind,]
  A <- A[,colSums(A)>0]
  A <- apply(A, 1, function(x) x/colSums(A)) %>% t()
  Xt <- X %*% A
  # Simulate data from a sparse factor regression model
  p <- ncol(Xt); n <- nrow(Xt)
  U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
  nrandx = function(n){
    c(stats::runif(round(n/2), min = 0.8, max = 1),
      stats::runif(n-round(n/2), min = -1, max = -0.8))
  }
  ref_table <- cbind(colSums(A !=0), 1/as.numeric(apply(Xt, 2, stats::sd)^2),
                     colnames(A)) %>% data.frame() %>% dplyr::arrange(X1,X2)
  if(setting == 2){
    ref_table <- cbind(colSums(A !=0), as.numeric(apply(Xt, 2, sd)^2), colnames(A)) %>%
      data.frame() %>% dplyr::arrange(X1,X2)
  }
  for (i in 1:nrank) {
    ssize <- round(p*sparsity)
    if (ssize%%2 > 0) ssize <- ssize + 1
    M <- cbind(colSums(A),0)
    col_sel <- c(); col_rm <- c()

    if(setting  <= 2){
      ## Setting 1: rare nodes are favoured
      ## Setting 2: most abundant nodes are favoured
      temp <- dplyr::filter(ref_table, X1 == 4)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 1, replace = FALSE, prob = probv/sum(probv)  )
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector()
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 3)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 1, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 2)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 3, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 1) %>% dplyr::filter(!(X3 %in% col_rm) )
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, ssize - 5, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)
    }
    #
    if(setting == 3){
      # Only leaf nodes are selected to have significant effect on the outcome
      temp <- dplyr::filter(ref_table, X1 == 1) %>%
        dplyr::filter(!(X3 %in% col_rm) )
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, ssize, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)
    }
    if(setting == 4){
      # Higher level taxa are favoured to have a significant effect
      temp <- dplyr::filter(ref_table, X1 == 4)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 4, replace = FALSE, prob = probv/sum(probv)  )
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector()
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 3)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 2, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 2)
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, 10, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)

      temp <- dplyr::filter(ref_table, X1 == 1) %>% dplyr:: filter(!(X3 %in% col_rm) )
      probv <- as.numeric(temp$X2)
      tem <- sample(temp$X3, ssize - 16, replace = FALSE, prob = probv/sum(probv))
      col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
        unlist() %>% as.vector() %>% c(.,col_rm)
      col_sel <- c(col_sel,tem)
    }
    #
    M[col_sel,2] <- c(rep(1,ssize/2), rep(-1,ssize/2)) #nrandx(ssize)
    M[M[,2]!=0,2] <- far::orthonormalization(M[M[,2]!=0,], basis=F, norm=F)[,2]
    U[,i] <- M[,2]
    V[5*i + (1:5),i] <- nrandx(q*0.1) # small example, not very sparse
  }
  # U <- far::orthonormalization(cbind(colSums(A),U), basis=F, norm=F)[,2:4]
  U <- apply(U,2,function(x)x/sqrt(sum(x^2)))
  V <- apply(V,2,function(x)x/sqrt(sum(x^2)))
  D <- diag(c(20,15,10))/2
  C <- U%*%D%*%t(V)
  C0 <- matrix(intercept, nrow = 1, ncol = q)
  # Xsigma <- xrho^abs(outer(1:p, 1:p,FUN="-"))
  sim.sample <- simulate_data_mbiome(X,A, U, D, V, n, C0, snr, rho)
  Y <- sim.sample$Y;
  X <- sim.sample$X
  return(list(X=X, Y =Y, A = A, tree = phy, U = U, D=D, V=V, C0 = C0))
}


#' Set of parameters that allows us to control the execution of the TARO algorithms.
#'
#' list of parameters for controlling taro fitting
#'
#' @param mu penalty parameter used in enforcing orthogonality
#' @param nu penalty parameter used in enforcing orthogonality (incremental rate of mu)
#' @param MMerr tolerance in the majorization maximization(MM) algorithm for computing initial values when missing value occurs
#' @param MMiter maximum number iterations in the MM algorithm
#' @param outTol tolerance of convergence of outer loop in CURE
#' @param outMaxIter maximum number of outer loop iteration in CURE
#' @param inMaxIter maximum number of inner loop iteration in CURE
#' @param inTol tolerance value required for convergence of inner loop in CURE
#' @param lamMaxFac a multiplier of calculated lambda_max
#' @param lamMinFac a multiplier of determining lambda_min as a fraction of lambda_max
#' @param gamma0 power parameter in the adaptive weights
#' @param elnetAlpha elastic net penalty parameter
#' @param alpha0 parameter allows us to promote aggregation in the tree guided regression
#' @param se1 set the parameter for promoting 1se rule (se1 = 1) in cross-validation
#' @param spU maximum proportion of nonzero elements in each column of U
#' @param spV maximum proportion of nonzero elements in each column of V
#'
#' @return a list of controlling parameter.
#' @export
#' @examples
#' #require(taro)
#'
#' control <- taro_control()
#'
#' @references
#' Mishra et al. (2023) \emph{Tree aggregated factor regression model}
taro_control = function(mu=1.0, nu=1.1,
                        MMerr=1e-3, MMiter=100,
                        outTol=1e-8, outMaxIter=500,
                        inMaxIter=200, inTol=1e-4,
                        lamMaxFac=1, lamMinFac=1e-10,
                        gamma0=2, elnetAlpha=0.95,
                        alpha0=1,se1 = 1,
                        spU=0.25, spV=0.25) {

  list(mu=mu, nu=nu, MMerr=MMerr, MMiter=MMiter, outTol=outTol, outMaxIter=outMaxIter,
       inMaxIter=inMaxIter, inTol=inTol, lamMaxFac=lamMaxFac, lamMinFac=lamMinFac,
       gamma0=gamma0, alpha0 = alpha0, se1 = se1, elnetAlpha=elnetAlpha, spU=spU ,spV=spV)
}




#' Tree aggregated factor regression model (TARO)
#'
#' A multivariate model for learning the association between continuous multiple
#' outcomes and high-dimensional predictors that are related by a
#' phylogenetic/taxonomic tree
#'
#'
#'
#' @param Y response matrix
#' @param X covariate matrix (Typically log transform of the microbial abundance data)
#' @param A tree guided adjacency matrix
#' @param Ac matrix (LHS) for imposing the linearity constraint on the parameters
#' of the coefficient matrix
#' @param Bc vector (RHS) for imposing the linearity constraint on the parameters
#' of the coefficient matrix
#' @param Z control variables
#' @param phylo reference phylogenetic tree
#' @param maxrank an integer specifying the maximum desired rank/number of factors
#' @param nlambda maximum number of lambda values used for generating the solution path
#' @param orthV if TRUE, orthogonality of the left singular vector is imposed or not
#' @param nfold number of folds for the k-fold cross validation for selecting the tuning parameter
#' @param verbose if TRUE, prints the message implementing the sequential procedure
#' @param control a list of internal parameters controlling the model fitting
#' @aliases taro
#' fit = fit.nlayer, lam = lamSel
#' @return
#'   \item{C}{estimated coefficient matrix}
#'   \item{Z}{estimated parameter corresponding to control variables}
#'   \item{U}{estimated left singular vector matrix }
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated right singular vector matrix (loading matrix)}
#'   \item{lam}{selected lambda values based on the chosen information criterion}
#'   \item{fit}{optimal parameters estimate including crossvalidation error estimate from each of the sequential step}
#' @export
#' @import magrittr
#' @importFrom  Rcpp evalCpp
#' @importFrom igraph set_vertex_attr
#' @importFrom ape as.phylo
#' @useDynLib taro
#' @examples
#' require(taro)
#' set.seed(123)
#'
#' #' # Simulate data from a sparse factor regression model
#' snr <- .25; xrho <- 0.5; nrank <- 3; q <- 50; n = 300; intercept = 0.5
#' rho = 0.5 # error correlation
#' \donttest{
#' input_data <- taro_sim(n, q, nrank, rho, snr, intercept,
#'                         taxa_table = NULL)
#'
#' Y <- input_data$Y
#' X <- input_data$X
#' A <- input_data$A
#' phylo <- input_data$tree
#' n <- nrow(Y); q <- ncol(Y)
#' set.seed(123)
#' # Implementation setting
#' maxrank = 5;
#' Z = NULL;
#' A = A;
#' Ac <- matrix(1,1,ncol = ncol(X)) %*% A
#' Bc <- matrix(0,1,1)
#' nfold = 5; trace = TRUE; verbose = TRUE
#' nlambda = 100;PATH = TRUE
#' control <- taro_control(alpha0 = 0.5, gamma0 = 1, spU = 0.5,
#'                        inTol = 1e-5, inMaxIter = 300,
#'                        outMaxIter = 1000,outTol = 1e-8,
#'                        spV=0.5, lamMaxFac = 1e2, se1 = 1)
#' # Weight Yes
#' # fit_seq <- taro_path(Y, X, A, Ac, Bc, Z = NULL, phylo = phylo,
#' #                       maxrank = 5, nlambda = 100,
#' #                       control = control,
#' #                       nfold = 5, orthV = TRUE,
#' #                       verbose = TRUE)
#' }
#'@references
#' Mishra et al. (2023) \emph{Tree aggregated factor regression model}
taro_path <- function(Y, X, A, Ac, Bc, Z = NULL, phylo = NULL,
                       maxrank = 10, nlambda = 40,
                       control = list(),
                       nfold = 5, orthV = TRUE,
                       verbose = TRUE) {
  if(verbose) cat("Initializing...", "\n")
  n <- nrow(Y)
  p <- ncol(A)
  q <- ncol(Y)
  if (ncol(Ac) != p) stop('incorrect dimension of Ac')
  if (nrow(Bc) != nrow(Ac)) stop('incorrect dimension of linear constraint')
  taxTree <- A_to_igraph(A)

  if (is.null(Z)) {
    cIndex <- 1
    Xt <- X %*% A
    X0 <- cbind(1, Xt)
  } else {
    cIndex <- c(1, (1:ncol(Z)) + 1)
    Xt <- X %*% A
    X0 <- cbind(1, Z, Xt)
  }
  wu <- as.matrix(apply(Xt, 2, stats::sd))

  Yk <- Y
  naind <- (!is.na(Y)) + 0
  misind <- any(naind == 0) + 0
  if (misind == 1) Yk[is.na(Y)] <- 0
  Z0 <- stats::coef(stats::lm(Yk~X0[,cIndex,drop = FALSE]-1))

  # Output parameters matrix
  U <- matrix(0, nrow = p, ncol = maxrank)
  V <- matrix(0, nrow = q, ncol = maxrank)
  D <- rep(0, maxrank)
  Z <- matrix(0, nrow = length(cIndex), ncol = q)
  lamSel <- rep(0, maxrank)

  ## Begin exclusive extraction procedure
  N.nna <- sum(!is.na(Y))
  naind <- (!is.na(Y)) + 0
  Yk <- Y
  ind.nna <- which(!is.na(Y))
  misind <- any(naind == 0) + 0
  ID <- sample(rep(1:nfold, len = N.nna), N.nna, replace = FALSE)
  fit.nlayer <- fit.nfold <- vector("list", maxrank)
  xx_in <- list(u = matrix(rep(1,p)), d = 1, v = matrix(rep(1,q)) )
  ortU <- t(rep(0, p)); Bu <- 0
  ortV <- t(rep(0, q)); Bv <- 0
  for (k in 1:maxrank) {  #  k=4
    message("Estimation of unit-rank unit: ", k, "\n")
    Yini <- Yk;  Yini[is.na(Yini)] <- 0
    Z0 <- stats::coef(stats::lm(Yini~X0[,cIndex,drop = FALSE]-1))
    lamMat <- crossprod(Xt,Yini) - crossprod(Xt,X0[,cIndex,drop = FALSE]%*% Z0)
    xx_in <- tarorrr_cpp_ur2(Yini, X0, cIndex,
                            Ac,  Bc,
                            Z0, matrix(0,p,q),
                            control)
    message(xx_in$d, "\n")
    initW <- list(wu = wu*(colSums(A!=0)^control$alpha0) ,
                  wd = abs(xx_in$d)^(-control$gamma0),
                  wv = abs(xx_in$v)^(-control$gamma0))
    lamMat <- lamMat/((initW$wd*initW$wu) %*% t(initW$wv))
    message("Cross validation of unit-rank component:", k, "\n")
    dev <- matrix(NA, nfold, nlambda)
    fitT <- vector("list", nfold)
    for (ifold in 1:nfold) { # ifold = 1
      message("Fold ", ifold, "\n")
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Yk
      Yte[-ind.test] <- NA
      Ytr <- Yk
      Ytr[ind.test] <- NA
      sste <- sum(!is.na(Yte))

      naindC <- (!is.na(Ytr)) + 0
      misindC <- any(naindC == 0) + 0
      Ytr[is.na(Ytr)] <- 0
      # xx_in <- list(u = matrix(rep(1,p)), d = 1, v = matrix(rep(1,q)) )

      # Rcpp::sourceCpp("taro/src/taro.cpp")
      fitT[[ifold]] <- taro_cpp_ur(Ytr,X0,Ac,Bc, nlambda, cIndex,
                                   initW, Dini = xx_in$d, Zini = Z0,
                                   Uini = xx_in$u, Vini = xx_in$v,
                                   ortV, Bv,
                                   control, misindC, 1, naindC, lamMat)

      for (im in 1:fitT[[ifold]]$nkpath) dev[ifold, im] <-
        sum((Yte - fitT[[ifold]]$etapath[,,im])^2,na.rm = T)/sste

    }
    dev.mean <- colMeans(dev, na.rm = FALSE)
    sderr <- apply(dev,2,stats::sd,na.rm = F)/sqrt(nfold)
    l.mean <- which.min(dev.mean)
    if (is.na(sderr[l.mean])) {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean])))
    } else {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean] +
                                         control$se1 * sderr[l.mean])))
    }
    lamS <- fitT[[ifold]]$lamseq[l.mean]
    # save(list = ls(), file = 'Aditya.rda')
    fit.nlayer[[k]] <- fit.layer <- taro_cpp_ur(Yini,X0,Ac,Bc, nlambda, cIndex,
                                                initW, Dini = xx_in$d, Zini = Z0,
                                                Uini = xx_in$u, Vini = xx_in$v,
                                                ortV, Bv,
                                                control, misind, 1, naind, lamMat)
    l.mean <- which(fit.layer$lamKpath == lamS)
    if (length(l.mean) == 0)  l.mean <- 1

    fit.nlayer[[k]]$initw <- initW
    fit.nlayer[[k]]$dev <- dev
    fit.nlayer[[k]]$lamS <- lamS
    fit.nlayer[[k]]$sderr <- sderr
    U[, k] <- fit.layer$ukpath[, l.mean]
    V[, k] <- fit.layer$vkpath[, l.mean]
    D[k] <- fit.layer$dkpath[l.mean, 1]
    message(D[k], "\n")
    lamSel[k] <- fit.layer$lamKpath[l.mean, 1]
    if (D[k] == 0) {
      if (k > 1) {
        U <- matrix(U[, 1:(k - 1)], ncol = k - 1)
        V <- matrix(V[, 1:(k - 1)], ncol = k - 1)
        D <- D[1:(k - 1)]
        lamSel <- lamSel[1:(k - 1)]
      }
      break
    }
    Ck <- D[k] * tcrossprod(U[, k], V[, k])
    Z <- fit.layer$zpath[, , l.mean, drop=F]
    Yk <- Yk - X0[, -cIndex] %*% Ck
    ortU <- t(rep(0, p)); Bu <- 0
    ortV <- t(V[,1:k,drop=F]); Bv <- rep(0,k)
    if(!orthV) ortV <- 0*ortV
  }
  message("Estimated rank =", sum(D != 0), "\n")
  if (sum(D != 0) == maxrank) {
    message("Increase maxrank value!")
  }
  # taxTree <- A_to_igraph(A)
  # all(igraph::V(taxTree)$name == c(colnames(A), 'Root'))
  # for (i in 1:ncol(U))
  #   taxTree %<>% igraph::set_vertex_attr(sprintf("u%s",i),
  #                                        value = c(U[,i],0) )

  if (is.null(phylo)) phylo <- ape::as.phylo(A_to_igraph(A))
  ft1 <- list(fit = fit.nlayer, C = U %*% (D * t(V)), Z = Z,
              U = U, V = V, D = D, lam = lamSel,tree = phylo)
  return(ft1)
}









#' Linear constrained reduced rank regression (rrr_lc)
#'
#' Compute the parameters estimate of a multivariate linear constrained reduced
#' rank regression model. The computational procedure estimate the rank using
#' cross-validation procedure.
#'
#' @param Y response matrix
#' @param X covariate matrix (Typically log transform of the microbial abundance data)
#' @param A tree guided adjacency matrix
#' @param Ac matrix (LHS) for imposing the linearity constraint on the parameters
#' of the coefficient matrix
#' @param Bc vector (RHS) for imposing the linearity constraint on the parameters
#' of the coefficient matrix
#' @param Z control variables
#' @param phylo reference phylogenetic tree
#' @param maxrank an integer specifying the maximum desired rank/number of factors
#' @param nfold number of folds for the k-fold cross validation for selecting the tuning parameter
#' @param verbose if TRUE, prints the message implementing the sequential procedure
#' @param control a list of internal parameters controlling the model fitting
#' @return
#'   \item{C}{estimated coefficient matrix}
#'   \item{Z}{estimated parameter corresponding to control variables}
#'   \item{U}{estimated left singular vector matrix }
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated right singular vector matrix (loading matrix)}
#' @export
#' @import magrittr
#' @importFrom graphics title
#' @importFrom stats lm
#' @importFrom stats coef
#' @useDynLib taro
#' @examples
#' require(taro)
#'
#' # Simulate data from a sparse factor regression model
#' snr <- .25; xrho <- 0.5; nrank <- 3; q <- 50; n = 300; intercept = 0.5
#' rho = 0.5 # error correlation
#' \donttest{
#' input_data <- taro_sim(n, q, nrank, rho, snr, intercept,
#'                         taxa_table = NULL)
#'
#' Y <- input_data$Y
#' X <- input_data$X
#' A <- input_data$A
#' phylo <- input_data$phylo
#' n <- nrow(Y); q <- ncol(Y)
#' set.seed(123)
#' # Implementation setting
#' maxrank = 5;
#' Z = NULL;
#' A = A;
#' Ac <- matrix(1,1,ncol = ncol(X)) %*% A
#' Bc <- matrix(0,1,1)
#' nfold = 5; trace = TRUE; verbose = TRUE
#' control <- taro_control(alpha0 = 0.5, gamma0 = 1, spU = 0.5,
#'                        inTol = 1e-5, inMaxIter = 300,
#'                        outMaxIter = 1000,outTol = 1e-8,
#'                        spV=0.5, lamMaxFac = 1e2, se1 = 1)
#'
#' # test_srrr <- rrr_lc(Y, X, A, Ac, Bc, Z = NULL, maxrank = 5,
#' # control, nfold = 5)
#' }
#' @references
#' Mishra et al. (2023) \emph{Tree aggregated factor regression model}
rrr_lc = function(Y, X, A, Ac, Bc, Z = NULL, phylo = NULL,
                  maxrank = 10, control = list(),
                  nfold = 5, verbose = TRUE) {
  rrr_lc_miss = function(Yr, Xr, k, cIndex,
                         Ac, Bc,
                         Zini, control,naind){
    out <- tarorrr(Yr, Xr, k, cIndex,
                   Ac, Bc, Zini, control)
    if(any(naind==0)){
      for (j in 1:(control$MMiter)) {
        C_temp <- out$C
        pred_val <- Xr %*% out$C
        miss_ind <- naind == 0
        Yr[miss_ind] <- pred_val[miss_ind]
        out <- tarorrr(Yr, Xr, k, cIndex,
                       Ac, Bc,  Zini, control)
        err <- norm(C_temp-out$C,'F')/norm(C_temp,'F')
        if (err < 100*control$MMerr) break;
      }
    }
    return(out)
  }

  if(verbose) cat("Initializing...", "\n")
  n <- nrow(Y)
  p <- ncol(A)
  q <- ncol(Y)
  if (ncol(Ac) != p) stop('incorrect dimension of A')
  if (nrow(Bc) != nrow(Ac)) stop('incorrect dimension of B')

  if (is.null(Z)) {
    cIndex <- 1
    X0 <- cbind(1, X %*% A)
  } else {
    cIndex <- c(1, (1:ncol(Z)) + 1)
    X0 <- cbind(1, Z, X %*% A)
  }

  Yin <- Y
  naind <- (!is.na(Y)) + 0
  misind <- any(naind == 0) + 0
  if (misind == 1) Yin[is.na(Y)] <- 0
  Z0 <- stats::coef(stats::lm(Yin~X0[,cIndex,drop = FALSE]-1))

  ## Begin exclusive extraction procedure
  N.nna <- sum(!is.na(Y))
  naind <- (!is.na(Y)) + 0
  ind.nna <- which(!is.na(Y))
  fit.nfold <- vector("list", maxrank)
  # generate kfold index
  ID <- rep(1:nfold, len = N.nna)
  ID <- sample(ID, N.nna, replace = FALSE)
  ## store the deviance of the test data
  dev <- matrix(NA, nfold, maxrank)
  dev0 <- rep(NA,nfold)
  # rank selection via cross validation
  for (k in 1:maxrank) { # desired rank extraction
    if(verbose) cat('Cross-validation for rank r: ', k, '\n')
    fitT <- vector("list", nfold)
    for (ifold in 1:nfold) { # ifold=3; k =1
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Y
      Yte[-ind.test] <- NA
      Ytr <- Y
      Ytr[ind.test] <- NA
      sste <- sum(!is.na(Yte))

      naind2 <- (!is.na(Ytr)) + 0
      misind <- any(naind2 == 0) + 0
      Ytr[is.na(Ytr)] <- 0

      fitT[[ifold]] <- rrr_lc_miss(Ytr, X0, k, cIndex,
                                   Ac, Bc, Z0, control, naind2)

      cat('Fold ', ifold, ': [Error,iteration] = [',
          with(fitT[[ifold]],diffobj[maxit]),
          with(fitT[[ifold]],maxit), ']', '\n')

      # compute test sample error
      if(k==1) dev0[ifold]  <- sum((Yte - X0[,cIndex,drop = F] %*% Z0)^2,
                                   na.rm = T)/sste
      dev[ifold, k] <- sum((Yte - X0 %*% fitT[[ifold]]$C)^2, na.rm = T)/sste;
    }
    fit.nfold[[k]] <- fitT
  }
  # select rank with lowest mean;
  dev.mean <- colMeans(cbind(dev0,dev), na.rm = T)
  rank_sel <- which.min(dev.mean)-1

  if (is.null(phylo)) phylo <- ape::as.phylo(A_to_igraph(A))

  # compute model estimate for the selected rank
  Y[is.na(Y)] <- 0
  if (rank_sel > 0){
    out <- rrr_lc_miss(Y, X0, rank_sel, cIndex,
                       Ac, Bc, Z0, control, naind)

    out$cv.err <- cbind(dev0,dev)
    out$C <- out$C[-cIndex, ,drop=F]
    out$Z <- out$C[cIndex, ,drop=F]
    out$V <- out$v
    out$D <- out$d[,1]
    out$U <- out$u
    out$Y <- Y;  out$X = X
  } else {
    C <- matrix(rep(0,ncol(X0)*q), ncol(X0),q)
    C[cIndex,] <- Z0
    out = list(C = C, diffobj =NA, converged = NA, ExecTimekpath = NA,
               maxit = NA, converge = NA, Z = Z0, D = 0,
               U = matrix(rep(0,ncol(X0) -length(cIndex))),
               V = matrix(rep(0,q)), Y = Y, X = X, tree = phylo)
  }
  return(out)
}






#
# taro_sim  = function (n = 300, q = 50, nrank = 3, rho = 0.5,
#                       snr = 0.25, intercept = 0.5, taxa_table = NULL){
#   X2 = X3 = X1 = NULL;
#   if (is.null(taxa_table)) {
#     # load('taro/data/sample_taxa_table.rda') #
#     # utils::data('sample_taxa_table')
#     taxa_table <- sample_taxa_table
#     nf <- nrow(taxa_table)
#   } else {nf <- nrow(taxa_table)}
#
#   # Process the data
#   Stool_simulation <- SparseDOSSA2::SparseDOSSA2(template = "Stool",
#                                                  n_sample = n,
#                                                  n_feature = nf,
#                                                  verbose = F)
#   otuTab <- Stool_simulation$simulated_data
#   rownames(otuTab) <- rownames(taxa_table)
#   physeq <- phyloseq::phyloseq( phyloseq::otu_table(otuTab, taxa_are_rows = T),
#                                 phyloseq::tax_table(as.matrix(taxa_table)) )
#   physeq <- phyloseq::subset_taxa(physeq,
#                                   rowSums(phyloseq::otu_table(physeq) > 0) > 2 )
#   df <- phyloseq::tax_table(physeq) %>% data.frame() %>%
#     lapply(factor) %>% as.data.frame()
#   formula_str <- paste(setdiff(colnames(df),'unique'),collapse = '/') %>%
#     paste('~',.,sep = '')
#   phy <- tax_table_to_phylo(eval(parse(text = formula_str)),
#                             data = df, collapse = F)
#   A <- phylo_to_A(phy) %>% as.matrix()
#   # all(phy$tip.label == rownames(A))
#   colnames(A) <- gsub("'",'',colnames(A))
#   ## Find columns to drop
#   drop_col <- c()
#   for (i in 3:(ncol(df)-1)) {
#     temp <- table(df[,c(i-1,i)])
#     temp <- rowSums(temp>0) # rowSums(temp)
#     drop_col %<>% c(names(temp[temp<2]))
#   }
#   # dim(A)
#   temp <- crossprod(as.matrix(A))
#   A <- A[,!(colnames(A) %in% drop_col)]
#   X <- phyloseq::otu_table(physeq) %>% t()
#   X <- log(X + (X==0))
#   tind <- colSums(X) !=0
#   X <- X[,tind]; A <- A[tind,]
#   A <- A[,colSums(A)>0]
#   A <- apply(A, 1, function(x) x/colSums(A)) %>% t()
#   Xt <- X %*% A
#   # Simulate data from a sparse factor regression model
#   p <- ncol(Xt); n <- nrow(Xt)
#   U <- matrix(0,ncol=nrank ,nrow=p);  V <- matrix(0,ncol=nrank ,nrow=q)
#   nrandx = function(n){
#     c(stats::runif(round(n/2), min = 0.8, max = 1),
#       stats::runif(n-round(n/2), min = -1, max = -0.8))
#   }
#   ref_table <- cbind(colSums(A !=0), 1/as.numeric(apply(Xt, 2, stats::sd)^2),
#                      colnames(A)) %>% data.frame() %>% dplyr::arrange(X1,X2)
#   for (i in 1:nrank) {
#     ssize <- round(p*0.05)
#     if (ssize%%2 > 0) ssize <- ssize + 1
#     M <- cbind(colSums(A),0)
#     col_sel <- c(); col_rm <- c()
#     #
#     temp <- dplyr::filter(ref_table, X1 == 1) %>%
#       dplyr::filter(!(X3 %in% col_rm) )
#     probv <- as.numeric(temp$X2)
#     tem <- sample(temp$X3, ssize, replace = FALSE, prob = probv/sum(probv))
#     col_rm <- apply(A[,tem, drop=F],2, function(x) colnames(A)[x != 0]) %>%
#       unlist() %>% as.vector() %>% c(.,col_rm)
#     col_sel <- c(col_sel,tem)
#     #
#     M[col_sel,2] <- c(rep(1,ssize/2), rep(-1,ssize/2)) #nrandx(ssize)
#     M[M[,2]!=0,2] <- far::orthonormalization(M[M[,2]!=0,], basis=F, norm=F)[,2]
#     U[,i] <- M[,2]
#     V[5*i + (1:5),i] <- nrandx(q*0.1) # small example, not very sparse
#   }
#   # U <- far::orthonormalization(cbind(colSums(A),U), basis=F, norm=F)[,2:4]
#   U <- apply(U,2,function(x)x/sqrt(sum(x^2)))
#   V <- apply(V,2,function(x)x/sqrt(sum(x^2)))
#   D <- diag(c(20,15,10))/2
#   C <- U%*%D%*%t(V)
#   C0 <- matrix(intercept, nrow = 1, ncol = q)
#   # Xsigma <- xrho^abs(outer(1:p, 1:p,FUN="-"))
#   sim.sample <- simulate_data_mbiome(X,A, U, D, V, n, C0, snr, rho)
#   Y <- sim.sample$Y;
#   X <- sim.sample$X
#   return(list(X=X, Y =Y, A = A, U = U, D=D, V=V, C0 = C0))
# }

