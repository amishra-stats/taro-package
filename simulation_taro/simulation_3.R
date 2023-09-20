# Setting 3:Only leaf nodes are selected to have significant effect on the outcome

library(taro)
set.seed(123)

nsim <- 50
snr <- .25; xrho <- 0.5; nrank <- 3;
q <- 50; n = 300; intercept = 0.5
rho = 0.5; setting = 3; sparsity = 0.05
simulation_3 <- vector('list', nsim)
ii = 1; se_id <- 1
while (ii <= nsim) {
  tryCatch({
    set.seed(se_id*1234)
    se_id <- se_id + 1
    sim.sample <- taro_sim(n, q, nrank, rho, snr, intercept,
                           taxa_table = NULL, setting, sparsity)
    if(ncol(sim.sample$X) > 200){
      simulation_3[[ii]] <- sim.sample
      ii <- ii + 1
    }
  },  error=function(error_message) {
    message("This is my custom message.")
    message(error_message)
    return(NA)
  }
  )
}
save(simulation_3, file = 'data/simulation3.RData')






## Running on cluster for estimation
rm(list = ls())
load('data/simulation3.RData')
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
library(magrittr)
library(Matrix)
library(phyloseq)
library(trac)
library(reticulate)
# upload python path and install classo package in python
use_python("~/anaconda3/envs/akmishra/bin/python")
classo = import("classo")

trac_marginal = function(X0,Y0,A_trac){
  C_trac <- matrix(0,ncol(A_trac),ncol(Y0) )
  for (i in 1:ncol(Y0)) {
    naind <- !is.na(Y0[,i])
    tem_trac <- trac(X0[naind,], Y0[naind,i], A_trac)
    cv_fit <- cv_trac(tem_trac, X0[naind,], Y0[naind,i], A_trac)
    C_trac[,i] <- tem_trac[[1]]$alpha[, cv_fit$cv[[1]]$ibest]
  }
  # C_trac <- fit_trac
  svdc <- svd(C_trac)
  D <- svdc$d
  V <- svdc$v[,abs(D) > 1e-5]
  U <- svdc$u[,abs(D) > 1e-5]
  D <- D[abs(D) > 1e-5]
  list(C = C_trac, U = U, D = D, V=V, Z = matrix(0.5,1,ncol(Y0)) )
}


index <- 1
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

input_data <- simulation_3[[index]]

miss <- 0.05
Y <- input_data$Y
X <- input_data$X
A <- input_data$A
n <- nrow(Y); q <- ncol(Y)
t.ind <- sample.int(n*q, size = miss*n*q)
y <- as.vector(Y); y[t.ind] <- NA;  Ym <- matrix(y,n,q)
set.seed(123)


##
library(gofar)

family <- list(gaussian(), binomial(), poisson())
control1 <- gofar_control()
nlam <- 40 # number of tuning parameter
SD <- 123
# Model fitting begins:
control1$epsilon <- 1e-7
control1$spU <- 0.5
control1$spV <- 0.5
control1$maxit <- 1000
familygroup <- rep(1,q)

set.seed(SD)
rank.est <- 5
fit.seq <- gofar_s(Y, X %*% A,
                   nrank = rank.est, family = family,
                   nlambda = nlam, familygroup = familygroup,
                   control = control1, nfold = 5)


set.seed(SD)
rank.est <- 5
fit.seqM <- gofar_s(Ym, X %*% A,
                    nrank = rank.est, family = family,
                    nlambda = nlam, familygroup = familygroup,
                    control = control1, nfold = 5)




# Marginal trac
A_trac <- (A!=0)+0
fit_trac <- trac_marginal(X,Y-0.5,A_trac)


# Implementation setting
maxrank = 5;
Z = NULL;
A = A;
Ac <- matrix(1,1,ncol = ncol(X)) %*% A
Bc <- matrix(0,1,1)
nfold = 5; trace = T; verbose = TRUE
nlambda = 100;PATH = T

control <- taro_control(alpha0 = 0.5, gamma0 = 1, spU = 0.5,
                        inTol = 1e-5, inMaxIter = 300,
                        outMaxIter = 1000,outTol = 1e-8,
                        spV=0.5, lamMaxFac = 1e2, se1 = 1)

# Weight Yes
fit_seqT <- taro_path(Y, X, A, Ac, Bc, Z = NULL,
                      maxrank = 5, nlambda = 100,
                      control = control,
                      nfold = 5, orthV = TRUE,
                      verbose = TRUE)


# Weight Yes with missing data
fit_seqTM <- taro_path(Ym, X, A, Ac, Bc, Z = NULL,
                       maxrank = 5, nlambda = 100,
                       control = control,
                       nfold = 5, orthV = TRUE,
                       verbose = TRUE)

# Marginal trac with missing data
A_trac <- (A!=0)+0
fit_tracM <- trac_marginal(X,Ym-0.5,A_trac)


# factor model as SRRR
test_srrr <- rrr_lc(Y, X, A, Ac, Bc,
                    Z = NULL, maxrank = 5,
                    control, nfold = 5)

# factor model as SRRR  with missing data
test_srrrM <- rrr_lc(Ym, X, A, Ac, Bc,
                     Z = NULL, maxrank = 5,
                     control, nfold = 5)


control <- taro_control(alpha0 = 0.5, gamma0 = 0, spU = 0.5,
                        inTol = 1e-5, inMaxIter = 300,
                        outMaxIter = 1000,outTol = 1e-8,
                        spV=0.5, lamMaxFac = 1e2, se1 = 1)

# Weight No
fit_seqF <- taro_path(Y, X, A, Ac, Bc, Z = NULL,
                      maxrank = 5, nlambda = 100,
                      control = control,
                      nfold = 5, orthV = TRUE,
                      verbose = TRUE)


# Weight No with missing data
fit_seqFM <- taro_path(Ym, X, A, Ac, Bc, Z = NULL,
                       maxrank = 5, nlambda = 100,
                       control = control,
                       nfold = 5, orthV = TRUE,
                       verbose = TRUE)

save(list = ls(), file = sprintf('output/simulation3_%s.RData',index))

# Code generating script to run simulation on server
# sink("simulation3.sh")
# for(i in 1:50){
#   cat(sprintf('bsub -W 3:00 -q short -n 1 -M 16 -R \'rusage[mem=6]\' -o out/out3_%s -e err/err3_%s \"Rscript simulation3.R %s 1>Rout/out3_%s 2>Rerr/err3_%s\"\n',i,i,i,i,i))
# }
# sink()









## --------------------------------------------------------

## Model comparisons
rm(list = ls())
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
# load('data/processed_data.rds')
library(magrittr)
library(Matrix)
library(phyloseq)
library(tidyverse)


## A function to compute the error measures for comparing different methods
## Error in C, Error in XC, FPR, FNR,  R%, rank measure
compare_taro  = function(U,D,V,C0,n,X, A, Z=NULL, fit){
  Xt <- X %*% A
  if (is.null(Z)) {
    cVar <- matrix(1,n,1)
  } else {
    cVar <- cbind(1, Z)
  }
  p <- nrow(U); q <- nrow(V)
  C <- U%*%D%*%t(V)
  nrank <- ncol(U)
  Xsigma <- crossprod(cbind(cVar,Xt))
  cy <- rbind(C0,C)-rbind(fit$Z,fit$C)
  Err.C <- 100*(norm(cy,"f"))^2 /(p*q)                ## Er(C)
  Err.Y <- 100*sum(diag(t(cy)%*%Xsigma%*%cy))/(n*q)   ## Er(XC)
  ord <- order(fit$D, decreasing = T)
  ord <- apply(abs(cor(fit$V,V)), 2, which.max)
  ord <- c(ord, setdiff(1:ncol(fit$V),ord))
  U.est <- fit$U[,ord]; V.est <- fit$V[,ord]; D.est <- fit$D[ord]
  D1.est <- diag(fit$D[ord])
  beta.sim <- c(as.vector(U),as.vector(V))
  truezero <- sum(beta.sim==0)
  truenonzero <- sum(beta.sim!=0)
  rank.est <- NA

  if(is.null(U.est) | is.null(V.est) ) {
    fp <- truezero;fn <- 0
  } else {
    Ue <- matrix(0,p,nrank); Ve <- matrix(0,q,nrank)
    rank.est <- ncol(U.est)
    if(nrank < rank.est){
      Ue <- U.est[,1:nrank]; Ve <- V.est[,1:nrank]
    } else{
      Ue[,1:rank.est] <- U.est; Ve[,1:rank.est] <- V.est
    }
    beta.est <- c(as.vector(Ue),as.vector(Ve))
    fp <- sum(beta.sim==0 & beta.est !=0)
    fn <- sum(beta.sim!=0 & beta.est ==0)
  }

  fp2 <- sum(C==0 & fit$C !=0)/sum(C==0)
  fn2 <- sum(C!=0 & fit$C !=0)/sum(C!=0)

  # test for rank and output
  if(is.na(rank.est)){
    return(c(Err.C,Err.Y,fp/truezero, fn/truenonzero,0))
  } else{
    if(nrank >= rank.est){
      return(c(Err.C,Err.Y,fp/truezero, fn/truenonzero,0))
    } else{
      return(c(Err.C,Err.Y,fp/truezero, fn/truenonzero,
               100*sum((diag(D1.est)[-(1:nrank)])^2)/sum(D1.est^2))  )
    }
  }
}



dfol <- 'output/simulation/'
fname <- list.files(dfol,full.names = T)
fname %<>% .[grep('simulation3_',.)]

out <- NULL
for (i in 1:length(fname)) {
  load(fname[i])
  U <- input_data$U
  D <- input_data$D
  V <- input_data$V
  C0 <- input_data$C0
  n <- nrow(X)
  # fit_seqT, fit_trac, fit_seqTM, fit_tracM,
  # test_srrr, test_srrrM, fit_seqF, fit_seqFM
  # compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_trac)
  if((sum(fit_seqT$D) >0)){
    compare.fit00 <- rbind(compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_seqT),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_seqF),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_trac),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, test_srrr),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit.seq),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_seqTM),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_seqFM),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit_tracM),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, test_srrrM),
                           compare_taro(U,D,V,C0,n,X, A, Z=NULL, fit.seqM))

    rbind(C0, U%*%D%*%t(V)) - rbind(fit_trac$Z, fit_trac$C)


    colnames(compare.fit00) <- c('Er(C)', 'Er(XC)', 'FPR', 'FNR', 'R%')
    rownames(compare.fit00) <- c('TARO[W+]','TARO','TRAC','CRRR', 'SeCURE',
                                 'TARO*[W+]','TARO*','TRAC*','CRRR*','SeCURE*')

    out <- compare.fit00 %>% data.frame() %>% rownames_to_column('Model') %>%
      mutate(Index = i) %>%
      bind_rows(out,.)
  }
}

pldf <- out %>%
  gather('Performance','value',c(-Model,-Index)) %>%
  mutate(Model = factor(Model, levels = rownames(compare.fit00))) %>%
  mutate(Performance = recode(Performance, 'Er.C.' = 'Er(C)',
                              'Er.XC.' = 'Er(XC)')) %>%
  filter(Performance != 'R.') %>%
  filter(Model != 'TARO[W+]') %>%
  filter(Model != 'TARO*[W+]')

out_table <- pldf %>%
  # filter(Performance != 'R.') %>%
  group_by(Model,Performance) %>%
  summarise(mean = mean(value,na.rm = T)) %>%
  mutate(mean = signif(mean,2)) %>%
  spread(Performance,mean) #%>%
# .[match(rownames(compare.fit00), .$Model),]


## Generate latex table for the manuscript
library(kableExtra)
kbl(out_table,format = 'latex',booktabs = T,  caption = "Setting 3", linesep = "") %>%
  pack_rows("Full data", 1,4) %>%
  pack_rows("Missing data", 5,8) %>%
  kable_styling(latex_options = "striped", stripe_index = c(5:8))


kbl(out_table,booktabs = T,  caption = "Setting 3", linesep = "") %>%
  pack_rows("Full data", 1,4) %>%
  pack_rows("Missing data", 5,8) %>%
  kable_styling(latex_options = "striped", stripe_index = c(5:8)) %>%
  row_spec(c(1,5),bold=T,hline_after = T)


library(ggplot2)
library(ggh4x)
pldf1 <- pldf %>%
  filter(Performance %in% c('Er(C)','Er(XC)'))
pldf1$Missing <- 'No'
pldf1$Missing[grep('[*]',pldf1$Model)] <- 'Yes'

pldf1 %<>% mutate(Method = gsub('[*]','',Model)) %>%
  mutate(Method = factor(Method, levels = c('TARO','TRAC','SeCURE','CRRR')))


plot_view <- ggplot(pldf1 %>% filter(Missing == 'No') , aes(x=Method,y=value, color = Method ))+
  geom_boxplot() + scale_y_log10()  +
  facet_wrap2(vars(Performance), scales = 'free_y', ncol=5) +
  scale_fill_manual(values = c('gray80','gray99')) +
  # facet_grid(key ~ Source, scales = 'free', space = 'fixed') +
  theme_bw() +
  # coord_flip() +
  xlab(NULL) +
  ylab('Error (log-scale)') +
  # ggtitle(sprintf("Comparison of mTraCX score in case of FMT with Donor1+Donor2")) +
  theme( panel.spacing=unit(0.2,"lines"),
         strip.background = element_blank(),
         legend.position = 'top',
         plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
         axis.text.x = element_text(size = 14, color = 'black',angle = 90,vjust = 0.5, hjust = 1.0),
         axis.text.y = element_text(size = 14, color = 'black',angle = 0),
         axis.title = element_text(size = 15, color = 'black'),
         strip.text = element_text(size = 13, color = 'black', angle = 0),
         legend.text = element_text(size = 10, color = 'black', angle = 0),
         legend.background = element_rect(fill="grey95",
                                          size=0.5, linetype="solid"),
         legend.title =  element_text(size = 13, color = 'black', #face = 'bold',
                                      angle = 0, hjust = 0.5) )


setEPS()
postscript('plots/setting3.eps', width = 10, height = 5)
plot_view
dev.off()


