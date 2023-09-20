

# TARO: Tree aggregated factor regression model for microbiome data analysis


We propose a Tree-Aggregated factor Regression mOdel (TARO) for integrating microbiome data with other high-dimensional data types, such as metabolomics. Technical limitations prevent us from obtaining the absolute count of the microbial species; hence the microbial-abundance profile of a sample is inherently compositional. In addition, microbial species are related by phylogeny. TARO treats the microbial abundance data as compositional data and suitably encodes the dependency among the ASV/OTUs through a phylogeny-inspired adjacency matrix.

The adjacency matrix allows us to consider leaf and node in the phylogenetic tree as predictors in the model and learn their association with the multivariate response in terms of a low-rank and sparse coefficient matrix. The required regularized structure of the coefficient matrix allows us to identify multiple latent factors (each represented as a subset of predictors) associated with only a subset of responses. 

We demonstrate through simulation studies that TARO can accurately recover the low-rank coefficient matrix and identify relevant features. 

![alt text](https://github.com/amishra-stats/taro-package/blob/main/misc/schema_taro.jpg)



Getting started  
--------------
The `taro` package is currently available on GitHub and can be installed as follows.

```
# Install packages
devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
# load library
library(taro)
```

We have implemented the procedure in the function *taro_path* that primarily requires a multivariate response matrix, log-transformed microbial abundance matrix as predictors and a phylogeny inspired adjacency matrix.

## Working examples

This is a basic example which shows you how to solve a common problem. The code below relies on SparseDOSSA2 to create the simulated data. SparseDOSSA2 can be installed as follows:

``` r
devtools::install_github("biobakery/SparseDOSSA2")
```

Once SparseDOSSA2 is installed, the following example code can be run to simulate the data and fit the TARO model:

``` r
library(taro)
library(SparseDOSSA2)
# Setting parameters for simulating the data
snr <- .25; xrho <- 0.5; nrank <- 3; 
q <- 50; n = 300; intercept = 0.5
rho = 0.5 # error correlation

# Simulate data 
input_data <- taro_sim(n, q, nrank, rho, snr, intercept,
                        taxa_table = NULL)
Y <- input_data$Y
X <- input_data$X
A <- input_data$A ## Phylogeny inspired adjacency matrix 

## Model fitting 
set.seed(123)
n <- nrow(Y); q <- ncol(Y)
maxrank = 5;
Z = NULL;
A = A;
Ac <- matrix(1,1,ncol = ncol(X)) %*% A
Bc <- matrix(0,1,1)
nfold = 5; trace = TRUE; verbose = TRUE
nlambda = 100;PATH = TRUE
# Control parameters for obtaining the model fit 
control <- taro_control(alpha0 = 0.5, gamma0 = 1, spU = 0.5,
                       inTol = 1e-5, inMaxIter = 300,
                       outMaxIter = 1000,outTol = 1e-8,
                       spV=0.5, lamMaxFac = 1e2, se1 = 1)
# Weight Yes seeting 
fit_seq <- taro_path(Y, X, A, Ac, Bc, Z = Z,
                      maxrank = maxrank, nlambda = nlambda,
                      control = control,
                      nfold = nfold, orthV = TRUE,
                      verbose = TRUE)
                    

```

We plot the model output using the parameters estimate. The loading matrix estimate of ***V*** and the left-singular vector components for constructing the latent factors  ***AU*** are presented below, respectively. 

``` r
library(magrittr)
library(tibble)
library(ComplexHeatmap)
df_plot <- fit_seq$V %>% data.frame() %>%
  dplyr::mutate(NodeID = 1:ncol(Y)) %>% .[rowSums(.[,1:3])!=0,] %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames('NodeID')
jpeg('../misc/sup_comp_loadings.jpg', width = 800, height = 300)
ht <- Heatmap(t(as.matrix(df_plot)), show_row_dend = F, show_column_dend = F,
              cluster_rows = F, cluster_columns = F,
              heatmap_legend_param = list(title = 'Loading Matrix\nComposition'))
draw(ht, heatmap_legend_side = "left",annotation_legend_side="right")
dev.off()


df_plot <- A %*% fit_seq$U %>% data.frame() %>%
  .[rowSums(.[,1:3])!=0,] 
jpeg('../misc/sup_comp_latent_fac.jpg', width = 1000, height = 300)
ht <- Heatmap(t(as.matrix(df_plot)), show_row_dend = F, show_column_dend = F,
              cluster_rows = F, cluster_columns = T,
              heatmap_legend_param = list(title = 'Latent Factor\nMatrix'))
draw(ht, heatmap_legend_side = "left",annotation_legend_side="right")
dev.off()
```

![alt text](https://github.com/amishra-stats/taro-package/blob/main/misc/sup_comp_latent_fac.jpg)

![alt text](https://github.com/amishra-stats/taro-package/blob/main/misc/sup_comp_loadings.jpg)





Community Guidelines
--------------------

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution guidelines](https://github.com/mingzehuang/latentcor/blob/master/CONTRIBUTING.md) prior
    to submitting a pull request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/mda-primetr/mtracx/issues).
3.  Contributors must adhere to the [Code of Conduct](https://github.com/amishra-stats/latentcor/blob/master/CODE_OF_CONDUCT.md).

Acknowledgments
--------------

We thank Wargo Lab members for useful comments on the project.

## Inquiries

You can also contact us via email

- [Aditya Mishra](mailto:akmishra@mdanderson.org)
- [Christine Peterson](mailto:CBPeterson@mdanderson.org)
