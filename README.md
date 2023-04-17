

# TARO: Tree aggregated factor regression model for microbiome data analysis

The abundance profiles observed in microbiome samples are compositional because of several technical limitations in 16S rRNA/metagenomics sequencing. The R package `TARO` uses Earth Movers Distance (EMD) to compare the microbiome compositions of any two states. 

Using this package, we have assessed the extent of microbiome colonization in a fecal microbiome transplant experiment. The underlying structure of microbiome transfer in an FMT experiment are illustrated in the Figure below:

![alt text](https://github.com/mda-primetr/mtracx/raw/main/misc/figure/schema.png)



Getting started  
--------------
The `taro` package is currently available on GitHub and can be installed as follows.

```
# Install packages
devtools::install_github('mda-primetr/mtracx/mtracx', force = TRUE)
# load library
library(mtracx)
```

We have inplmented the procedure in the function *mtracx_phyloseq* that operates on *phyloseq object* (storing microbial abundandance profiles of all relevant samples) and  *experiment_file* (a data frame providing details of the FMT experiment). We have provided an example of the experimental file in the figure below: 

![alt text](https://github.com/mda-primetr/mtracx/raw/main/misc/figure/experiment_file.png)


## Working examples

This is a basic example which shows you how to solve a common problem:

``` r
library(taro)
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
