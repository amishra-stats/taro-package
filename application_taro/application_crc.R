# Data preprocessing for the analysis
rm(list = ls())
# data folder
dfol <- 'application_taro/'
load(file.path(dfol,'processed_data_crc.rds'))
# Load data
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
library(magrittr)
library(Matrix)
library(phyloseq)
library(trac)
library(reticulate)
# use_python("/Users/amishra/opt/anaconda3/bin/python")
# classo = import("classo")
library(microbiome)
Choice <- 'Genus'
temp <- aggregate_taxa(microbiome_metag,Choice)
microbiome_filt <- subset_taxa(temp,rowSums(otu_table(temp)!=0) > 70)
physeq <- microbiome_filt

## Create A matrix for the analysis
df <- tax_table(physeq) %>%
  data.frame()
col_names <- names(df)
df[col_names] <- lapply(df[col_names] , factor)
collapse= F
formula_str <- paste(setdiff(col_names,'unique'),collapse = '/') %>% paste('~',.,sep = '')
phy <- tax_table_to_phylo(eval(parse(text = formula_str)),
                          data = df, collapse = collapse)
A <- trac::phylo_to_A(phy) %>%
  as.matrix()
all(phy$tip.label == rownames(A))
colnames(A) <- gsub("'",'',colnames(A))
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
X <- otu_table(physeq) %>% t()
X <- log(X + (X==0))
tind <- colSums(X) !=0
X <- X[,tind]; A <- A[tind,]
A <- A[,colSums(A)>0]
A <- apply(A, 1, function(x) x/colSums(A)) %>% t()
Xt <- X %*% A
sum(is.na(cor(Xt))) # sanity check
# Extract metabolomics data as the outcome of the model
mtb_filt <- subset_taxa(metabolom_metag,
                        rowSums(otu_table(metabolom_metag)!=0) > 70)
Y <- otu_table(mtb_filt) %>%
  data.frame() %>% as.matrix() %>% t()
Y <- log(Y + (Y==0))
## Save the data for the analysis
input_data <- list()
input_data$X <- X
input_data$Y <- Y
input_data$A <- A
meta_data <- sample_data(physeq) %>% data.frame()
# Save data for the full data analysis
save(X,Y,A,meta_data,phy,file= file.path(dfol,'processed_data_crc_full_data.rds') )




## --------------------------------------------- Server run full model --------
# We set orthogonality always true
# bsub -W 24:00 -q medium -n 1 -M 16 -R 'rusage[mem=16]' -o out/out_full -e err/err_full "Rscript application-crc.R 50 1>Rout/out_full 2>Rerr/err_full"
# bsub -W 240:00 -q long -n 1 -M 32 -R 'rusage[mem=32]' -o out/out_full -e err/err_full "Rscript application-crc.R 50 1>Rout/out_full 2>Rerr/err_full"
rm(list = ls())
dfol <- 'application_taro'
load(file.path(dfol,'processed_data_crc_full_data.rds'))
# Load data
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
library(magrittr)
library(Matrix)
library(phyloseq)
library(trac)
library(reticulate)
# Setting for calling TARO for model estimation
maxrank = 5;
Z = NULL;
A = A;
Ac <- matrix(1,1,ncol = ncol(X)) %*% A
Bc <- matrix(0,1,1)
nfold = 5; trace = T; verbose = TRUE
nlambda = 30;PATH = T
control <- taro_control(alpha0 = 0, gamma0 = 1, spU = 0.8,
                        outTol=1e-4, outMaxIter=200,
                        inMaxIter=10, inTol=1e-4,
                        spV=0.8, lamMaxFac = 10, se1 = 1)

# Weight Yes
fit_seqT <- taro_path(Y, X, A, Ac, Bc, Z = NULL,
                      maxrank = 10, nlambda = nlambda,
                      control = control,
                      nfold = 5, orthV = TRUE,
                      verbose = TRUE)
save(list = ls(), file = 'output/crc_full_data.RData')









## --------------------------------------------- Output analysis --------
rm(list = ls())
# data folder
dfol <- 'application_taro/'
load(file.path(dfol,'processed_data_crc_full_data.rds'))
# Load data
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
library(magrittr)
library(Matrix)
library(phyloseq)
library(trac)
library(reticulate)
library(tidyverse)
# Implementation setting
maxrank = 5;
Z = NULL;
A = A;
Ac <- matrix(1,1,ncol = ncol(X)) %*% A
Bc <- matrix(0,1,1)
nfold = 5; trace = T; verbose = TRUE
nlambda = 30;PATH = T
control <- taro_control(alpha0 = 0, gamma0 = 1, spU = 0.8,
                        outTol=1e-4, outMaxIter=200,
                        inMaxIter=10, inTol=1e-4,
                        spV=0.8, lamMaxFac = 10, se1 = 1)
# # Weight Yes
# fit_seqT <- taro_path(Y, X, A, Ac, Bc, Z = NULL,
#                       maxrank = 10, nlambda = nlambda,
#                       control = control,
#                       nfold = 5, orthV = TRUE,
#                       verbose = TRUE)
#
# save(list = ls(), file = 'output/crc_full_data.RData')
load('output/crc_full_data.RData')
## Create these latent factors:
LF_MB0 <- X %*% A %*% fit_seqT$U %*% diag(fit_seqT$D)
LF_MTB <- Y %*% fit_seqT$V

pldf <- LF_MB0 %>% data.frame() %>%
  rownames_to_column('ID') %>%
  gather('Components','Microbiome',-ID) %>%
  left_join(.,LF_MTB %>% data.frame() %>%
              rownames_to_column('ID') %>%
              mutate(ID = gsub('X','',ID)) %>%
              gather('Components','Metabolites',-ID))

plot_view <- ggplot(pldf, aes(x=Microbiome, y=Metabolites, color=Components)) +
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed", fill="blue", size = 0.5) +
  theme_bw() +
  guides(color=guide_legend(title="Supervised Latent Factor \n (Microbiome vs Metabolites)")) +
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
         legend.title =  element_text(size = 13, color = 'black', face = 'bold',
                                      angle = 0, hjust = 0.5) )

setEPS()
postscript('application_taro/plots/application/sup_LF.eps', width = 6, height = 6)
plot_view
dev.off()


## Access significance of latent factor constructed from microbiome data
df <- meta_data %>%
  mutate(outcome = ifelse(Study.Group == 'Healthy',1,0)) %>%
  select(outcome,Age,Gender,BMI) %>%
  cbind(LF_MB0) %>% data.frame()

fit.glm <- glm(outcome ~ ., data = df, family = binomial())
summary(fit.glm)

## Access significance of latent factor constructed from metabolites data
df <- meta_data %>%
  mutate(outcome = ifelse(Study.Group == 'Healthy',1,0)) %>%
  select(outcome,Age,Gender,BMI) %>%
  cbind(LF_MTB) %>% data.frame()

fit.glm <- glm(outcome ~ ., data = df, family = binomial())
summary(fit.glm)




## ----------------------------
## Make circular plot for the metabolites contribution
## 2 5 7 components
df_plot <- fit_seqT$V %>%
  data.frame() %>%
  select_at(c(2,5,7)) %>%
  mutate(NodeID = colnames(Y)) %>%
  .[rowSums(.[,1:3])>0,]
df_plot$X5 <- df_plot$X5 * -1
df_plot$NodeID <- sub(".*_","",df_plot$NodeID)
pldf <- df_plot[order(rowSums(abs(df_plot[,1:3])),decreasing = T)[1:60],] %>%
  remove_rownames() %>%
  tibble::column_to_rownames('NodeID')
library(ComplexHeatmap)


## Prepare file for metabolites enrichment analysis using metaboanalyst
fileConn<-file("application_taro/plots/application/output.txt")
temp <- c()
for (i in 1:ncol(pldf) ) {
  temp <- c(temp,colnames(pldf)[i])
  temp <- c(temp,sub(".*_","",rownames(pldf)[pldf[,i] !=0]))
  temp <- c(temp,"\n")
}
writeLines(temp, fileConn)
close(fileConn)


set.seed(123)
setEPS()
postscript('application_taro/plots/application/sup_comp_mtb.eps', width = 12, height = 5)
ht <- Heatmap(t(as.matrix(pldf)), show_row_dend = F, show_column_dend = F,
              cluster_rows = F,
              heatmap_legend_param = list(title = 'Latent factor\nComposition'))
draw(ht, heatmap_legend_side = "left",annotation_legend_side="right")
dev.off()

set.seed(123)
# add circular plot of the analysis
library(circlize)
col_fun1 = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))
setEPS()
postscript('plots/application/sup_comp_mtb_circ.eps', width = 8, height = 8)
circos.par(gap.after = c(12))
circos.heatmap(pldf,  col = col_fun1,track.height = 0.4,
               bg.border = "green", bg.lwd = 2, bg.lty = 2,
               show.sector.labels = F,rownames.cex = 0.7,
               # dend.side = "inside", #rownames.cex = 1,
               rownames.side = "outside", cluster = TRUE)

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = rev(colnames(pldf))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)
circos.clear()
dev.off()







## ----------------------------
## Make circular plot for the microbiome data selected
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(phytools)

phy$node.label <- gsub("'",'',phy$node.label)
sel_comp <- c(2,5,7)
U_est <- fit_seqT$U[,sel_comp, drop = F]
rownames(U_est) <- colnames(A)
# Find the selected nodes and tips
lef <- rownames(A)
nds <- setdiff(colnames(A),lef)
sel_id <- rownames(U_est)[rowSums(abs(U_est)) > 0.1]
node_sel <- intersect(nds,sel_id)
tip_sel <- intersect(lef,sel_id)
nodeids <- nodeid(phy, node_sel)
node_leaf <- NULL
for (i in nodeid(phy, node_sel)) {
  node_leaf %<>%
    c(.,c(phy$tip.label,phy$node.label)[getDescendants(phy, i)])
}
tip_sel <- c(tip_sel, intersect(lef,node_leaf)) %>% unique()


pl_mb <- U_est[sel_id,]
colnames(pl_mb) <- paste('X',sel_comp,sep='')
pl_mb[,'X5'] <- pl_mb[,'X5'] * -1


setEPS()
postscript('plots/application/sup_comp_mbo.eps', width = 16, height = 5)
ht <- Heatmap(t(pl_mb), show_row_dend = F, show_column_dend = F,
              cluster_rows = F,
              heatmap_legend_param = list(title = 'Latent factor\nComposition'))
draw(ht, heatmap_legend_side = "left",annotation_legend_side="right")
dev.off()


## find tips of all the selected node
phy_ref <- ape::drop.tip(phy,setdiff(lef,tip_sel))
nodeids <- nodeid(phy_ref, node_sel)
nodedf <- data.frame(node=nodeids)
labdf <- data.frame(node=nodeids, label=node_sel, pos=1)
tipdf <- data.frame(node=nodeid(phy_ref, tip_sel), label=tip_sel, pos=1)



## Prepare data for heatmap
U_hat <- A %*% U_est
colnames(U_hat) <- paste('X',sel_comp,sep = '')
U_hat[,2] <- -1 * U_hat[,2]
U_plot <- matrix(0,length(phy_ref$tip.label),ncol(U_hat))
rownames(U_plot) <- phy_ref$tip.label
U_plot <- U_hat[rownames(U_plot),,drop = F]
colnames(U_plot) <- colnames(U_hat)
node_dat <- data.frame(ID = rownames(U_plot),
                       Size = rowSums(abs(U_plot))) %>%
  mutate(Selected = ifelse(Size == 0, 'No','Yes')) %>%
  mutate(Size = exp(Size))

U_plot %<>% data.frame() %>%
  rownames_to_column('ID') %>%
  gather('Components','Value',-ID)
labdf$pos = -1# seq(3,4,length.out = nrow(labdf))
tipdf$pos = 2

# The circular layout tree.
p <- ggtree(phy_ref, layout="radial", size=0.15, open.angle=5) +
  # geom_tiplab2() +
  geom_tiplab(size = 2.4,offset = 2) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6, alpha=0.5, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf,
                mapping=aes(node=node,
                            label=label,
                            offset.text=pos),
                hjust= 0.5,
                angle=10,
                barsize=NA,
                horizontal=FALSE,
                fontsize=2.4,
                fontface="italic",
                extend=-0.2,
                offset=0
  )
p <- p + new_scale_fill() +
  geom_fruit(data=U_plot, geom=geom_bar,
             mapping=aes(y=ID, x=Value, fill=Components),
             alpha = 0.5,
             pwidth=0.38,
             orientation="y",
             stat="identity",
  ) +
  guides(fill=guide_legend(title="Selected\nFactor")) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3,
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(.02, 0.1),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=10.5),
        legend.text=element_text(size=8),
        legend.spacing.y = unit(0.02, "cm"),
  )
pdf('plots/application/sup_comp_mb_circ.pdf',
    width = 10, height = 10)
p
dev.off()



## -----------------------------------------------------------------
## Indipendent PCA component analysis of the microbiome and metabolites data
rm(list = ls())
# data folder
dfol <- 'application_taro/'
load(file.path(dfol,'processed_data_crc_full_data.rds'))
# Load data
# devtools::install_github('amishra-stats/taro-package/taro', force = TRUE)
library(taro)
library(magrittr)
library(Matrix)
library(phyloseq)
library(trac)
library(reticulate)
library(tidyverse)
# Implementation setting
maxrank = 5;
Z = NULL;
A = A;
Ac <- matrix(1,1,ncol = ncol(X)) %*% A
Bc <- matrix(0,1,1)
nfold = 5; trace = T; verbose = TRUE
nlambda = 30;PATH = T

## Create these latent factors:
LF_MB0 <- X %*% A
LF_MTB <- Y
library(sparsepca)
# Compute SPCA of microbiome
out_mb <- spca(LF_MB0, k=8, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
print(out_mb)
summary(out_mb)

df_mb <- LF_MB0 %*% out_mb$loadings %>%
  data.frame()


## Principal components of metabolites
out_mt <- spca(LF_MTB, k=8, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
print(out_mt)
summary(out_mt)
df_mt <- LF_MTB %*% out_mt$loadings %>%
  data.frame()
temp <- df_mt

## Align principal components using colinear way
for (i in 1:ncol(df_mt)) {
  ind <- which.max(abs(cor(df_mb[,i,drop=F],temp)))
  df_mt[,i] <- temp[,ind]; temp <- temp[,-ind,drop=F]
}

## make a plot of the components
pldf <- df_mb %>% data.frame() %>%
  rownames_to_column('ID') %>%
  gather('Components','Microbiome',-ID) %>%
  left_join(.,df_mt %>% data.frame() %>%
              rownames_to_column('ID') %>%
              mutate(ID = gsub('X','',ID)) %>%
              gather('Components','Metabolites',-ID))

plot_view <- ggplot(pldf, aes(x=Microbiome, y=Metabolites, color=Components)) +
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed", fill="blue", size = 1) +

  theme_bw() +
  guides(color=guide_legend(title="Independent PCA component \n (Microbiome vs Metabolites)")) +
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
         legend.title =  element_text(size = 13, color = 'black', face = 'bold',
                                      angle = 0, hjust = 0.5) )

setEPS()
postscript('application_taro/plots/application/independent_pca.eps', width = 6.5, height = 6.5)
plot_view
dev.off()


## Sparse PCA plots
df <- meta_data %>%
  mutate(outcome = ifelse(Study.Group == 'Healthy',1,0)) %>%
  select(outcome,Age,Gender,BMI) %>%
  cbind(df_mt) %>% data.frame()
fit.glm <- glm(outcome ~ ., data = df, family = binomial())
summary(fit.glm)


df <- meta_data %>%
  mutate(outcome = ifelse(Study.Group == 'Healthy',1,0)) %>%
  select(outcome,Age,Gender,BMI) %>%
  cbind(df_mb) %>% data.frame()
fit.glm <- glm(outcome ~ ., data = df, family = binomial())
summary(fit.glm)

