###
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)
##
library(Seurat)
## library(ggplot2)
## library(cowplot)
## library(grid)
## library(gridExtra)
## library(ggExtra)
## library(RColorBrewer)
## library(ComplexHeatmap)
## library(circlize)
## library(viridis)
## library(ggrastr)
## theme_set(theme_grey())

rm(list=ls())

outdir <- "./3_LDA.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


###
### seurat object
sc <- read_rds("./0_data.outs/0_Tcell.PHA.chem3.rds") 
rn <- rownames(sc)
sc <- sc%>%NormalizeData()%>%
    ScaleData(features=rn, verbose=T)

X <- sc@assays$RNA@data
mu <- apply(X, 1, mean)
Xc <- sweep(X, 1, mu, "-") 
ncell <- ncol(X)
var <- rowSums(Xc*Xc)/(ncell-1)
names(var) <- rn

meta <- sc@meta.data


###
### differential analysis
res <- read_rds("../SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
res <- res%>%mutate(gamma=ifelse( (abs(beta)>0.5)&(qval<0.1), 1, 0))
###
res0 <- res%>%filter(MCls=="Tcell", contrast=="PHA-DEX")
gamma <- res0$gamma
names(gamma) <- res0$gene
gamma[is.na(gamma)] <- 0
s0 <- gamma[rn]
##


###
###
iid <- read.csv("./0_data.outs/0.1_iid.csv")


###
### Simulation 
corr <- mclapply(1:50, function(i){
###
   cat("id", i, "\n")
   id <- iid[,i]
   X0 <- X[,id]
   meta0 <- meta[id,]
   ncell <- ncol(X0)
    
   mu <- rowMeans(X0)
   xc <- sweep(X0, 1, mu, "-") 
   var <- rowSums(xc*xc)/(ncell-1)
    
   ### matrix 1 
   ii <- meta0%>%filter(treat2=="PHA")%>%dplyr::pull(NEW_BARCODE)
   X1 <-X0[,ii]
   mu1 <- apply(X1, 1, mean)

   ### matrix 2
   ii <- meta0%>%filter(treat2=="PHA-DEX")%>%dplyr::pull(NEW_BARCODE)
   X2 <-X0[,ii]
   mu2 <- apply(X2, 1, mean)


   Diff <- as.matrix((mu2-mu1)*(1/var)*s0)
   Diff[is.na(Diff)] <- 0

   z1 <- as.vector(crossprod(xc, Diff))
   ##
   xt <- t(as.matrix(X0)) 
   rr <- cor(xt, z1, method="pearson")
   rr
},mc.cores=1)
 
corr <- do.call(cbind, corr)
corr <- as.data.frame(corr)
corr <- cbind(rownames(corr), corr)
names(corr) <- c("gene", paste("id", 1:50, sep=""))

opfn <- "./3_LDA.outs/1_corr_pseudotime.csv"
write.csv(corr, opfn, row.names=F)
