##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
## library(SeuratDisk)
## library(SeuratWrappers)
## library(monocle3)


sc <- read_rds("./0_data.outs/0_Tcell.PHA.chem3.rds") 
rn <- rownames(sc)
sc <- sc%>%NormalizeData()%>%
    ScaleData(features=rn, verbose=T)

Norm <- sc@assays$RNA@data

iid <- read.csv("./0_data.outs/0.1_iid.csv")

pseu <- read.csv("./2_scanpy.outs/1_pseudotime.csv")

corr <- lapply(1:50, function(i){
    ##
    cat("id", i, "\n")
    id <- iid[,i]
    norm <- as.matrix(Norm[,id])
    t0 <- pseu[,i]
    norm2 <- norm[,is.finite(t0)] 
    t1 <- t0[is.finite(t0)]
    ##
    rr <- cor(t(norm2), t1, method="pearson")
    rr
})

###
corr <- do.call(cbind, corr)
corr <- as.data.frame(corr)
corr <- cbind(rownames(corr), corr)
names(corr) <- c("gene", paste("id", 1:50, sep=""))

write.csv(corr, "./2_scanpy.outs/2_corr_pseudotime.csv", row.names=F)
                      
