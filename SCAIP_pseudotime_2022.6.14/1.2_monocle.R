##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(monocle3)

##
rm(list=ls())

outdir <- "./1.2_monocle.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)



###
### run monocle
sc <- read_rds("./0_data.outs/0_Tcell.PHA.chem3.rds")

count <- sc@assays$RNA@counts
meta <- sc@meta.data
gene_anno <- data.frame(gene_short_name=rownames(count))
rownames(gene_anno) <- rownames(count)

## cell_data_set object
cds <- new_cell_data_set(count, cell_metadata=meta, gene_metadata=gene_anno)

### pre-process the data
cds <- preprocess_cds(cds, num_dim=30)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition=F)

figfn <- "./1.2_monocle.outs/Figure0.1_monocle.png"
png(figfn, width=480, height=480, res=100)
plot_cells(cds, color_cells_by="treat2", label_cell_groups=F,
           label_leaves=T, label_branch_points=T, graph_label_size=1.5)
dev.off()

##
get_earliest_principal_node <- function(cds, treat0="PHA"){
   ### 
   cell_ids <- which(colData(cds)[, "treat2"] == treat0)
    
   closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
   root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
###
   root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes="Y_365")

figfn <- "./1.2_monocle.outs/Figure0.3_pseudotime.png"
png(figfn, width=480, height=480, res=100)
plot_cells(cds, color_cells_by="pseudotime", label_cell_groups=F,
           label_leaves=T, graph_label_size=1.5)
dev.off()


##################
### simulation ###
##################

###
##
get_earliest_principal_node <- function(cds, treat0="PHA"){
   ### 
   cell_ids <- which(colData(cds)[, "treat2"] == treat0)
    
   closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
   root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
###
   root_pr_nodes
}

iid <- read.csv("./0_data.outs/0.1_iid.csv")
sc <- read_rds("./0_data.outs/0_Tcell.PHA.chem3.rds")
 
sc <- NormalizeData(sc)
sc <- ScaleData(sc, features=rownames(sc), verbose=T)

count <- sc@assays$RNA@counts
NormData <- sc@assays$RNA@data    
meta <- sc@meta.data
## gene meta data
gene_anno <- data.frame(gene_short_name=rownames(count))
rownames(gene_anno) <- rownames(count)


### 50 replicates
corr <- mclapply(1:50, function(i){
   ##
   cat("id", i, "\n") 
   id <- iid[,i]
   count0 <- count[,id]
   meta0 <- meta[id,]
   norm <- as.matrix(NormData[,id])
  ## cell_data_set object
  cds <- new_cell_data_set(count0, cell_metadata=meta0, gene_metadata=gene_anno)
  ##
  cds <- cds%>%preprocess_cds(num_dim=30)%>%
     reduce_dimension()%>%
     cluster_cells()%>%
     learn_graph(use_partition=F)
  ###
  root <- get_earliest_principal_node(cds)
  cds <- order_cells(cds, root_pr_nodes=root)
  pseudo <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  ##
  rr <- cor(t(norm), pseudo, method="pearson")
  rr
},mc.cores=20)

corr <- do.call(cbind, corr)
corr <- as.data.frame(corr)
corr <- cbind(rownames(corr), corr)
names(corr) <- c("gene", paste("id", 1:50, sep=""))

opfn <- "./1.2_monocle.outs/1_corr_pseudotime.csv"
write.csv(corr, opfn, row.names=F)
    
  
