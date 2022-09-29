##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)

##
outdir <- "./0_data.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### generate data
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2_output/Filter2_DEG6571/2_Batch1456.DEG.rds"
sc <- read_rds(fn)


### extract data, Tcells chem3, PHA or PHA-DEX 
sc2 <- subset(sc, subset=treats%in%c("PHA-EtOH", "PHA-DEX")&MCls=="Tcell"&chem=="V3")
## counts
count <- sc2@assays$RNA@counts
rownames(count) <- gsub("S-|\\..*", "", rownames(count))

## cell meta data
meta <- sc2@meta.data%>%mutate(treat2=gsub("-EtOH", "", treats))%>%
   dplyr::select(nCount_RNA, nFeature_RNA, total_counts, NEW_BARCODE, SNG.BEST.GUESS, treat2, MCls)
 
## gene meta data
gene_anno <- data.frame(gene_short_name=rownames(count))
rownames(gene_anno) <- rownames(count)

## seurat object
scNew <- CreateSeuratObject(counts=count, project="SCAIP", meta.data=meta)
## scNew <- ScaleData(scNew, features=rownames(count), verbose=T) 
opfn <- "./0_data.outs/0_Tcell.PHA.chem3.rds"
write_rds(scNew, opfn)


### index of subsets
iid <- map_dfc(1:100, function(i){
   ##
   ncell <- ncol(scNew) 
   subi <- sample(1:ncell, size=round(0.5*ncell))
   subi
})
##
iid <- as.matrix(iid)
colnames(iid) <- paste("id", 1:100, sep="")
write.csv(iid, file="./0_data.outs/0.1_iid.csv", row.names=F)



### convert data into AnnData
sc <- read_rds("0_data.outs/0_Tcell.PHA.chem3.rds") 
SaveH5Seurat(sc, filename="./0_data.outs/1_Tcell.PHA.h5Seurat")
Convert("./0_data.outs/1_Tcell.PHA.h5Seurat", dest="h5ad")




