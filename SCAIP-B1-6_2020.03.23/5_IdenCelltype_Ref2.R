###
library(Rcpp)
library("Seurat")
library("SeuratData")
library(harmony)
library(tidyverse)
library(knitr)
library(future)
library(cowplot)
library(grid)
library(gridExtra)
theme_set(theme_grey())
library(data.table)
library(annotables)
library(Matrix)
library(RColorBrewer)
###
outdir <- "./5_IdenCelltype_Ref2_output"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

##################################################
### 2020.5.17, reference panel pbmc3k, pbmcsca ###
##################################################

future::plan(strategy = 'multicore', workers=1)
options(future.globals.maxSize = 10 * 1024 ^ 3)
plan()


######################################################
### 1, Transfer label of reference onto query data ###
######################################################

### 2, pbmc3k
sc <- read_rds("./5_IdenCelltype_output/1_SCAIP.spliced.rds")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method="vst", nfeatures=2000)
sc <- ScaleData(sc, verbose=T)
###
ref <- read_rds("./5_IdenCelltype_output/Reference_pbmc3k.rds")
DefaultAssay(sc) <- "RNA"
DefaultAssay(ref) <- "RNA"

### annotate cell type of query data             
anchors <- FindTransferAnchors(reference=ref, query=sc, dims=1:30)
pred <- TransferData(anchorset=anchors, refdata=ref$celltype2)
###
meta <- sc@meta.data
metaNew2 <- cbind(meta, pred)
###
###output 2
opfn2 <- "./5_IdenCelltype_Ref2_output/Ref2_pbmc3k.cell2.rds" ##
write_rds(metaNew2, opfn2)


### 3, pbmcsca
query <- read_rds("./5_IdenCelltype_output/2_SCAIP.spliced.NormChem.rds")
ref <- read_rds("./5_IdenCelltype_output/Reference_pbmcsca.rds")
#DefaultAssay(query) <- "RNA"
#DefaultAssay(ref) <- "RNA"

### annotate cell type of query data             
anchors <- FindTransferAnchors(reference=ref, query=query, dims=1:30,npcs=30)
pred <- TransferData(anchorset=anchors, refdata=ref$celltype2, dims=1:30)
###
meta <- query@meta.data
metaNew3 <- cbind(meta,pred)

###output 1
opfn3 <- "./5_IdenCelltype_Ref2_output/Meta3.pbmcsca_cell2.rds" ##
write_rds(metaNew3, opfn3)



###########################
### 2, summary results  ###
###########################

### figure 0, score
meta <- read_rds("./5_IdenCelltype_Ref2_output/Meta2.pbmc3k_cell2.rds")
dd <- meta[,c("predicted.id","prediction.score.max")]
fig0 <- ggplot(dd,aes(x=prediction.score.max,col=predicted.id))+
        geom_density()+
        theme_bw()+
        scale_colour_brewer(palette="Set3")
figfn <- "./5_IdenCelltype_Ref2_output/Figure2.0.prediction.Score.pbmc3k_cell2.png"
png(figfn, width=600, height=600, res=105)
fig0
dev.off()


### figure 1, umap
meta <- read_rds("./5_IdenCelltype_Ref2_output/Ref2_pbmc3k.cell2.rds")
dd <- meta[,c("predicted.id","prediction.score.max")]

sc <- read_rds("./4_Harmony_output/2_Norm.Chem.dims50.Cl.rds")
umap <- Embeddings(sc, reduction="umap")
mydf <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                      UMAP_2=as.numeric(umap[,2]),
                      BATCH=sc$BATCH, 
                      chem=sc$chem,
                      treats=sc$treats,
                      cluster=sc$seurat_clusters,
                      MCls=dd$predicted.id)

fig1 <- ggplot(mydf,aes(x=UMAP_1,y=UMAP_2, colour=MCls))+
        geom_point(size=0.1)+
        guides(col=guide_legend(override.aes=list(size=3)))+
        scale_colour_brewer(palette="Set3")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
figfn <- "./5_IdenCelltype_Ref2_output/Figure2.1_umap.pbmc3k.cell2.png"
png(figfn, width=600, height=600, res=105)
fig1
dev.off()

## figure 2, heatmap,
df0 <- mydf%>%group_by(MCls, cluster)%>%
       summarize(Freq=n())%>%
       group_by(cluster)%>%
       mutate(Perc=Freq/sum(Freq)*100)
       
fig2 <- ggplot(df0, aes(x=cluster, y=MCls,fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("Cluster")+ylab("")+
            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
figfn <- "./5_IdenCelltype_Ref2_output/Figure2.2.heatmap.pbmc3k_cell2.png"
png(figfn, width=600, height=600, res=105)
fig2
dev.off()

### heatmap2
## figure 2, heatmap,
df1 <- mydf%>%group_by(MCls, cluster2)%>%
       summarize(Freq=n())%>%
       group_by(cluster2)%>%
       mutate(Perc=Freq/sum(Freq)*100)
       
fig2 <- ggplot(df1, aes(x=cluster2, y=MCls, fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("Cluster")+ylab("")+
            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
figfn <- "./5_IdenCelltype_Ref2_output/Figure2.2.heatmap2.pbmc3k_cell2.png"
png(figfn, width=600, height=600, res=105)
fig2
dev.off()


