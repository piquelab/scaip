##
rm(list=ls())
source("./Bin/LibraryPackage.R")

###  
outdir <- "./5_IdenCelltype_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F) 


################################################################
### 2020.11.18, Using seurat to deconvolution, by Julong wei ### 
###                  used in downstream analysis             ###
################################################################


#future::plan(strategy = 'multicore', workers = 10)
#options(future.globals.maxSize = 40 * 1024 ^ 3)
#plan()




##################### 
### 1 query data  ###
##################### 
          
if(FALSE){
cat("1.", "rebuild query data", "\n")
grchUnq <- grch38%>%
           distinct(ensgene,.keep_all=T)%>%
           dplyr::select(ensgene, symbol, chr, start, end) ##63697
### ensembel id have same symbol name
###spliced genes
sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
meta <- sc@meta.data
x <- sc@assays$RNA@counts

anno <- tibble(rn=rownames(x))%>% 
           mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), uns=grepl("S-",rn), rnz=rowSums(x))%>%
           filter(uns, rnz>0)%>%
           left_join(grchUnq, by="ensgene")%>%
           drop_na(symbol)
                      
xs <- x[anno$rn,]
rownames(xs) <- anno$symbol
scNew <- CreateSeuratObject(xs, min.cells=0, min.features=0, project="SCAIP.spliced")
scNew <- AddMetaData(object=scNew, metadata=sc@meta.data)

write_rds(scNew, "./5_IdenCelltype_output/1_SCAIP.spliced.rds")

###
#scNew <- read_rds("./5_IdenCelltype_output/1_SCAIP.spliced.rds")
#sp <- SplitObject(scNew, split.by = "chem")
#sp <- lapply(X = sp, FUN = function(x) {
#    x <- NormalizeData(x)
#    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})

#anchors <- FindIntegrationAnchors(object.list=sp, dims=1:30)
#combined <- IntegrateData(anchorset = anchors, dims=1:30)
#sc2 <- ScaleData(combined,verbose=T)
sc <- read_rds("./5_IdenCelltype_output/1_SCAIP.spliced.rds")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method="vst", nfeatures=2000)
sc <- ScaleData(sc, verbose=T)
sc <- RunPCA(sc, npcs=100, verbose=T)
sc <- RunHarmony(sc, "chem", reduction="pca")
###


###output1
opfn1 <- "./5_IdenCelltype_output/2_SCAIP.spliced.NormChem.rds"
write_rds(sc2, opfn1)

} ###




#########################
### 2, Reference data ###
#########################

if(FALSE){

cat("2.", "construct reference data", "\n")
###reference 1, Zheng68k
ref.data <- Read10X(data.dir="../SCAIP-ALL-2019.10.24/PBMCs/filtered_matrices_mex/hg19/")
ref <- CreateSeuratObject(ref.data, min.cells=0, min.features=0, project="pbmc68k")
###
cell <- read.table(file="../SCAIP-ALL-2019.10.24/PBMCs/68k_pbmc_barcodes_annotation.tsv", sep="\t", header=T)  ###

metaNew <- cbind(ref@meta.data,cell)
x2 <- as.character(metaNew$celltype)
x2 <- gsub("CD4\\+.*", "CD4+", x2)
x2 <- gsub("CD8\\+.*", "CD8+", x2)
metaNew$celltype2 <- x2

ref2 <- AddMetaData(object=ref, metadata=metaNew)
 
### normalized data
ref2 <- NormalizeData(object=ref2)
ref2 <- FindVariableFeatures(object=ref2, selection.method="vst", nfeatures=2000)
ref2 <- ScaleData(ref2)
ref2 <- RunPCA(ref2, verbose=T, npcs=100)
###
opfn <- "./5_IdenCelltype_output/Reference_Zheng68k.rds"
write_rds(ref2, opfn)

#ref <- read_rds("../SCAIP-ALL-2019.10.24/5_Identify_CelltypeNew_output/Reference_Zheng64k.rds")
#write_rds(ref, "./5_IdenCelltype_output/Reference_Zheng68k.rds")
#ref <- read_rds("./5_IdenCelltype_output/Reference_Zheng68k.rds")

###reference 2
#InstallData("pbmc3k")
#data("pbmc3k")
#ref <- NormalizeData(object=pbmc3k)
#ref <- FindVariableFeatures(object=ref, selection.method="vst", nfeatures=2000)
#ref <- ScaleData(ref)  
###    
#opfn <- "./5_IdenCelltype_output/Reference_pbmc3k.rds"
#write_rds(ref,opfn)

#fn <- "./5_IdenCelltype_output/Reference_pbmc3k.rds"
#ref <- read_rds(fn)
#meta <- ref@meta.data
#x2 <- as.character(meta$seurat_annotations)
#x2 <- gsub(".*CD4.*", "CD4+", x2)
#x2 <- gsub(".*CD8.*", "CD8+", x2)
#x2 <- gsub(".*Mono.*", "Mono", x2)

#meta$celltype2 <- x2
#ref2 <- AddMetaData(object=ref, metadata=meta)

#write_rds(ref2, "./5_IdenCelltype_output/Reference_pbmc3k.rds")
#meta <- ref@meta.data
#dd <- meta%>%group_by(seurat_annotations)%>%summarize(ncell=n())
#write.csv(dd, file="./5_IdenCelltype_output/tmp/Ref2.pbmc3k.celltype.csv", row.names=F)


### reference 3
##InstallData("pbmcsca")
#data("pbmcsca")
#sp <- SplitObject(pbmcsca, split.by = "Method")
#sp <- lapply(X = sp, FUN = function(x) {
#    x <- NormalizeData(x)
#    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})
#anchors <- FindIntegrationAnchors(object.list=sp, dims=1:30)
#combined <- IntegrateData(anchorset = anchors, dims=1:30)

#ref <- ScaleData(combined,verbose=T)
#opfn <- "./5_IdenCelltype_output/Reference_pbmcsca.rds"
#write_rds(ref,opfn)
##
#meta <- ref@meta.data
#dd <- meta%>%group_by(CellType)%>%summarize(ncell=n())
#write.csv(dd, file="./5_IdenCelltype_output/tmp/Ref3.pbmcsca.celltype.csv", row.names=F) 

#ref <- read_rds("./5_IdenCelltype_output/Reference_pbmcsca.rds")
#meta <- ref@meta.data 
#x2 <- as.character(meta$CellType)
#x2 <- gsub(".*monocyte", "Mono", x2)
#x2 <- gsub("B.*", "B", x2)
#x2 <- gsub(".*dritic.*","DC", x2)
#x2 <- gsub("Natural.*", "NK", x2)
#x2 <- gsub("CD4\\+.*", "CD4+", x2)
#x2 <- gsub("Cyto.*", "CD8+", x2)
#x2 <- gsub("Mega.*", "Mega", x2)

#meta$celltype2 <- x2
#ref2 <- AddMetaData(object=ref, metadata=meta)

#write_rds(ref2, "./5_IdenCelltype_output/Reference_pbmcsca.rds")

} ###


#########################
### 3, transfer label ###
#########################
if(FALSE){

cat("3.", "transfer label", "\n")

### read query data and reference data
#query <- read_rds("./5_IdenCelltype_output/2_SCAIP.spliced.NormChem.rds")
query <- read_rds("./5_IdenCelltype_output/1_SCAIP.spliced.rds")
ref <- read_rds("./5_IdenCelltype_output/Reference_Zheng68k.rds")
#DefaultAssay(query) <- "RNA"


### annotate cell type of query data             
anchors <- FindTransferAnchors(reference=ref, query=query, dims=1:30)
pred <- TransferData(anchorset=anchors, refdata=ref$celltype2, dims=1:30)
###
meta <- query@meta.data
metaNew <- cbind(meta,pred)

###output 1
opfn2 <- "./5_IdenCelltype_output/3_Meta.Zheng68k.rds" ## default
write_rds(metaNew, opfn2)


###option-2, transfer label
rm(list=ls())
query <- read_rds("./5_IdenCelltype_output/1_SCAIP.spliced.rds")
query.ls <- SplitObject(query, split.by = "chem")
ref <- read_rds("./5_IdenCelltype_output/Reference_Zheng68k.rds")

meta <- lapply(query.ls,function(x){
   anchors <- FindTransferAnchors(reference=ref, query=x, dims=1:30)
   pred <- TransferData(anchorset=anchors, refdata=ref$celltype2, dims=1:30)
   tmp <- cbind(x@meta.data,pred)
})
metaNew <- do.call(rbind, meta)

###
opfn2 <- "./5_IdenCelltype_output/3_Meta2.Zheng68k.rds" ##
write_rds(metaNew, opfn2) 


#query <- AddMetaData(query, metadata=metaNew)
#opfn2 <- "./5_IdenCelltype_output/3_SCAIP.anno.Zheng68k.rds" ##
#write_rds(query, opfn2)
### integrate new annotation into previous data  
#query <- read_rds("./5_IdenCelltype_output/3_SCAIP.anno.Zheng68k.rds")
#sc <- read_rds("./4_OldNormAndChem_output/4_OldNorm.Chem.pc50.Cl.rds")
#meta$MCls <- query$predicted.id
#sc2 <- AddMetaData(sc,metadata=meta)
#opfn3 <- "./5_IdenCelltype_output/4_SCAIP.MCls.Zheng68k.rds"
#write_rds(sc2, opfn3)

###
} ####



#######################
### 4, show figures ###
#######################
if(FALSE){

cat("4.", "Summary annotaion results", "\n")

### figure 1, distribution score
meta <- read_rds("./5_IdenCelltype_output/3_Meta.Zheng68k.rds")
dd <- meta[,c("predicted.id","prediction.score.max","NEW_BARCODE")]
fig1 <- ggplot(dd,aes(x=prediction.score.max,col=predicted.id))+
        geom_density()+
        theme_bw()+
        ggtitle("Prediction score across cell type")+
        scale_colour_brewer(palette="Set3")+
        theme(plot.title=element_text(hjust=0.5))
figfn <- "./5_IdenCelltype_output/Figure1.1.prediction.Score.png"
png(figfn, width=700, height=600, res=110)
fig1
dev.off()

### figure 2, umap by identified cell type
sc <- read_rds("./4_Harmony_output/2_Norm.Chem.dims50.Cl.rds")
umap <- Embeddings(sc, reduction="umap")

mydf <- data.frame(NEW_BARCODE=as.character(rownames(meta)),
                   UMAP_1=as.numeric(umap[,1]), 
                      UMAP_2=as.numeric(umap[,2]),
                      BATCH=sc$BATCH, 
                      chem=sc$chem,
                      treats=sc$treats,
                      cluster=sc$seurat_clusters,
                      MCls=dd$predicted.id)

fig2 <- ggplot(mydf,aes(x=UMAP_1,y=UMAP_2, colour=MCls))+
        geom_point(size=0.1)+
        guides(col=guide_legend(override.aes=list(size=3)))+
        scale_colour_brewer(palette="Set3")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
figfn <- "./5_IdenCelltype_output/Figure1.2.umap.Zheng68k.png"
png(figfn, width=700, height=600, res=110)
fig2
dev.off()

### figure for split
#mydf <- data.frame(NEW_BARCODE=as.character(rownames(umap)),
#                   UMAP_1=as.numeric(umap[,1]), 
#                   UMAP_2=as.numeric(umap[,2]),
#                   BATCH=sc$BATCH, 
#                   chem=sc$chem,
#                   treats=sc$treats,
#                   cluster=sc$seurat_clusters)%>%
#                   inner_join(dd,by="NEW_BARCODE")
#mydf <- mydf%>%rename(MCls=predicted.id)

#fig2 <- ggplot(mydf,aes(x=UMAP_1,y=UMAP_2, colour=MCls))+
#        geom_point(size=0.1)+
#        guides(col=guide_legend(override.aes=list(size=3)))+
#        scale_colour_brewer(palette="Set3")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA))
              
#figfn <- "./5_IdenCelltype_output/FigureSplit1.2.umap.Zheng68k.png"
#png(figfn, width=700, height=600, res=110)
#fig2
#dev.off()
#figfn <- "./5_IdenCelltype_output/Figure1.1.Zheng68k.png"
#png(figfn, width=1000,height=600, res=150)
#plot_grid(fig1, fig2, ncol=2)
#dev.off()


## figure 2, heatmap,
df0 <- mydf%>%group_by(MCls, cluster)%>%
       summarize(Freq=n())%>%
       group_by(cluster)%>%
       mutate(Perc=Freq/sum(Freq)*100)
       
fig3 <- ggplot(df0, aes(x=cluster, y=MCls,fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("")+ylab("")+
            theme_bw()+theme(axis.text.x=element_text(hjust=1))
figfn <- "./5_IdenCelltype_output/Figure1.3.heatmap.Zheng68k.png"
png(figfn, width=800, height=600, res=120)
fig3
dev.off()


###figure 4 
#library(ggalluvial)
#library(ggrepel)

#fig4 <- ggplot(mydf,aes(y=Freq, axis1=seurat, axis2=singleR))+
#            geom_alluvium(aes(fill=seurat),width=0.04,knot.pos=0, reverse=FALSE)+
#            guides(fill=FALSE)+
#            geom_stratum(width = 0.04, reverse = FALSE) +
#            scale_x_continuous(breaks=0:3, labels=c("","seurat","singleR",""))+
#            xlim(-0.3,3.5)+
#            geom_label_repel(stat = "stratum", reverse = FALSE,size=3,label.strata=T,
#                                                      nudge_x=c(rep(-0.3,6),rep(0.6,6)),direction="y",
#                                                      hjust=c(rep(1,6),rep(0,6)),vjust=0.5,segment.size=0.3)+
#            theme(legend.position="none",
#                                panel.grid.major=element_blank(),
#                                panel.grid.minor=element_blank(),
#                                panel.background=element_blank(),
#                                axis.title=element_blank(),
#                                axis.ticks=element_blank(),
#                                axis.text=element_blank())
              #axis.text.y=element_blank(),
              #axis.text.x=element_text(size=12, face = "bold"))
###
#figfn <- "./5_Identify_CelltypeNew_output/SingleR.2.alluvial.pdf"
#pdf(figfn,height=7,width=7)
#fig4
#dev.off()
} ###


###################################
### 5, new cell type annotation ###
###################################

if(FALSE){ 

cat("5.1.", "re-annotate cell type", "\n")

sc <- read_rds("./4_Harmony_output/2_Norm.Chem.dims50.Cl.rds")
#Idents(sc) <- sc$integrated_snn_res.0.15
meta <- sc@meta.data

ii <- Idents(sc)
MCls <- rep("Tcell", nrow(meta))
MCls[ii==1] <- "NKcell"
MCls[ii==2] <- "Bcell"
MCls[ii%in%c(3,6)] <- "Monocyte"

meta$MCls <- MCls
sc <- AddMetaData(sc,metadata=meta)

opfn <- "./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds"
write_rds(sc,opfn)

sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
meta <- sc@meta.data
write.csv(meta, file="./5_IdenCelltype_output/5_Harmony.meta")

} ###

####################
### show figures ###
#################### *** used for paper
if(FALSE){
rm(list=ls())

cat("5.2.", "UMAP", "\n")

### (1), UMAP   

cat("(1).", "UMAP colored by cell type", "\n")
sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")

#col0 <- c("Tcell"="#ff7f00", "NKcell"="#a65628", 
#          "Bcell"="#4daf4a", "Monocyte"="#984ea3")

col0 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #NK for"#a63728" alternative
fig1 <- DimPlot(sc, label=T, group.by="MCls", cols=col0)+
        theme_bw()+
        theme(legend.position="none",
              axis.title=element_text(size=12),
              axis.text=element_text(size=12),
              #legend.text=element_text(size=10), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(),            
              panel.border=element_rect(colour="black", fill=NA))
png("./5_IdenCelltype_output/Figure2.1.umap.png", width=500, height=600, res=100)
print(fig1)
dev.off()

###
### (2) combine UMAP
#sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")

cat("(2).", "UMAP colored by Cell type and treatments", "\n")

umap <- sc@reductions$umap@cell.embeddings
df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                  UMAP_2=as.numeric(umap[,2]),
                  MCls=sc@meta.data$MCls, treats=sc@meta.data$treats)
                  
df2$treat2 <- gsub("-EtOH", "", df2$treats)

### umap.1
col0 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
p0 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2))+
        geom_point(aes(colour=MCls), size=0.1)+
        scale_colour_manual(values=col0,guide=guide_legend(override.aes=list(size=1)))+
        #guides(colour=guide_legend(override.aes=list(size=1)))+
        theme_bw()+
        theme(axis.title=element_text(size=8),
              axis.text=element_text(size=8),
              legend.position=c(0.2, 0.85),
              legend.title=element_blank(),
              legend.background=element_blank(),
              legend.key=element_blank(),  #unit(3,"cm") colour="transparent",
              legend.key.size=grid::unit(1,"lines"),
              legend.text=element_text(size=8),
              legend.box.background=element_blank())        


### umap.2
col1 <- c("CTRL"="#828282", 
           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")


Labalpha <- c("LPS"="a","LPS-DEX"="b", "PHA"=c, "PHA-DEX"="d", "CTRL"="e") 
df2$treat1 <- Labalpha[df2$treat2]
Labtreat <- c("a"="LPS", "b"="LPS+DEX", "c"="PHA", "d"="PHA+DEX", "e"="CTRL")

p1 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2))+      
           stat_density_2d(aes(fill=stat(level)), geom="polygon", contour=T)+
           facet_wrap(~factor(treat1), nrow=3, labeller=as_labeller(Labtreat))+
           scale_fill_viridis_c(direction=-1)+theme_bw()+
           theme(axis.text=element_text(size=8),
                 axis.title=element_text(size=8),
                 legend.position=c(0.8,0.15),
                 legend.title=element_blank(),
                 legend.background=element_rect(colour=NA, fill=NA),
                 legend.key=element_rect(fill=NA),
                 legend.key.size=grid::unit(0.25,"cm"),
                 legend.text=element_text(size=6),
                 strip.text=element_text(size=8))

              
png("./5_IdenCelltype_output/Figure2.2.umap.png", width=1000, height=600,res=130)
print(plot_grid(p0, p1, ncol=2)) 
dev.off()

} ###

###5.3
if(FALSE){
sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
umap <- sc@reductions$umap@cell.embeddings
df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                  UMAP_2=as.numeric(umap[,2]),
                  MCls=sc@meta.data$MCls, treats=sc@meta.data$treats, Cluster=Idents(sc))

### umap.1
col0 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for(oneMCl in MCls){
   cat(oneMCl, "\n")
   fig0 <- ggplot(df2%>%filter(MCls==oneMCl), aes(x=UMAP_1, y=UMAP_2))+
           geom_point(colour=col0[oneMCl], size=0.1)+
           theme_void()
        
### output
   figfn <- paste("./5_IdenCelltype_output/Figure2.3_", oneMCl, ".umap.png", sep="")
   png(figfn, width=250, height=300, res=150)         
   print(fig0)
   dev.off()       
   Sys.sleep(5)
}

}###


###
### 5.4
if (FALSE){
sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
umap <- sc@reductions$umap@cell.embeddings
df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                  UMAP_2=as.numeric(umap[,2]),
                  MCls=sc@meta.data$MCls, 
                  treats=sc@meta.data$treats, 
                  BATCH=sc@meta.data$BATCH, chem=sc@meta.data$chem,Cluster=Idents(sc))
df2$treat2 <- gsub("-EtOH","", df2$treats)

col0 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

### 1,           
fig1 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(MCls)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(chem),ncol=2)+
        scale_colour_manual(values=col0,guide=guide_legend(override.aes=list(size=2)))+
        #guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA),
              legend.key.size=grid::unit(1,"lines"))
              
png("./5_IdenCelltype_output/Figure2.4_umap.chem.png", width=750, height=450, res=120)
print(fig1)
dev.off()    
       
### 2,
Labalpha <- c("LPS"="a","LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e") 
df2$treat1 <- Labalpha[df2$treat2]
Labtreat <- c("a"="LPS", "b"="LPS+DEX", "c"="PHA", "d"="PHA+DEX", "e"="CTRL")

fig2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(MCls)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(treat1), ncol=3, labeller=as_labeller(Labtreat))+
        scale_colour_manual(values=col0, guide=guide_legend(override.aes=list(size=2)))+
        #guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.85,0.25),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
png("./5_IdenCelltype_output/Figure2.4_umap.treat.png", width=1000, height=900, res=150)
print(fig2)
dev.off() 

### 3,
fig3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(MCls)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(BATCH), ncol=3)+
        scale_colour_manual(values=col0, guide=guide_legend(override.aes=list(size=2)))+
        #guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
png("./5_IdenCelltype_output/Figure2.4_umap.BATCH.png", width=1100, height=900, res=150)
print(fig3)
dev.off()
} 


#x0 <- c("IL7R", "CCR7", "S100A4", "CD8A", "GNLY", "NKG7", "MS4A1", "CD14")
#x0 <- c("CD3D", "CD3E", "CD3G", "GNLY", "MS4A1", "CD14")
######################################
### 6. show marker gene expression ###
######################################

#x0 <- c("IL7R", "CCR7", "S100A4", "CD8A", "GNLY", "NKG7", "MS4A1", "CD14")
## CD4+, IL7R, CCR7, S100A4
## CD8+, CD8A, 
## NK cell, GNLY, NKG7
## B cell, MS4A1
## Monocyte, CD14, LYZ, FCGR3A, MS4A7
## DC, FCER1A, CST3  
#x0 <- c("CD3D", "CD3E", "CD3G", "GNLY", "MS4A1", "CD14")
#x0 <- c("CD3D", "GNLY", "MS4A1", "CD14")

if (FALSE){
sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
DefaultAssay(sc) <- "RNA"
rn <- rownames(sc)

MCls <- c("B cells", "Monocytes", "NK cells", "T cells")
x0 <- c("MS4A1", "CD14", "GNLY", "CD3D")
names(x0) <- MCls

gene0 <- grch38%>%filter(symbol%in%x0)   #%>%select(symbol,ensgene)
ens0 <- paste("S-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(MCls,function(oneMCl){
   ii <- x0[oneMCl]
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           scale_color_gradient(ii,low="lightgrey",high="blue")+
           ggtitle(oneMCl)+
           theme_bw()+
           theme(legend.title=element_text(size=8),
                 legend.key.size=grid::unit(0.7,"lines"),
                 plot.title=element_text(size=12, hjust=0.5))           
   fig0
})

png("./5_IdenCelltype_output/Figure3.0.feature.png", width=750, height=750, res=120)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], figs_ls[[4]], ncol=2)
dev.off()

} ###

###
if(FALSE){
sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
DefaultAssay(sc) <- "RNA"
rn <- rownames(sc)

## figure3.1, Bcell
x0 <- c("MS4A1", "CD79A")
gene0 <- grch38%>%filter(symbol%in%x0)%>%dplyr::select(symbol,ensgene)
ens0 <- paste("S-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
          scale_color_gradient("", low="lightgrey",high="blue")+
           ggtitle(ii)+
           theme_bw()+
           theme(legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5)) 
                      
   fig0
})

png("./5_IdenCelltype_output/Figure3.1.Bcell.feature.png", width=650, height=350, res=100)
plot_grid(figs_ls[[1]], figs_ls[[2]])
dev.off()


###3.2
x0 <- c("CD14", "LYZ", "FCGR3A", "MS4A7","FCER1A", "CST3", "S100A8", "S100A9")
gene0 <- grch38%>%filter(symbol%in%x0)%>%dplyr::select(symbol,ensgene)
ens0 <- paste("S-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           scale_color_gradient("", low="lightgrey",high="blue")+
           ggtitle(ii)+
           theme_bw()+
           theme(axis.title=element_text(size=10),
                 legend.key.size=grid::unit(0.5,"lines"),
                 legend.text=element_text(size=6),
                 plot.title=element_text(size=10, hjust=0.5))           
   fig0
})

png("./5_IdenCelltype_output/Figure3.2.Monocyte.feature.png", width=950, height=500, res=100)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], 
          figs_ls[[4]], figs_ls[[5]], figs_ls[[6]], 
          figs_ls[[7]], figs_ls[[8]], ncol=4)
dev.off()
####

### figure 3.3,  NK cell
x0 <- c("GNLY", "NKG7")
gene0 <- grch38%>%filter(symbol%in%x0)%>%dplyr::select(symbol,ensgene)
ens0 <- paste("S-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           scale_color_gradient("", low="lightgrey",high="blue")+
           ggtitle(ii)+
           theme_bw()+
           theme(legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5)) 
                      
   fig0
})

png("./5_IdenCelltype_output/Figure3.3.NKcell.feature.png", width=650, height=350, res=100)
plot_grid(figs_ls[[1]], figs_ls[[2]], ncol=2)
dev.off()

### figure 3.4,
x0 <- c("CD3D", "CCR10", "TNFRSF18", "IL7R", "CCR7", "S100A4", "CD8A", "ID3")
gene0 <- grch38%>%filter(symbol%in%x0,chr%in%as.character(1:22))%>%dplyr::select(symbol,ensgene)
ens0 <- paste("S-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           scale_color_gradient("", low="lightgrey",high="blue")+
           ggtitle(ii)+
           theme_bw()+
           theme(axis.title=element_text(size=10),
                 legend.key.size=grid::unit(0.5,"lines"),
                 legend.text=element_text(size=6),
                 plot.title=element_text(size=10, hjust=0.5))                       
   fig0
})

png("./5_IdenCelltype_output/Figure3.4.Tcell.feature.png", width=950, height=500, res=100)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], figs_ls[[4]],
          figs_ls[[5]], figs_ls[[6]], figs_ls[[7]], figs_ls[[8]], ncol=4)
dev.off()
} ###
##

###
### (3)
if (TRUE){

sc <-  read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
DefaultAssay(sc) <- "RNA"
rn <- rownames(sc)

### marker gene information
x0 <- c("MS4A1", "CD79A", "CD14", "MS4A7", "GNLY", "NKG7", "CD3D", "CD8A")
MCls <- rep(c("B cells", "Monocytes", "NK cells", "T cells"), each=2)
gene0 <- grch38%>%filter(symbol%in%x0)
ens <- sapply(x0, function(ii){
  ens_ii <- gene0%>%filter(symbol==ii)%>%dplyr::pull(ensgene)
  ens_ii <- paste("S-", ens_ii, sep="")
  rn0 <- rn[grepl(ens_ii, rn)]
  rn0
}) 
anno <- data.frame(MCls=MCls, symbol=x0, ens=ens)

#### plot
figs_ls <- lapply(1:nrow(anno), function(i){
   oneMCl <- anno$MCls[i]
   ens_ii <- anno$ens[i]
   symbol <- anno$symbol[i]
   fig0 <- FeaturePlot(sc, features=ens_ii)+
           scale_color_gradient(symbol,low="lightgrey",high="blue")+
           ggtitle(oneMCl)+
           theme_bw()+
           theme(legend.title=element_text(size=8),
                 legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5))           
   fig0
})

png("./5_IdenCelltype_output/Figure4_MCls.feature.png", width=950, height=500, res=100)
print(plot_grid(figs_ls[[1]], figs_ls[[2]], 
          figs_ls[[3]], figs_ls[[4]],
          figs_ls[[5]], figs_ls[[6]],
          figs_ls[[7]], figs_ls[[8]], nrow=2, byrow=F))
dev.off()

}###


###
###
#sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.New.rds")
#meta <- sc@meta.data
#my <- meta%>%
#      group_by(treats, MCls)%>%
#       summarize(Freq=n())%>%
#       group_by(treats)%>%
#       mutate(Perc=Freq/sum(Freq)*100)
#write.csv(my, file="./5_IdenCelltype_output/5.MCls.summary.csv", row.names=F)

#######################
### unspliced genes ###
#######################

if(FALSE){
### figure 3.1, CD4+
x0 <- c("IL7R", "CCR7", "S100A4", "CD8A")
gene0 <- grch38%>%filter(symbol%in%x0)%>%select(symbol,ensgene)
ens0 <- paste("U-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           ggtitle(ii)           
   fig0
})

png("./5_IdenCelltype_output/Figure5.1_U.Tcell.png", width=700, height=600, res=90)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], figs_ls[[4]], ncol=2)
dev.off()


### figure 3.2,  NK cell
x0 <- c("GNLY", "NKG7")
gene0 <- grch38%>%filter(symbol%in%x0)%>%select(symbol,ensgene)
ens0 <- paste("U-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           ggtitle(ii)           
   fig0
})

png("./5_IdenCelltype_output/Figure5.2_U.NKcell.png", width=600, height=400, res=90)
plot_grid(figs_ls[[1]], figs_ls[[2]], ncol=2)
dev.off()

## figure3.3
x0 <- c("MS4A1")
gene0 <- grch38%>%filter(symbol%in%x0)%>%select(symbol,ensgene)
ens0 <- paste("U-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           ggtitle(ii)           
   fig0
})

png("./5_IdenCelltype_output/Figure5.3_U.Bcell.feature.png", width=300, height=400, res=60)
plot_grid(figs_ls[[1]])
dev.off()

###3.4
x0 <- c("CD14", "LYZ", "FCGR3A", "MS4A7","FCER1A", "CST3")
gene0 <- grch38%>%filter(symbol%in%x0)%>%select(symbol,ensgene)
ens0 <- paste("U-", gene0$ensgene, sep="")
ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
names(ens1) <- gene0$symbol

figs_ls <- lapply(x0,function(ii){
   fig0 <- FeaturePlot(sc, features=ens1[ii])+
           ggtitle(ii)           
   fig0
})

png("./5_IdenCelltype_output/Figure5.4_U.Monocyte.png", width=800, height=600, res=80)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], 
          figs_ls[[4]], figs_ls[[5]], figs_ls[[6]], ncol=3)
dev.off()
} ###



#############################
### compare resulst rpca  ###
#############################
##dir.create("./5_IdenCelltype_output/SCT2/")
####
#if(FALSE){
#meta <- read_rds("./5_IdenCelltype_output/3_Meta.Zheng68k.rds")
#dd <- meta[,c("predicted.id","prediction.score.max","NEW_BARCODE")]
#### figure 2, umap by identified cell type
#sc <- read_rds("./4_SCTransform_output/2_SCT2.Chem.dims50.Cl.rds")
#umap <- Embeddings(sc, reduction="umap")
#
#mydf <- data.frame(NEW_BARCODE=as.character(rownames(meta)),
#                   UMAP_1=as.numeric(umap[,1]), 
#                   UMAP_2=as.numeric(umap[,2]),
#                   BATCH=sc$BATCH, 
#                   chem=sc$chem,
#                   treats=sc$treats,
#                   cluster=sc$seurat_clusters,
#                   MCls=dd$predicted.id)
#
#fig2 <- ggplot(mydf,aes(x=UMAP_1,y=UMAP_2, colour=MCls))+
#        geom_point(size=0.1)+
#        guides(col=guide_legend(override.aes=list(size=3)))+
#        scale_colour_brewer(palette="Set3")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA))
#              
#figfn <- "./5_IdenCelltype_output/SCT2/Figure2.2.umap.Zheng68k.png"
#png(figfn, width=700, height=600, res=110)
#fig2
#dev.off()
#
#
### figure 3, heatmap,
#df0 <- mydf%>%group_by(MCls, cluster)%>%
#       summarize(Freq=n())%>%
#       group_by(cluster)%>%
#       mutate(Perc=Freq/sum(Freq)*100)
#       
#fig3 <- ggplot(df0, aes(x=cluster, y=MCls,fill=Perc))+
#            geom_tile()+
#            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
#            xlab("")+ylab("")+
#            theme_bw()+theme(axis.text.x=element_text(hjust=1))
#figfn <- "./5_IdenCelltype_output/SCT2/Figure2.3.heatmap.Zheng68k.png"
#png(figfn, width=800, height=600, res=120)
#fig3
#dev.off()
#
#
####figure 4
#sc <- read_rds("./4_SCTransform_output/2_SCT2.Chem.dims50.Cl.rds")
##Idents(sc) <- sc$integrated_snn_res.0.15
#meta <- sc@meta.data
#
#ii <- Idents(sc)
#MCls <- rep("Tcell", nrow(meta))
#MCls[ii==2] <- "NKcell"
#MCls[ii==3] <- "Bcell"
#MCls[ii%in%c(4,11)] <- "Monocyte"
#
#meta$MCls <- MCls
#sc <- AddMetaData(sc,metadata=meta)
#
#opfn <- "./5_IdenCelltype_output/SCT2/4_SCAIP.MCls.New.rds"
#write_rds(sc,opfn)
#
##sc <-  read_rds("./5_IdenCelltype_output/SCT2/4_SCAIP.MCls.New.rds")
#col0 <- c("Tcell"="#e41a1c", "NKcell"="#377eb8", 
#          "Bcell"="#4daf4a", "Monocyte"="#984ea3")
#fig4 <- DimPlot(sc, label=F, group.by="MCls", cols=col0)+
#        theme(legend.position=c(0.1, 0.85),
#              legend.text=element_text(size=12),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA),            
#              panel.border=element_rect(colour="black", fill=NA))
#png("./5_IdenCelltype_output/SCT2/Figure2.4.umap.MCls.png", width=500, height=600, res=100)
#fig4
#dev.off()
#
#### 5, Expression of marker gene 
###
#sc <-  read_rds("./5_IdenCelltype_output/SCT2/4_SCAIP.MCls.New.rds")
#DefaultAssay(sc) <- "RNA"
#rn <- rownames(sc)
#x0 <- c("CD3D", "GNLY", "MS4A1", "CD14")
#gene0 <- grch38%>%filter(symbol%in%x0)#%>%select(symbol,ensgene)
#ens0 <- paste("S-", gene0$ensgene, sep="")
#ens1 <- sapply(ens0, function(ii) rn[grepl(ii,rn)])
#names(ens1) <- gene0$symbol
#
#MCls <- c("T cells", "NK cells", "B cells", "Monocytes")
#names(MCls) <- x0
#figs_ls <- lapply(x0,function(ii){
#   fig0 <- FeaturePlot(sc, features=ens1[ii])+
#           scale_color_gradient(ii,low="lightgrey",high="blue")+
#           ggtitle(MCls[ii])          
#   fig0
#})
#
####
#png("./5_IdenCelltype_output/SCT2/Figure2.5.umap.feature.png", width=1000, height=900, res=120)
#plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], figs_ls[[4]], ncol=2)
#dev.off()
#}
