
rm(list=ls())
source("./Bin/LibraryPackage.R")
outdir <- "./4_Harmony_output/"
if (!file.exists(outdir)) dir.create(outdir)


##################################################
### 1, Integrate data by chemistry (V2 and V3) ###
##################################################
### setting parallel 

future::plan(strategy = 'multicore', workers = 5)
options(future.globals.maxSize = 20 * 1024 ^ 3)
plan()

sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method="vst", nfeatures=2000)
sc <- ScaleData(sc, verbose=TRUE)
sc <- RunPCA(sc, npcs=100, verbose=TRUE)
sc <- RunHarmony(sc, "chem", reduction="pca")
###
opfn <- "./4_Harmony_output/1_Norm.Chem.pca.rds"
write_rds(sc, opfn)


sc <- RunUMAP(sc, dims=1:50, reduction="harmony", verbose=TRUE)
sc <- FindNeighbors(sc, dims=1:50, reduction="harmony", verbose=TRUE)
sc <- FindClusters(sc, resolution=0.15, verbose=TRUE)        ##we used 0.15
###
write_rds(sc,"./4_Harmony_output/2_Norm.Chem.dims50.Cl.rds")

### 5 pcs for 53%; 9 pcs for 60%
### 19 pcs for 70.3%; 20 pcs for 71%
### 37 pcs for 80.00%; 50 pcs for 85%
### 65 pcs for 90%;
#sc <- read_rds("./4_Harmony_output/1_Norm.Chem.pca.rds")
#lapply(c(20,37,65), function(ii){
#   sci <- RunUMAP(sc,dims=1:ii, reduction="harmony", verbose=TRUE)
#   
#})



##option (2)
#rm(list=ls())
#sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
#sc <- NormalizeData(sc)
#sc <- FindVariableFeatures(sc, selection.method="vst", nfeatures=2000)
#sc <- ScaleData(sc, verbose=TRUE)
#sc <- RunPCA(sc, npcs=100, verbose=TRUE)
#sc <- RunHarmony(sc, c("chem","treats"), reduction="pca")
###
#opfn <- "./4_Harmony_output/1_Treats.Chem.pca.rds"
#write_rds(sc, opfn)

###
#sc <- read_rds("./4_Harmony_output/1_Treats.chem.pca.rds")
#sc <- RunUMAP(sc, dims=1:50, reduction="harmony", verbose=TRUE)
#sc <- FindNeighbors(sc, dims=1:50, reduction="harmony", verbose=TRUE)
#sc <- FindClusters(sc, resolution=0.15, verbose=TRUE)        ##we used 0.15
###
#write_rds(sc,"./4_Harmony_output/2_Treats.Chem.dims50.Cl.rds")

###option (3)
sc <- read_rds("./4_Harmony_output/1_Norm.Chem.pca.rds")
sc <- RunUMAP(sc, dims=1:50, reduction="harmony", verbose=TRUE)
sc <- FindNeighbors(sc, dims=1:50, reduction="harmony", verbose=TRUE)
sc <- FindClusters(sc, resolution=0.5, verbose=TRUE)        ##we used 0.15
###
write_rds(sc,"./4_Harmony_output/2_Norm.Chem.dims50.res0.5.rds")


#######################
### 2, show figures ###
#######################

### (1) 
sc <- read_rds("./4_Harmony_output/1_Norm.Chem.pca.rds")

##1.0, pca
y1 <- Stdev(sc, reduction="pca")
y2 <- cumsum(y1^2)
y2 <- y2/sum(y1^2) 

dd <- data.frame(x=1:length(y1), y1=y1, y2=y2)
fig0 <- ggplot(dd)+
        geom_point(aes(x=x,y=y1), color="red")+
        scale_y_continuous("Standard Deviation")+
        theme_bw()
fig1 <- ggplot(dd)+
        geom_point(aes(x=x,y=y2), color="blue")+
        scale_y_continuous("Cumulative proportion")+
        theme_bw()
figfn0 <- "./4_Harmony_output/Figure1.0.elbow.png"        
png(figfn0, width=800, height=500, res=120)
plot_grid(fig0, fig1, labels=c("A","B"))
dev.off()

###1.1, harmony
### 5 pcs for 53%; 9 pcs for 60%
### 19 pcs for 70.3%; 20 pcs for 71%
### 37 pcs for 80.00%; 50 pcs for 85%
### 65 pcs for 90%;

y1 <- Stdev(sc, reduction="harmony")
y2 <- cumsum(y1^2)
y2 <- y2/sum(y1^2) 

dd <- data.frame(x=1:length(y1), y1=y1, y2=y2)
fig0 <- ggplot(dd)+
        geom_point(aes(x=x,y=y1), color="red")+
        scale_y_continuous("Standard Deviation")+
        theme_bw()
fig1 <- ggplot(dd)+
        geom_point(aes(x=x,y=y2), color="blue")+
        scale_y_continuous("Cumulative proportion",limits=c(0,1))+
        theme_bw()
figfn0 <- "./4_Harmony_output/Figure1.1.elbow.png"        
png(figfn0, width=800, height=500, res=120)
plot_grid(fig0, fig1, labels=c("A","B"))
dev.off()
}##

###

####1.2, comparison of umap of different pcs
#fig_ls <- lapply(c(30, 40, 50, 67), function(ii){
#   fn <- paste("./4_OldNormAndChem_output/3_OldNorm.Chem.pc", ii, ".rds", sep="")
#   cat(ii,"\n")
#   sc0 <- read_rds(fn)
#   fig0 <- DimPlot(sc0, label=TRUE)+
#           NoLegend()+
#           ggtitle(paste(ii, "PCs", sep=""))+
#           theme(plot.title=element_text(hjust=0.5),
#                 panel.border=element_rect(fill=NA,colour="black"))
#   fig0
#})
#png("./4_Harmony_output/Figure1.2.umap.png", width=1200, height=1200, res=150)
#plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], labels="AUTO")
dev.off()   
   


####################################################
### (2), show umap group by different conditions ###
####################################################
if(FALSE){
sc <- read_rds("./4_Harmony_output/2_Norm.Chem.dims50.Cl.rds")
umap <- sc@reductions$umap@cell.embeddings
df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                  UMAP_2=as.numeric(umap[,2]), 
                  treats=sc$treats,
                  BATCH=sc$BATCH, chem=sc$chem,
                  cluster=sc$seurat_clusters) 
                  #cluster1=sc$integrated_snn_res.0.1,
                  #cluster2=sc$integrated_snn_res.0.2)
df2$treat2 <- gsub("-EtOH","", df2$treats)

### 1,           
fig1 <- ggplot(df2, aes(x=UMAP_1,y=UMAP_2, colour=factor(cluster)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(chem),ncol=2)+
        guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
        theme_bw()+
        theme(legend.title=element_blank(), legend.position="none",
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
png("./4_Harmony_output/Figure2.1.umap.chem.png", width=700, height=450, res=100)
fig1
dev.off()    
       
### 2,
Labalpha <- c("LPS"="a","LPS-DEX"="b", "PHA"=c, "PHA-DEX"="d", "CTRL"="e") 
df2$treat1 <- Labalpha[df2$treat2]
Labtreat <- c("a"="LPS", "b"="LPS+DEX", "c"="PHA", "d"="PHA+DEX", "e"="CTRL")

fig2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(cluster)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(treat1), ncol=2, labeller=as_labeller(Labtreat))+
        guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.8,0.2),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
png("./4_Harmony_output/Figure2.2.umap.treat.png", width=700, height=900, res=140)
fig2
dev.off() 

### 3,
fig3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(cluster)))+
        geom_point(size=0.1)+
        facet_wrap(~factor(BATCH),ncol=3)+
        guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA))
              
png("./4_Harmony_output/Figure2.3.umap.BATCH.png", width=1000, height=800, res=140)
fig3
dev.off() 

### UMAP by cluster
sc <- read_rds("./4_Harmony_output/2_Norm.Chem.dims50.res0.5.rds")
fig4 <- DimPlot(sc, label=TRUE)+
   NoLegend()+
   theme(plot.title=element_text(hjust=0.5),
         panel.border=element_rect(fill=NA, colour="black"))              
png("./4_Harmony_output/Figure3.4_umap.Cluster.png", width=500, height=600, res=120)
fig4
dev.off() 


umap <- sc@reductions$umap@cell.embeddings
df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                  UMAP_2=as.numeric(umap[,2]), 
                  cluster=sc$RNA_snn_res.0.2)

### 1,           
fig1 <- ggplot(df2, aes(x=UMAP_1,y=UMAP_2, colour=factor(cluster)))+
   geom_point(size=0.1)+
   guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
   theme_bw()+
   theme(legend.title=element_blank(), legend.position="none",
         legend.background=element_rect(colour=NA, fill=NA),
         legend.key=element_rect(fill=NA))
              
png("./4_Harmony_output/Figure3.5_umap.cluster.png", width=500, height=600, res=120)
fig1
dev.off()    

#######################################
### (3), harmony by chem and treats ###
#######################################

#sc <- read_rds("./4_Harmony_output/2_Treats.Chem.dims50.Cl.rds")
#umap <- sc@reductions$umap@cell.embeddings
#df2 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
#                  UMAP_2=as.numeric(umap[,2]), 
#                  treats=sc$treats,
#                  BATCH=sc$BATCH, chem=sc$chem,
#                  cluster=sc$seurat_clusters) 
                  #cluster1=sc$integrated_snn_res.0.1,
                  #cluster2=sc$integrated_snn_res.0.2)
#df2$treat2 <- gsub("-EtOH","", df2$treats)

### 1,           
#fig1 <- ggplot(df2, aes(x=UMAP_1,y=UMAP_2, colour=factor(cluster)))+
#        geom_point(size=0.1)+
#        facet_wrap(~factor(chem),ncol=2)+
#        guides(col=guide_legend(override.aes=list(size=3),ncol=2))+
#        theme_bw()+
#        theme(legend.title=element_blank(), legend.position="none",
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA))
              
#png("./4_Harmony_output/Figure4.1.umap.chem.png", width=700, height=450, res=100)
#fig1
#dev.off()    
       
### 2,
#x1 <- as.character(df2$treat2)
#x1[x1=="CTRL"] <- "e" 
#x1[x1=="LPS"] <- "a"
#x1[x1=="LPS-DEX"] <- "b"
#x1[x1=="PHA"] <- "c"
#x1[x1=="PHA-DEX"] <- "d"

#df2$treat1 <- x1
#treat.lab <- c("CTRL", "LPS", "LPS-DEX", "PHA", "PHA-DEX")
#names(treat.lab) <- c("e", "a", "b", "c", "d")

#fig2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(cluster)))+
#        geom_point(size=0.1)+
#        facet_wrap(~factor(treat1), ncol=2, labeller=as_labeller(treat.lab))+
#        guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.position=c(0.8,0.2),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA))
#              
#png("./4_Harmony_output/Figure4.2.umap.treat.png", width=1000, height=1000, res=140)
#fig2
#dev.off() 

### 3,
#fig3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(cluster)))+
#        geom_point(size=0.1)+
#        facet_wrap(~factor(BATCH),ncol=3)+
#        guides(col=guide_legend(override.aes=list(size=2),ncol=2))+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA))
#              
#png("./4_Harmony_output/Figure4.3.umap.BATCH.png", width=1000, height=800, res=140)
#fig3
#dev.off() 







#sc2 <- sc
#Idents(sc2) <- sc2$integrated_snn_res.0.5
#fig4 <- DimPlot(sc2, label=TRUE)+
#           NoLegend()+
#           theme(plot.title=element_text(hjust=0.5),
#                 panel.border=element_rect(fill=NA,colour="black"))              
#png("./4_OldNormAndChem_output/Figure2.4.umap.Cluster.res0.5.png", width=500, height=600, res=120)
#fig4
#dev.off()
