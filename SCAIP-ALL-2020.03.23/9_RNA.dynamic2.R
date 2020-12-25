#
rm(list=ls())
source("./Bin/LibraryPackage.R")

outdir <- "./9_RNA.dynamic2_output/Old/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F)

#future::plan(strategy = 'multicore', workers = 5)
#options(future.globals.maxSize = 20 * 1024 ^ 3)
#plan()


###########################################################
### 2020-12-11, psedotime definition based on 6,571 DEGs ###
###          By Julong wei                              ###
###########################################################
## In final we use this results


#########################################################
### 1, New SCAIP object using 6,571 DEG, spliced gene  ###
#########################################################

if(FALSE){   
cat("1.", "Extract DEG and Batch 1,4,5 and 6", "\n")
#####################################
### (1), extract 6413 spliced DEG ###
#####################################
## single cell data
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
counts <- sc@assays$RNA@counts
meta <- sc@meta.data 
DefaultAssay(sc) <- "RNA"

### differetial expressed genes
load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
DEG <- sigs

### 
anno <- data.frame(rn=rownames(sc))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), uns=grepl("S-", rn))%>%
        mutate(sel=ensgene%in%DEG)%>%
        filter(sel,uns)

### Generating new SCAIP object, only containing DEG          
count0 <- counts[anno$rn,]
sc2 <- CreateSeuratObject(count0, project="SCAIP")
sc2@meta.data <- meta

### output
opfn0 <- "./9_RNA.dynamic2_output/1_SCAIP.DEG.rds"
write_rds(sc2, opfn0)
 

##################################################
### (2), keeping Batch 1,4,5,6  ###
##################################################
sc2 <- read_rds("./9_RNA.dynamic2_output/1_SCAIP.DEG.rds")
#rn <- rownames(sc2)

###remain Batch 1, 4 and 5
meta <- sc2@meta.data
cellNew <- meta[meta$BATCH %in% c("SCAIP1", "SCAIP4", "SCAIP5", "SCAIP6"),"NEW_BARCODE"]
scNew <- subset(sc2, cells=cellNew)
###
opfn1 <- "./9_RNA.dynamic2_output/2_Batch1456.DEG.rds"
write_rds(scNew, opfn1)

} ##End, 1



#########################################################
### (3), split by cell type and correct batch effects ###
#########################################################
if(FALSE){
cat("2.", "split data then cluster analysis", "\n") 
rm(list=ls())

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

sc2 <- read_rds("./9_RNA.dynamic2_output/2_Batch1456.DEG.rds")
rn <- as.character(rownames(sc2))
sp <- SplitObject(sc2, split.by = "MCls")


correct <- "Old"  ##defaul Old
outdir2 <- paste("./9_RNA.dynamic2_output/", correct, "/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F)

for (oneMCl in MCls){
   cat(oneMCl, "\n")
   sc0 <- sp[[oneMCl]]
   if (correct=="Old"){
      sp0 <- SplitObject(sc0, split.by = "chem")
      sp0 <- lapply(X=sp0, function(x){
         x <- NormalizeData(x)
         x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      })
      anchors <- FindIntegrationAnchors(object.list=sp0, dims=1:50)
      combined <- IntegrateData(anchorset=anchors, dims=1:50)  
  
      sc0 <- ScaleData(combined, features=rn, verbose=T)
      sc0 <- RunPCA(sc0, features=rn, verbose=T, npcs=100)
      sc0 <- RunUMAP(sc0, dims=1:50, verbose=T)
      opfn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="") 
      write_rds(sc0, opfn)
   }else{
      sc0 <- NormalizeData(sc0)
      sc0 <- FindVariableFeatures(sc0, selection.method="vst", nfeatures=2000)
      sc0 <- ScaleData(sc0, features=rn, verbose=T)
      sc0 <- RunPCA(sc0, features=rn, npcs=100, verbose=T)
      sc0 <- RunHarmony(sc0, "chem", reduction="pca")     
      sc0 <- RunUMAP(sc0, dims=1:50, reduction="harmony", verbose=T)
      opfn <- paste("./9_RNA.dynamic2_output/Harmony/3_MCl.", oneMCl, ".harmony.rds", sep="")
      write_rds(sc0, opfn) 
      cat("Harmony Done", "\n")  
   }
   
}###

} ### End, 2


#######################################################
### 3, show the results, PCA and UMAP as pseudotime ###
####################################################### 
if(FALSE){

rm(list=ls())
cat("3.", "show PCA and UMAP as psedotime", "\n")
################
### (1). pca ###
################  
cat("(1).", "PCA", "\n")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
df2 <- map_dfr(MCls, function(oneMCl){
   #cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
   sc0 <- read_rds(fn)
   meta <- sc0@meta.data
   pca <- Embeddings(sc0, reduction="pca")          
   dd <- data.frame(pca[,1:2], treats=meta$treats, MCls=oneMCl, chem=meta$chem)
   dd
})
df2$treat2 <- gsub("-EtOH", "", df2$treats)

fig1 <- ggplot(df2, aes(x=PC_1, y=PC_2))+
        geom_point(aes(colour=factor(treat2)), size=0.1)+
        facet_wrap(~factor(MCls), ncol=2, scales="free")+
        scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=3)))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA),
              strip.background=element_blank(),
              strip.text.x=element_text(size=12))
              
png("./9_RNA.dynamic2_output/Old/Figure1.1.pca.png", width=700, height=600, res=120)
print(fig1)
dev.off() 
Sys.sleep(5)

###
fig1.2 <- ggplot(df2, aes(x=PC_1, y=PC_2))+
        geom_point(colour="#fb9a99", size=0.1)+
        facet_grid(factor(chem)~factor(MCls),scales="free_x")+
        theme_bw()
              
png("./9_RNA.dynamic2_output/Old/Figure1.2.pca.chem.png", width=700, height=600, res=120)
print(fig1.2)
dev.off() 


#################
### (2). umap ###
#################
cat("(2).", "umap", "\n")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
df2 <- map_dfr(MCls, function(oneMCl){
###
   #cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
   sc0 <- read_rds(fn)
   meta <- sc0@meta.data
   umap <- Embeddings(sc0, reduction="umap")          
   dd <- data.frame(umap[,1:2], treats=meta$treats, MCls=oneMCl, chem=meta$chem)
   dd
})
df2$treat2 <- gsub("-EtOH", "", df2$treats)

fig2.1 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2))+
        geom_point(aes(colour=factor(treat2)), size=0.1)+
        facet_wrap(~factor(MCls), ncol=2, scales="free")+
        scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=3)))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA),
              strip.background=element_blank(),
              strip.text.x=element_text(size=12))
              
png("./9_RNA.dynamic2_output/Old/Figure2.1.UMAP.png", width=700, height=600, res=120)
print(fig2.1)
dev.off() 
Sys.sleep(5)

fig2.2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2))+
        geom_point(colour="#fb9a99", size=0.1)+
        facet_grid(factor(chem)~factor(MCls),scales="free_x")+
        theme_bw()
        #scale_colour_manual(values=col0, guide=guide_legend(override.aes=list(size=3)))+
        #guides(colour=guide_legend(override.aes=list(size=3)))+
        #theme_bw()
        #theme(legend.title=element_blank(),
        #      legend.background=element_rect(colour=NA, fill=NA),
        #      legend.key=element_rect(fill=NA))
              
png("./9_RNA.dynamic2_output/Old/Figure2.2.UMAP.chem.png", width=700, height=600, res=120)
print(fig2.2)
dev.off() 
Sys.sleep(5)

} ###

######################################################
### (3), density UMAP plots cell type by condition ###
######################################################
if(FALSE){
rm(list=ls())
cat("(3).", "Show density UMAP plots, facet_grid(MCls~treats)", "\n")
treat2lab <- c("CTRL"="CTRL",
               "LPS"="LPS", "LPS-DEX"="LPS+DEX", 
               "PHA"="PHA", "PHA-DEX"="PHA+DEX")

MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell")
## data for UMAP plot
df2 <- map_dfr(MCls, function(oneMCl){           
   ##
   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
   sc0 <- read_rds(fn)
   umap <- sc0@reductions$umap@cell.embeddings
   df0 <- data.frame(UMAP_1=as.numeric(umap[,1]), 
                     UMAP_2=as.numeric(umap[,2]), 
                     treats=sc0@meta.data$treats, MCls=oneMCl)
   df0
})

df2$treat2 <- gsub("-EtOH", "", df2$treats)

fig0 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2))+      
           stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
           facet_grid(MCls~treat2, scales="free_y", labeller=labeller(treat2=treat2lab))+
           scale_fill_viridis_c(direction=-1)+
           theme_bw()
figfn <- "./9_RNA.dynamic2_output/Old/Figure3.1.umapDensity.png"
png(figfn, width=900, height=800, res=150)
print(fig0)
dev.off()
}

####################
### (3), harmony ###
####################
#cat("3.3.", "harmony", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
#df2 <- lapply(MCls, function(oneMCl){
###
   #cat(oneMCl,"\n")
#   fn <- paste("./9_RNA.dynamic2_output/Harmony/3_MCl.", oneMCl, ".harmony.rds", sep="")
#   sc0 <- read_rds(fn)
#   meta <- sc0@meta.data
#   X <- Embeddings(sc0, reduction="harmony")          
#   dd <- data.frame(X[,1:2], treats=meta$treats, MCls=oneMCl)
#   dd
#})
#df2 <- do.call(rbind,df2)
#df2$treat2 <- gsub("-EtOH", "", df2$treats)

#fig2.3 <- ggplot(df2, aes(x=harmony_1, y=harmony_2))+
#        geom_point(aes(colour=factor(treat2)), size=0.1)+
#        facet_wrap(~factor(MCls),ncol=2,scales="free")+
#        scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=3)))+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA),
#              strip.background=element_blank(),
#              strip.text.x=element_text(size=12))
              
#png("./9_RNA.dynamic2_output/Harmony/Figure2.3.harmony.png", width=700, height=600, res=120)
#print(fig2.3)
#dev.off() 
#Sys.sleep(5)
 

#######################################################################
### 4, dimensional reduction using Linear Determinat Analysis (LDA) ###
#######################################################################

##########################
### 4.1, calculate LDA ###
##########################

if (FALSE){
rm(list=ls())
cat("4.1.", "Calculate LDA based on gene expression", "\n")

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
tmp <- mclapply(MCls, function(oneMCl){

   s1 <- Sys.time()
   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
   sc <- read_rds(fn)

   X <- sc@assays$RNA@data 
   ### here we directly used counts data. Maybe we need to try to normalize data first. 
   #X <- t(Embeddings(sc, reduction="pca"))
   meta <- sc@meta.data
   treats <- unique(meta$treats)

   ### Sw, within classes scatter matrix 
   xx_ls <- lapply(treats, function(one){
      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
      Xi <- X[,ii]
      xc <- scale(t(Xi),scale=FALSE)
      xx <- crossprod(xc)
      return(xx)
   })
   Sw <- Reduce("+", xx_ls)
   cat("Within", 1, "\n")
   
   ### Sb, between classes scatter matrix 
   mu <- apply(X, 1, mean)
   xx_ls <- lapply(treats, function(one){
      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
      Xi <- X[,ii]
      ni <- ncol(Xi)
      mu_k <- apply(Xi, 1, mean)
      xk <- mu_k-mu
      xx <- tcrossprod(xk)*ni
      return(xx)   
   })
   Sb <- Reduce("+",xx_ls)
   cat("Between", 2, "\n")

   ### pseudo inverse of matrix, Sw
   sw.eigen <- eigen(Sw)
   subi <- sw.eigen$values>1e-02
   uw <-sw.eigen$vectors[,subi]
   dw <- sw.eigen$values[subi]
   sw2 <- uw %*% diag(1/dw) %*% t(uw)

   ### generazied eigen decomposition
   sb.eigen <- eigen(Sb,symmetric=T)
   ub <- sb.eigen$vectors[,1:4]
   db <-sb.eigen$values[1:4]
   sb2 <- ub %*% diag(sqrt(db)) %*% t(ub)
   sb2i <- ub %*% diag(1/sqrt(db)) %*% t(ub)
   
   ss <- sb2 %*% sw2 %*% sb2
   ss.eigen <- eigen(ss,symmetric=T) 
   ###  
   w <- sb2i %*% ss.eigen$vectors[,1:4]
   values <- ss.eigen$values
   
   ### the top 4 LDA
   zk <- crossprod(X,w)
   meta$LDA_1 <- zk[,1]
   meta$LDA_2 <- zk[,2]
   meta$LDA_3 <- zk[,3]
   meta$LDA_4 <- zk[,4]
   opfn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
   write_rds(meta, opfn)
   
   ### loading information, weights of genes for the top 4 LDA
   w <- as.matrix(w)
   rownames(w) <- gsub("S-", "", rownames(X))
   colnames(w) <- paste("LDA", 1:4, sep="_")
   opfn <- paste("./9_RNA.dynamic2_output/Old/4.2_LDA.", oneMCl, ".loading.rds", sep="")
   write_rds(w, opfn)
      
   s2 <- Sys.time()
   d12 <- difftime(s2,s1,units="mins")
   cat(oneMCl, ":", d12, "\n")
   return(values)
},mc.cores=1)
save(tmp, file="./9_RNA.dynamic2_output/Old/4.3_LDA.eigenvalue.RData")

} ##4.1, End 


########################
### 4.2 show figures ###
########################

if (FALSE){
cat("4.2.", "Show results of LDA", "\n")

load("./9_RNA.dynamic2_output/Old/4.3_LDA.eigenvalue.RData")
tmp <- do.call(cbind,tmp)
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

### (1). show proportion of LDA
cat("(1).", "Show proportion of LDA", "\n")
x <- sweep(tmp, 2, colSums(tmp),"/")
colnames(x) <- MCls

###
fig_ls <- lapply(MCls, function(ii){
   xi <- x[,ii]
   di <- data.frame(x=1:10, y=xi[1:10])
   p <- ggplot(di,aes(x=x,y=y))+
        geom_point()+xlab("Nums of LDA")+ylab("Proportion")+
        ggtitle(ii)+
        theme_bw()+
        theme(plot.title=element_text(hjust=0.5))
   p
})

###
figfn <- paste("./9_RNA.dynamic2_output/Old/Figure4.1.decay.png", sep="")
png(figfn,width=750, height=700, res=120)
print(plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], ncol=2))
          #labels="AUTO", label_fontface="plain", label_fontfamily="serif", ncol=2)
dev.off() 


### (2). Show LDA1 and LDA2 ### 
cat("(2).", "Show LDA plots", "\n")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
           
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- map_dfr(MCls, function(oneMCl){
   #cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
   meta <- read_rds(fn)      
   dd <- data.frame(LDA_1=meta$LDA_1, LDA_2=meta$LDA_2, treats=meta$treats, MCls=oneMCl, chem=meta$chem)
   dd
})

df2$treat2 <- gsub("-EtOH","",df2$treats)
lab2 <- c("CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fig3.2 <- ggplot(df2, aes(x=LDA_1, y=LDA_2))+
        geom_point(aes(colour=factor(treat2)), size=0.1)+
        facet_wrap(~factor(MCls),ncol=2,scales="free")+
        scale_colour_manual(values=col1, labels=lab2, guide=guide_legend(override.aes=list(size=2)))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_rect(colour=NA, fill=NA),
              legend.key=element_rect(fill=NA),
              strip.background=element_blank(),
              strip.text.x=element_text(size=12))
              
png("./9_RNA.dynamic2_output/Old/Figure4.2.LDA12.png", width=800, height=650, res=120)
print(fig3.2)
dev.off() 

###
cat("(3).", "Show LDA plots by chem", "\n") 
fig3.3 <- ggplot(df2, aes(x=LDA_1, y=LDA_2))+
        geom_point(colour="#fb9a99", size=0.1)+
        facet_grid(factor(chem)~factor(MCls),scales="free_x")+
        theme_bw()
              
png("./9_RNA.dynamic2_output/Old/Figure4.3.LDA12.chem.png", width=700, height=600, res=120)
print(fig3.3)
dev.off() 

Sys.sleep(5)

} ###


#########################################################
### 4.3, show LDA_1 and LDA_2 by cell type separately ###
######################################################### ***
if(FALSE){
cat("4.3", "Show LDA plots with density plots", "\n")

cat("(1).", "By cell type together", "\n")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
           
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- map_dfr(MCls, function(oneMCl){
##
   cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
   meta <- read_rds(fn)      
   dd <- data.frame(LDA_1=meta$LDA_1, LDA_2=meta$LDA_2, treats=meta$treats, MCls=oneMCl)
   dd
})

df2$treat2 <- gsub("-EtOH", "", df2$treats)

## Bcell
fig_ls <- lapply(MCls, function(oneMCl){
   df0 <- df2%>%filter(MCls==oneMCl)
   #if( oneMCl=="NKcell") df0$LDA_2 <- -df0$LDA_2
   fig0 <- ggplot(df0, aes(x=LDA_1, y=LDA_2))+
           geom_point(aes(colour=factor(treat2)), size=0.1)+
           scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
           ggtitle(oneMCl)+
           theme_bw()+
           theme(legend.position="none",
                 plot.title=element_text(hjust=0.5, size=10), 
                 axis.title=element_text(size=6),
                 axis.text=element_text(size=6))
   fig0 <- ggMarginal(fig0, groupColour=T, groupFill=F, size=2)
   fig0
})

lab2 <- c("CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS+DEX",
           "PHA"="PHA", "PHA-DEX"="PHA+DEX")
legend2 <- get_legend(
   ggplot(df2%>%filter(MCls=="Bcell"), aes(LDA_1,LDA_2))+
   geom_point(aes(colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_rect(colour=NA, fill=NA),
         legend.text=element_text(size=8),
         legend.key=element_rect(fill=NA),
         legend.key.size=grid::unit(1,"lines"))
   )  
      
##
png("./9_RNA.dynamic2_output/Old/Figure5.X.png", width=900, height=700, res=130)  

fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]], 
                  fig_ls[[3]], fig_ls[[4]], 
                  nrow=2, ncol=2, align="hv",axis="tb")#+
#        draw_plot_label(c("LDA_1","LDA_2"), 
#                        x=c(0.4,0), y=c(0.025,0.5), angle=c(0,90),
#                        fontface="plain", size=8)
#fig2 <- plot_grid(legend2, NULL, nrow=2, rel_heights=c(1,3))        
                       
print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


cat("(2).", "By cell type separately", "\n")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
           
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- map_dfr(MCls, function(oneMCl){
###
   #cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
   meta <- read_rds(fn)      
   dd <- data.frame(LDA_1=meta$LDA_1, LDA_2=meta$LDA_2, treats=meta$treats, MCls=oneMCl)
   dd
})
df2$treat2 <- gsub("-EtOH", "", df2$treats)

### loop for plotting cell type respectively
for (i in 1:length(MCls)){
   oneMCl <- MCls[i]
   cat(i, oneMCl, "\n")
   df0 <- df2%>%filter(MCls==oneMCl)
   #if ( oneMCl=="NKcell") df0$LDA_2 <- -df0$LDA_2
   fig0 <- ggplot(df0, aes(x=LDA_1, y=LDA_2))+
           geom_point(aes(colour=factor(treat2)), size=0.1)+
           scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
           theme_bw()+
           theme(legend.position="none",
                 axis.title=element_text(size=8),
                 axis.text=element_text(size=6))
   fig0 <- ggMarginal(fig0, groupColour=T, groupFill=F, size=3)

   figfn <- paste("./9_RNA.dynamic2_output/Old/Figure5.", i, ".", oneMCl, ".png", sep="")
   png(figfn, width=500, height=500, res=150)         
   print(fig0)
   dev.off()
   Sys.sleep(60)

} ###

} ###



##################################
### 4.4, boxplot showing LDA 1 ###
##################################
if(FALSE){
cat("4.4.", "\n")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- map_dfr(MCls, function(oneMCl){
   #cat(oneMCl,"\n")
   fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
   meta <- read_rds(fn)      
   dd <- data.frame(LDA_1=meta$LDA_1, LDA_2=meta$LDA_2, 
                    treats=meta$treats, MCls=oneMCl, sampleID=meta$BEST.GUESS)
   dd
})
df2$treat2 <- gsub("-EtOH","",df2$treats)

df2 <- df2%>%group_by(sampleID, MCls, treat2)%>%summarise(y1=mean(LDA_1),y2=mean(LDA_2))

cols1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
            
###figure 1
fig1 <- ggplot(df2)+
        geom_boxplot(aes(x=treat2, y=y1, color=treat2))+
        scale_color_manual(values=cols1)+ylab("Average LDA_1")+
        facet_wrap(~MCls, nrow=2, scales="free")+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(), 
              strip.background=element_blank(),
              strip.text.x=element_text(size=12))
###
figfn <- "./9_RNA.dynamic2_output/Old/Figure5.1_LDA1.boxplot.png"
png(figfn, width=700, height=600, res=120)
print(fig1)
dev.off()

### figure 2
fig2 <- ggplot(df2)+
        geom_boxplot(aes(x=treat2, y=y2, color=treat2))+
        scale_color_manual(values=cols1)+ylab("Average LDA_2")+
        facet_wrap(~MCls, nrow=2, scales="free")+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(), 
              strip.background=element_blank(),
              strip.text.x=element_text(size=12))
###
figfn <- "./9_RNA.dynamic2_output/Old/Figure5.2_LDA2.boxplot.png"
png(figfn, width=700, height=600, res=120)
print(fig2)
dev.off() 

}






#####################################
### 5, calculate LDA used 100 PCs ###
#####################################
#if(FALSE){
#
#rm(list=ls())
#cat("5.1.", "calculate LDA used 100 PCs", "\n")
#### (1), calculate LDA 
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#tmp <- mclapply(MCls, function(oneMCl){
#   s1 <- Sys.time()
#   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
#   sc <- read_rds(fn)
#
#   #X <- sc@assays$RNA@data
#   X <- t(Embeddings(sc, reduction="pca"))
#   meta <- sc@meta.data
#   treats <- unique(meta$treats)
#
#   ### Sw, within classes scatter matrix 
#   xx_ls <- lapply(treats, function(one){
#      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
#      Xi <- X[,ii]
#      xc <- scale(t(Xi),scale=FALSE)
#      xx <- crossprod(xc)
#      return(xx)
#   })
#   Sw <- Reduce("+", xx_ls)
#
#   ### Sb, between classes scatter matrix 
#   mu <- apply(X, 1, mean)
#   xx_ls <- lapply(treats, function(one){
#      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
#      Xi <- X[,ii]
#      ni <- ncol(Xi)
#      mu_k <- apply(Xi, 1, mean)
#      xk <- mu_k-mu
#      xx <- tcrossprod(xk)*ni
#      return(xx)   
#   })
#   Sb <- Reduce("+",xx_ls)
#
#
#   ### pseudo inverse of matrix, Sw
#   sw.eigen <- eigen(Sw)
#   subi <- sw.eigen$values>1e-02
#   uw <-sw.eigen$vectors[,subi]
#   dw <- sw.eigen$values[subi]
#   sw2 <- uw %*% diag(1/dw) %*% t(uw)
#
#   ### generazied eigen decomposition
#   sb.eigen <- eigen(Sb,symmetric=T)
#   ub <- sb.eigen$vectors[,1:4]
#   db <-sb.eigen$values[1:4]
#   sb2 <- ub %*% diag(sqrt(db)) %*% t(ub)
#   sb2i <- ub %*% diag(1/sqrt(db)) %*% t(ub)
#   
#   ss <- sb2 %*% sw2 %*% sb2
#   ss.eigen <- eigen(ss,symmetric=T) 
#   ###  
#   w <- sb2i %*% ss.eigen$vectors[,1:4]
#   values <- ss.eigen$values
#
#   zk <- crossprod(X,w)
#   meta$LDA_1 <- zk[,1]
#   meta$LDA_2 <- zk[,2]
#   meta$LDA_3 <- zk[,3]
#   meta$LDA_4 <- zk[,4]
#   opfn <- paste("./9_RNA.dynamic2_output/Old/5.1_LDApca.", oneMCl, ".meta.rds", sep="")
#   write_rds(meta, opfn)
#   s2 <- Sys.time()
#   d12 <- difftime(s2,s1,units="mins")
#   cat(oneMCl, ":", d12, "\n")
#   return(values)
#}, mc.cores=5)
#
#save(tmp, file="./9_RNA.dynamic2_output/Old/5.2_tmpPCA.eigenvalue.RData")
#
#} ##5.1 , End
#
#############################
#### (2) show decay curve ###
#############################
#if (FALSE){
#rm(list=ls())
#cat("5.2.", "Show decay curve", "\n")
#load("./9_RNA.dynamic2_output/Old/5.2_tmpPCA.eigenvalue.RData")
#tmp <- do.call(cbind,tmp) 
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#
####show proportion of LDA
#x <- sweep(tmp,2,colSums(tmp),"/")
#colnames(x) <- MCls
####
#fig_ls <- lapply(MCls, function(ii){
#   xi <- x[,ii]
#   di <- data.frame(x=1:10, y=xi[1:10])
#   p <- ggplot(di,aes(x=x,y=y))+
#        geom_point()+xlab("Nums of LDA")+ylab("Proportion")+
#        ggtitle(ii)+
#        theme_bw()+
#        theme(plot.title=element_text(hjust=0.5))
#   p
#})
#
####
#figfn <- paste("./9_RNA.dynamic2_output/Old/Figure7.1.pca.decay.png", sep="")
#png(figfn,width=750, height=700, res=120)
#print(plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], ncol=2))
#dev.off() 
#
#
#
###############################
#### (3). Show LDA1 and LDA2 ###
###############################  
##meta$zk <- zk
####
#cat("5.3.", "Show LDA1 and LDA2", "\n")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#           
#           
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#df2 <- lapply(MCls, function(oneMCl){
####
#   #cat(oneMCl,"\n")
#   fn <- paste("./9_RNA.dynamic2_output/Old/5.1_LDApca.", oneMCl, ".meta.rds", sep="")
#   meta <- read_rds(fn)      
#   dd <- data.frame(LDA_1=meta$LDA_1, LDA_2=meta$LDA_2, treats=meta$treats, MCls=oneMCl)
#   dd
#})
#df2 <- do.call(rbind,df2)
#df2$treat2 <- gsub("-EtOH","",df2$treats)
#
#fig0 <- ggplot(df2, aes(x=LDA_1, y=LDA_2))+
#        geom_point(aes(colour=factor(treat2)), size=0.1)+
#        facet_wrap(~factor(MCls), ncol=2, scales="free")+
#        scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=3)))+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.background=element_rect(colour=NA, fill=NA),
#              legend.key=element_rect(fill=NA),
#              strip.background=element_blank(),
#              strip.text.x=element_text(size=12))
#              
#png("./9_RNA.dynamic2_output/Old/Figure7.2.pca.LDA12.png", width=800, height=650, res=120)
#print(fig0)
#dev.off() 
#
#} ##5, End
#


#####################################################
### 6. integrating data for dynamical QTL mapping ###
#####################################################

##
## average LDA
#tmp <- lapply(MCls, function(oneMCl){
###
#fn <- paste("./9_RNA.dynamical3_output/9_LDA.", oneMCl, ".meta", sep="")
#meta <- read_rds(fn)
#dd <- meta%>%select(BEST.GUESS, treats, LDA_1)
#dd <- dd %>%group_by(BEST.GUESS,treats)%>%summarize(LDA_ave=mean(LDA_1))

#opfn <- paste("./9_RNA.dynamical3_output/9_AVE.",oneMCl, ".txt",sep="")
#write.table(dd, file=opfn, row.names=F, quote=F, sep="\t")
#})

##################################################
### 6.1. defined Three bins according to LDA_1 ###
##################################################
if(FALSE){
rm(list=ls())

cat("6.1.", "Define 3 bins\n")
###Fun, class Bin   
Binfun <- function(LDA){
   x <- quantile(LDA,probs=c(1/3,2/3))
   Bin <- rep(1,length(LDA))
   Bin[(LDA>=x[1]&LDA<x[2])] <- 2 
   Bin[LDA>=x[2]] <- 3
   Bin
}

####
#cutBin <- function(meta, breaks=3){
#   meta <- meta%>%arrange(LDA)
#   L1<- cut(x, breaks, labels=1:breaks)
#   d1 <- data.frame(x, L1)     
#   meta <- meta%>%mutate(Bin=d1[,2])
#}

### class bin
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){
###
fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
meta <- read_rds(fn)
#if( oneMCl!="NKcell") meta$LDA_2 <- -meta$LDA_2
treats <- unique(meta$treats)
tmp <- lapply(treats, function(ii){
   meta0 <- meta%>%filter(treats==ii)
   meta0$Bin1 <- Binfun(meta0$LDA_1)
   meta0$Bin2 <- Binfun(meta0$LDA_2)
   meta0
})
meta <- do.call(rbind,tmp)
###
opfn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
write_rds(meta,opfn)
}

} ##6.1, End



############################################
### 6.2, average GE for each combination ###
############################################
if(TRUE){
rm(list=ls())

### (1). get counts data
cat("(1).", "get counts data", "\n")
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
count <- sc@assays$RNA@counts         

grch38_unq <- grch38%>%
              distinct(ensgene, .keep_all=T)%>%
              dplyr::select(ensgene, symbol, chr, biotype)
vars <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn),
               ensgene2=gsub("[SU]-", "", rn), 
               uns=grepl("S-",rn), rnz=rowSums(count))%>%
        left_join(grch38_unq, by="ensgene")

autosome <- as.character(1:22)        
varsSel <- vars%>%filter(uns, rnz>0, chr%in%autosome)          

count <- count[varsSel$rn,]
rownames(count) <- varsSel$ensgene2


### 
prefix <- "./9_RNA.dynamic2_output/Old/LDA2Bin/"
dir.create(prefix, showWarnings=F)

### (2), average expression value
cat("(2).", "average expression value by bin", "\n")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){
fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
meta0 <- read_rds(fn)
count0 <- count[,meta0$NEW_BARCODE]

bti <- meta0%>%transmute(bti=paste(BEST.GUESS, MCls, treats, BATCH, Bin2, sep="_")) %>% unlist %>% factor
ncell <- data.frame(x=bti)%>%
         group_by(x)%>%
         summarise(ncell=n())
X <- model.matrix(~0+bti)
colnames(X) <- gsub("bti","",colnames(X))
rownames(X) <- meta0$NEW_BARCODE

x1 <- as.character(ncell$x)
x2 <- as.character(colnames(X))
x3 <- as.character(colnames(count0))
x4 <- as.character(rownames(X))
cat(oneMCl, "bti", identical(x1,x2), "\n")
cat(oneMCl, "cells", identical(x3,x4), "\n")
#
####
YtX <- count0 %*% X 
YtX <- as.matrix(YtX)
opfn1 <- paste(prefix, "YtX.", oneMCl, ".sum.RData", sep="")
save(YtX, file=opfn1)
###
YtX_ave <- sweep(YtX, 2, ncell$ncell, "/")
opfn2 <- paste(prefix, "YtX.", oneMCl, ".ave.RData", sep="")
save(YtX_ave, file=opfn2)

opfn3 <- paste(prefix, "0_ncell.", oneMCl, ".ave.RData", sep="")
save(ncell, file=opfn3)
} ## (2)

} ###6.2
###
###


#########
### 7 ###
#########
##
if(FALSE){
getData <- function(MCls, index=1, top=50){

   df2 <- map_dfr(MCls, function(oneMCl){
      fn <- paste("./9_RNA.dynamic2_output/Old/4.2_LDA.", oneMCl, ".loading.rds", sep="")
      w <- read_rds(fn)
      wi <- sort(abs(w[,index]), decreasing=T)
      gene0 <- gsub("\\.[0-9]*", "", names(wi))
      gene0 <- gene0[1:top]
      df0 <- bitr(gene0, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
      df0 <- df0%>%
             mutate(MCls=oneMCl)%>%
             filter(!grepl("LOC",SYMBOL))%>%distinct(ENSEMBL,.keep_all=T)
      df0$MCls <- oneMCl
      df0
   })   
   df2
}

###get data
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- getData(MCls, index=2, top=100)
opfn <- "./9_RNA.dynamic2_output/Old/7.2.1_top100.rds"
write_rds(df2, opfn)

### enrich results
load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
geneBG <- gsub("\\.[0-9].*", "", rownames(YtX))
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


fn <- "./9_RNA.dynamic2_output/Old/7.2.1_top100.rds"
df2 <- read_rds(fn)
cg <- compareCluster(ENTREZID~MCls, 
                     data=df2, 
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     universe=BgDf$ENTREZID, 
                     ont="ALL",
                     pvalueCutoff=0.1,
                     qvalueCutoff=0.5,
                     minGSSize=0,
                     maxGSSize=nrow(BgDf))
opfn <- "./9_RNA.dynamic2_output/Old/7.2.2_enrichGO.rds"
write_rds(cg,opfn)

### plots
cg <- read_rds("./9_RNA.dynamic2_output/Old/7.2.2_enrichGO.rds")                 
p1 <- dotplot(cg, x=~MCls, showCategory=5)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))

figfn <- "./9_RNA.dynamic2_output/Old/Figure8.2.1.png"
png(figfn,width=1000, height=1000, res=150)
print(p1)
dev.off() 

cg2 <- cg%>%
       filter(grepl("glucocorticoid|corticosteroid|lipopolysaccharide", Description))
p2 <- dotplot(cg2, x=~MCls, showCategory=5)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))

figfn <- "./9_RNA.dynamic2_output/Old/Figure8.2.2.png"
png(figfn,width=1300, height=1000, res=150)
print(p2)
dev.off() 

### number of cells
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell") 
#prefix <- "./9_RNA.dynamic2_output/Old/LDA2Bin/"
#for (oneMCl in MCls){
#fn <- paste("./9_RNA.dynamic2_output/Old/4_LDA.", oneMCl, ".meta.rds", sep="")
#meta0 <- read_rds(fn)
#bti <- meta0 %>% transmute(bti=paste(BEST.GUESS, MCls, treats, BATCH, Bin2, sep="_")) %>% unlist %>% factor
#ncell <- data.frame(x=bti)%>%
#         group_by(x)%>%
#         summarise(ncell=n())
#opfn <- paste(prefix, "ncell.", oneMCl, ".RData", sep="")
#save(ncell,file=opfn)
#}
} ##

###
###


########################
### 7, new scale data ##
########################
#if (FALSE){
#outdir <- "./9_RNA.dynamic2_output/tmp1/"
#if (!file.exists(outdir)) dir.create(outdir, showWarnings=F) 
#
#countFun <- function(X, layer="S"){
#   layer <- paste(layer, "-", sep="")
#   rn <- rownames(X)
#   rn1 <- rn[grepl(layer, rn)]
#   x <- X[rn1,]
#   counts <- colSums(x)
#   counts_after <- median(counts)
#   counts_after <- ifelse(counts_after==0, counts_after+1, counts_after)
#   counts <- ifelse(counts!=0, counts/counts_after, counts+1)
#   counts 
#}
#
#nomalize_log1p <- function(X, counts){
#   counts <- counts[colnames(X)]
#   X <- sweep(X, MARGIN=2, STATS=counts, FUN="/")
#   X <- log2(X+1)
#}
#
#fn <- "./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds"
#sc <- read_rds(fn)
#X <- sc@assays$RNA@data
#Xs <- X[grepl("S-",rownames(X)),]
#counts <- countFun(Xs)
#
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#tmp <- mclapply(MCls, function(oneMCl){
#
#   s1 <- Sys.time()
#   fn <- paste("./9_RNA.dynamic2_output/Old/3_MCl.", oneMCl, ".old.rds", sep="")
#   sc <- read_rds(fn)
#
#   X <- sc@assays$RNA@data 
#   ### here we directly used counts data. Maybe we need to try to normalize data first. 
#   #X <- t(Embeddings(sc, reduction="pca"))
#   meta <- sc@meta.data
#   treats <- unique(meta$treats)
#
#   ### Sw, within classes scatter matrix 
#   xx_ls <- lapply(treats, function(one){
#      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
#      Xi <- X[,ii]
#      xc <- scale(t(Xi),scale=FALSE)
#      xx <- crossprod(xc)
#      return(xx)
#   })
#   Sw <- Reduce("+", xx_ls)
#   cat("Within", 1, "\n")
#   
#   ### Sb, between classes scatter matrix 
#   mu <- apply(X, 1, mean)
#   xx_ls <- lapply(treats, function(one){
#      ii <- as.character(meta[meta$treats==one,"NEW_BARCODE"])
#      Xi <- X[,ii]
#      ni <- ncol(Xi)
#      mu_k <- apply(Xi, 1, mean)
#      xk <- mu_k-mu
#      xx <- tcrossprod(xk)*ni
#      return(xx)   
#   })
#   Sb <- Reduce("+",xx_ls)
#   cat("Between", 2, "\n")
#
#   ### pseudo inverse of matrix, Sw
#   sw.eigen <- eigen(Sw)
#   subi <- sw.eigen$values>1e-02
#   uw <-sw.eigen$vectors[,subi]
#   dw <- sw.eigen$values[subi]
#   sw2 <- uw %*% diag(1/dw) %*% t(uw)
#
#   ### generazied eigen decomposition
#   sb.eigen <- eigen(Sb,symmetric=T)
#   ub <- sb.eigen$vectors[,1:4]
#   db <-sb.eigen$values[1:4]
#   sb2 <- ub %*% diag(sqrt(db)) %*% t(ub)
#   sb2i <- ub %*% diag(1/sqrt(db)) %*% t(ub)
#   
#   ss <- sb2 %*% sw2 %*% sb2
#   ss.eigen <- eigen(ss,symmetric=T) 
#   ###  
#   w <- sb2i %*% ss.eigen$vectors[,1:4]
#   values <- ss.eigen$values
#   
#   ### the top 4 LDA
#   #zk <- crossprod(X,w)
#   #meta$LDA_1 <- zk[,1]
#   #meta$LDA_2 <- zk[,2]
#   #meta$LDA_3 <- zk[,3]
#   #meta$LDA_4 <- zk[,4]
#   #opfn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="")
#   #write_rds(meta, opfn)
#   
#   ### loading information, weights of genes for the top 4 LDA
#   #w <- as.matrix(w)
#   #rownames(w) <- gsub("S-", "", rownames(X))
#   #colnames(w) <- paste("LDA", 1:4, sep="_")
#   #opfn <- paste("./9_RNA.dynamic2_output/Old/4.2_LDA.", oneMCl, ".loading.rds", sep="")
#   #write_rds(w, opfn)
#      
#   s2 <- Sys.time()
#   d12 <- difftime(s2,s1,units="mins")
#   cat(oneMCl, ":", d12, "\n")
#   return(values)
#},mc.cores=5)
#save(tmp, file="./9_RNA.dynamic2_output/Old/4.3_tmp.eigenvalue.RData")
#
#} ##4.1, End 






