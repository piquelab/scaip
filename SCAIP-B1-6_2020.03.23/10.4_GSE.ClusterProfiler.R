##
library(tidyverse)
library(purrr)
library(qvalue)
##
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
theme_set(theme_grey())


rm(list=ls())

outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


######################################################
### 1, Enrichment for up and down genes separately ###
######################################################

prefix <- "4_mu"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
 
### background gene list
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### geneCluster for enrichment results
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0, "Up", "Down"))
 
### GO enrichment
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- paste(outdir, prefix, "_enrichGO.rds", sep="")
write_rds(cg, opfn1)
                                                 
### (4)
### KEGG enrichment analysis
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     minGSSize=0,
                     maxGSSize=1000) #pval=0.2;qval=0.5
ck <- ck%>%mutate(Cluster1=Cluster)
opfn2 <- paste(outdir, prefix, "_enrichKEGG.rds", sep="")
write_rds(ck, opfn2)


#######################
### 2. Show figures ###
#######################

###3.1, MCl+treatment
# Show figures dotplots ordered by celltype+treatment

## prefix <- "4_mu"
## outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
       
## CL <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
##             rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep="_")
## lab2 <- setNames(gsub("-","+",CL), CL)
                      
## ### showing GO
## fn <- paste(outdir, prefix, "_enrichGO.rds", sep="")
## cg <- read_rds(fn)  
## cg2 <- cg%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))%>%filter(p.adjust<0.1)                 
## fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
##         scale_x_discrete(labels=lab2)+
##         theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
##               axis.text.y=element_text(size=15))
        
## figfn <- paste(outdir, prefix, "_Figure1.1_GO.png", sep="")
## png(figfn, width=2500, height=2000, res=180)
## print(fig1)
## dev.off()
       
## ### showing KEGG
## fn <-  paste(outdir, prefix, ".enrichKEGG.rds", sep="")
## ck <- read_rds(fn)
## ck2 <- ck%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))%>%filter(p.adjust<0.1)     
## fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=5)+
##         scale_x_discrete(labels=lab2)+
##         theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
##               axis.text.y=element_text(size=12))
        
## figfn <- paste(outdir, prefix, ".Figure1.2_KEGG.png", sep="")
## png(figfn, width=1500, height=1500, res=150)
## print(fig2)
## dev.off()       



###3.2, cell type+treatment+direction
### Show figures dotplots ordered by MCls+treat+direction"

prefix <- "4_mu"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")
ii <- paste(rep(c("A","B","C","D"),each=4), rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
cluster2 <- setNames(cl2, cl) 
lab2 <- setNames(gsub("-","+",cl),cl2)
                      
### GO
fn <- paste(outdir, prefix, "_enrichGO.rds", sep="")
cg <- read_rds(fn) 
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)] 
cg2 <- cg%>%mutate(ClusterNew=x,
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
   filter(maxGSSize<500, p.adjust<0.1)                

### png
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=12),
              axis.text.y=element_text(size=10))
        
figfn <- paste(outdir, prefix, "_Figure2.1_GO.png", sep="")
png(figfn, width=3000, height=2000, res=180)
print(fig1)
dev.off()

### pdf
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=12),
              axis.text.y=element_text(size=10))
        
figfn <- paste(outdir, prefix, "_Figure2.1_GO.pdf", sep="")
pdf(figfn, width=15, height=10)
print(fig1)
dev.off()
#figfn <- paste(outdir, prefix, "_Figure2.1_GO.pdf", sep="")
#pdf(figfn, width=15, height=10)


### KEGG
fn <-  paste(outdir, prefix, "_enrichKEGG.rds", sep="") 
ck <- read_rds(fn)
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x,
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%filter(maxGSSize<500, p.adjust<0.1)     

### png
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=10)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=10)) ##12 phiNew
        
figfn <- paste(outdir, prefix, "_Figure2.2_KEGG.png", sep="")
png(figfn, width=2000, height=1500, res=150) #180 phiNew
print(fig2)
dev.off()

## pdf
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=10)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=10)) ##12 phiNew
        
figfn <- paste(outdir, prefix, "_Figure2.2_KEGG.pdf", sep="")
pdf(figfn, width=13, height=10) #180 phiNew
print(fig2)
dev.off()


#################################
### 3. specific pathway terms ###
#################################

###
odds.fun <-  function(df){
###    
   res <- map_dfr(1:nrow(df), function(i){
      Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
      Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
      dat <- data.frame(Diff=Diff, Bg=Bg)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=as.numeric(fish$estimate),
                         CI.low=fish$conf.int[1],
                         CI.high=fish$conf.int[2])
      res0
   })
###
  df$odds <- res$odds
  df$CI.low <- res$CI.low
  df$CI.high <- res$CI.high
  df  
}

###
ExampleGOplot <- function(cg){

### prepare data    
   Drt2 <- c("Up"=1, "Down"=2) 
   cg <- cg%>%mutate(Direction2=Drt2[direction.x], 
      contrast2=paste(Direction2, contrast.x, sep="."))%>%
      mutate(contrast2=gsub("-", "+", contrast2))
   ##
   cg <- cg%>%drop_na(odds) 
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls.x))+
      geom_point(aes(size=odds, colour=p2))+
      scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
         "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
         "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
         "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      scale_colour_gradient(name="p.adjust",                           
         low="blue", high="red", trans="reverse", na.value=NA,
         guide=guide_colourbar(order=1), n.breaks=4)+    #"#ffa500"
      scale_size_binned("odds ratio",
         guide=guide_bins(show.limits=TRUE, axis=TRUE,
            axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
            n.breaks=4)+
      theme_bw()+
      theme(axis.title=element_blank(),
            axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
            axis.text.y=element_text(size=10),
            legend.background=element_blank(),
            legend.title=element_text(size=8),
            legend.text=element_text(size=6),
            legend.key.size=grid::unit(0.5, "lines"))
   fig0
}

contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
rn <- paste(rep(contrast, each=8), rep(rep(MCls, each=2), times=4),
                        rep(rep(c("Down", "Up"),times=4), times=4), sep=".")
tmp <- data.frame(contrast=rep(contrast, each=8),
                     MCls=rep(rep(MCls, each=2), times=4),
                     direction=rep(rep(c("Down", "Up"),times=4), times=4))%>%
       mutate(rn=paste(contrast, MCls, direction, sep="."))

### enrichment for DVG 
cg <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew_enrichGO.rds")
cg <- cg%>%as.data.frame()%>%
   mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
          Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
          Diff.not=Diff.total-Diff.in,
          Bg.in=as.numeric(gsub("/.*","", BgRatio)),
          Bg.total=as.numeric(gsub(".*/","", BgRatio)),
          Bg.not=Bg.total-Bg.in)
##
cg_mu <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/4_mu_enrichGO.rds")
cg_mu <- cg_mu%>%as.data.frame()%>%
   mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
          Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
          Diff.not=Diff.total-Diff.in,
          Bg.in=as.numeric(gsub("/.*","", BgRatio)),
          Bg.total=as.numeric(gsub(".*/","", BgRatio)),
          Bg.not=Bg.total-Bg.in)


###shared DEG and DVG
### DVG
pathway_ls <- c("type I interferon signaling pathway",
   "cytokine-mediated signaling pathway",
   "response to lipopolysaccharide",
   "innate immune response")

fig_ls <- lapply(pathway_ls, function(pathway){

### DEG
   mu2 <- cg_mu%>%filter(Description==pathway)
   mu2 <- odds.fun(mu2)
   mu2 <- mu2%>%full_join(tmp, by=c("Cluster"="rn"))
   mu2$p2 <- mu2$p.adjust
   mu2$p2[mu2$p2>0.1] <- NA

   fig1 <- ExampleGOplot(mu2)+
      ggtitle(paste(pathway, "(DEG)"))+
      theme(plot.title=element_text(hjust=0.5, size=10),
            legend.key.size=grid::unit(0.4, "lines"))
    
### DVG   
   cg2 <- cg%>%filter(Description==pathway)
   cg2 <- odds.fun(cg2)
   cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
   cg2$p2 <- cg2$p.adjust
   cg2$p2[cg2$p2>0.1] <- NA

   fig2 <- ExampleGOplot(cg2)+
      ggtitle(paste(pathway, "(DVG)"))+
      theme(plot.title=element_text(hjust=0.5, size=12),
            legend.key.size=grid::unit(0.6, "lines"))

   list(fig1, fig2)      
## figfn <- paste("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.",
##    ii, "_", pathway, "_DVG.png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


###
## figfn <- paste("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.",
##    ii, "_", pathway, "_DEG.png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig2)
## dev.off()
})

figfn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.0_shared.pdf"
pdf(figfn, width=10, height=10)
print(plot_grid(fig_ls[[1]][[1]], fig_ls[[1]][[2]],
   fig_ls[[2]][[1]], fig_ls[[2]][[2]],
   fig_ls[[3]][[1]], fig_ls[[3]][[2]],
   ## fig_ls[[4]][[1]], fig_ls[[4]][[2]],
   nrow=3, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()

############################
### DVG specific pathway ###
############################

## pathway <- "ribosome"
## pathway <- "protein targeting to ER"
## pathway <- "protein targeting to membrane"
## pathway <- "protein localization to endoplasmic reticulum"
## pathway <- "mRNA catabolic process"
pathway_ls <- c("ribosome",
   "protein targeting to membrane",
   "protein targeting to ER",
   "translation",
   "protein localization to endoplasmic reticulum",
   "mRNA catabolic process")
fig_ls <- lapply(pathway_ls, function(pathway){
## ### DEG
## mu2 <- cg_mu%>%filter(Description==pathway)
## mu2 <- odds.fun(mu2)
## mu2 <- mu2%>%full_join(tmp, by=c("Cluster"="rn"))
## mu2$p2 <- mu2$p.adjust
## mu2$p2[mu2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(mu2)+
##    ggtitle(paste(pathway, "(DEG)"))+
##    theme(plot.title=element_text(hjust=0.5, size=10))

###DVG
  cg2 <- cg%>%filter(Description==pathway)  
  cg2 <- odds.fun(cg2)
  cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
  cg2$p2 <- cg2$p.adjust
  cg2$p2[cg2$p2>0.1] <- NA

  fig2 <- ExampleGOplot(cg2)+
     ggtitle(pathway)+
     theme(plot.title=element_text(hjust=0.5, size=12))
  fig2
})

## figfn <- paste("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.",
##    ii, "_", pathway, "_DVG.png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 



### print
figfn <- paste("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure6.0_DVG.only.pdf", sep="")
pdf(figfn, width=10, height=10)
print(plot_grid(fig_ls[[1]], fig_ls[[2]],
   fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]],
   nrow=3, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()



###
###
pathway <- "ribosome"
cg2 <- cg%>%filter(Description==pathway)  
cg2 <- odds.fun(cg2)
cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.1] <- NA

p0 <- ExampleGOplot(cg2)+
   ggtitle("ribosome")+
   theme(plot.title=element_text(hjust=0.5, size=12))

figfn <-"./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure6.0.0_DVG.only.png"
png(figfn, width=500, height=400, res=120)
print(p0)
dev.off()






       
##############################################
### splite DEG and DVG enrichment analysis ###
##############################################

###
###
### my fun 1, generate data frame used for plots 
myDFxy <- function(dfx, dfy){
###
   dfx <- dfx%>%dplyr::select(beta, qval, rn, contrast, MCls, gene)
   dfy <- dfy%>%dplyr::select(beta, qval, rn)       
   dfxy <- dfx%>%inner_join(dfy, by="rn")
###        
   x <- dfxy$qval.x
   y <- dfxy$qval.y
   Bx <- abs(dfxy$beta.x)
   By <- abs(dfxy$beta.y)
   gr <- rep(1, nrow(dfxy))
   gr[(x<0.1&Bx>0.5)&((y>=0.1)|(y<0.1&By<=0.5))] <- 2
   gr[((x>=0.1)|(x<0.1&Bx<=0.5))&(y<0.1&By>0.5)] <- 3
   gr[(x<0.1&Bx>0.5)&(y<0.1&By>0.5)] <- 4
   dfxy$gr <- gr
###
   dfxy
}

outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
###
df1 <- read.table("./10_RNA.Variance_output/tmp9/3_phiNew.meta", header=TRUE)
df2 <- read.table("./10_RNA.Variance_output/tmp9/4_mu.meta", header=TRUE)
gr <- myDFxy(df1, df2)%>%dplyr::select(rn, gr)

### background gene list
res <- read.table("./10_RNA.Variance_output/tmp9/3_phiNew.meta", header=T)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### geneCluster for enrichment results 
res <- read.table("./10_RNA.Variance_output/tmp9/3_phiNew.meta", header=T)%>%
   left_join(gr, by="rn")%>%
   drop_na(beta,qval)%>%
   filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%
   inner_join(df0, by=c("gene"="ENSEMBL"))%>%
   mutate(direction=ifelse(beta>0, "Up", "Down"))
 
### GO enrichment
geneCluster0 <- geneCluster%>%filter(gr==2)
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster0,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- paste(outdir,  "5.1_unqDVG_enrichGO.rds", sep="")
write_rds(cg, opfn1)
                                                 
### (4)
### KEGG enrichment analysis
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster0,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     minGSSize=0,
                     maxGSSize=1000) #pval=0.2;qval=0.5
ck <- ck%>%mutate(Cluster1=Cluster)
opfn2 <- paste(outdir, "5.1_unqDVG_enrichKEGG.rds", sep="")
write_rds(ck, opfn2)


###############
### figures ###
###############

outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")
ii <- paste(rep(c("A","B","C","D"), each=4), rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
cluster2 <- setNames(cl2, cl) 
lab2 <- setNames(gsub("-","+",cl),cl2)
                      
### GO
fn <- paste(outdir, "5.1_unqDVG_enrichGO.rds", sep="")
cg <- read_rds(fn) 
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)] 
cg2 <- cg%>%mutate(ClusterNew=x,
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
   filter(maxGSSize<500, p.adjust<0.1)                

### png
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=12),
              axis.text.y=element_text(size=10))
        
figfn <- paste(outdir, "5.1_unqDVG_Figure1_enrichGO.png", sep="")
png(figfn, width=3000, height=2000, res=180)
print(fig1)
dev.off()


### KEGG
fn <-  paste(outdir, "5.1_unqDVG_enrichKEGG.rds", sep="") 
ck <- read_rds(fn)
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x,
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%filter(maxGSSize<500, p.adjust<0.1)     

### png
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=10)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=10)) ##12 phiNew
        
figfn <- paste(outdir, "5.1_unqDVG_Figure2_enrichKEGG.png", sep="")
png(figfn, width=2000, height=1500, res=150) #180 phiNew
print(fig2)
dev.off()


########################################
### compare shared GO by DVG and DMG ###
########################################
###calculate oods ration

## odds.fun <-  function(df){
## ###   
##    res <- map_dfr(1:nrow(df), function(i){
##       Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
##       Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
##       dat <- data.frame(Diff=Diff, Bg=Bg)
##       rownames(dat) <- c("in.category", "not.category")
##       fish <- fisher.test(dat)
##       res0 <- data.frame(odds=log2(fish$estimate),
##                          CI.low=log2(fish$conf.int[1]),
##                          CI.high=log2(fish$conf.int[2]))
##       res0
##    })
## ###
##   df <- cbind(df, res)
##   df  
## }


## cg_phi <-  read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew_enrichGO.rds")

## cg_mu <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/4_mu_enrichGO.rds")

## ###
## cg2_phi <- cg_phi%>%as.data.frame()%>%
##     filter(Description=="cytokine-mediated signaling pathway")%>%
##    dplyr::select(-ONTOLOGY, -ID, -Description, -geneID, -Count, -Cluster1)%>%
##    separate(GeneRatio, c("Diff.in", "Diff.t"), convert=T)%>%
##    mutate(Diff.not=Diff.t-Diff.in)%>%
##    separate(BgRatio, c("Bg.in", "Bg.t"), convert=T)%>%
##    mutate(Bg.not=Bg.t-Bg.in)

## cg2_phi <- odds.fun(cg2_phi)%>%mutate(gr=1) 

## ###
## cg2_mu <- cg_mu%>%as.data.frame()%>%
##    filter(Description=="cytokine-mediated signaling pathway")%>%
##    dplyr::select(-ONTOLOGY, -ID, -Description, -geneID, -Count, -Cluster1)%>%
##    separate(GeneRatio, c("Diff.in", "Diff.t"), convert=T)%>%
##    mutate(Diff.not=Diff.t-Diff.in)%>%
##    separate(BgRatio, c("Bg.in", "Bg.t"), convert=T)%>%
##    mutate(Bg.not=Bg.t-Bg.in)

## cg2_mu <- odds.fun(cg2_mu)%>%mutate(gr=2)

## df <- rbind(cg2_phi, cg2_mu)
## rownames(df) <- NULL

## col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##           "NKcell"="#aa4b56", "Tcell"="#ffaa00")

## fig <- ggplot(df, aes(x=odds, y=Cluster))+
##    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
##    geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, color=MCls, alpha=factor(gr)),
##                   size=0.5, height=0.2)+
##    geom_point(size=2.5, color="grey40")+
##    scale_colour_manual(values=col2)+
##    scale_alpha_manual(values=c("1"=0.3, "2"=0.7),
##       labels=c("1"="RNA varibility", "2"="RNA abundance"))+
##    guides(colour="none")+
##    ggtitle("cytokine-mediated signaling pathway")+ 
##    scale_y_discrete("")+#, labels=Cluster)+
##    scale_x_continuous("Odds ratio (log2 scale)")+
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5)) 

## figfn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.go.png"
## png(figfn, width=700, height=600, res=120)
## print(fig)
## dev.off()
    






