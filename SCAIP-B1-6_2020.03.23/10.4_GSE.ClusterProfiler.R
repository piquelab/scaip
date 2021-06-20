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
### 3. Show figures ###
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


###
###        



########################################
### compare shared GO by DVG and DMG ###
########################################
###calculate oods ration
odds.fun <-  function(df){
###   
   res <- map_dfr(1:nrow(df), function(i){
      Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
      Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
      dat <- data.frame(Diff=Diff, Bg=Bg)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=log2(fish$estimate),
                         CI.low=log2(fish$conf.int[1]),
                         CI.high=log2(fish$conf.int[2]))
      res0
   })
###
  df <- cbind(df, res)
  df  
}


cg_phi <-  read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew_enrichGO.rds")

cg_mu <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/4_mu_enrichGO.rds")

###
cg2_phi <- cg_phi%>%as.data.frame()%>%
    filter(Description=="cytokine-mediated signaling pathway")%>%
   dplyr::select(-ONTOLOGY, -ID, -Description, -geneID, -Count, -Cluster1)%>%
   separate(GeneRatio, c("Diff.in", "Diff.t"), convert=T)%>%
   mutate(Diff.not=Diff.t-Diff.in)%>%
   separate(BgRatio, c("Bg.in", "Bg.t"), convert=T)%>%
   mutate(Bg.not=Bg.t-Bg.in)

cg2_phi <- odds.fun(cg2_phi)%>%mutate(gr=1) 

###
cg2_mu <- cg_mu%>%as.data.frame()%>%
   filter(Description=="cytokine-mediated signaling pathway")%>%
   dplyr::select(-ONTOLOGY, -ID, -Description, -geneID, -Count, -Cluster1)%>%
   separate(GeneRatio, c("Diff.in", "Diff.t"), convert=T)%>%
   mutate(Diff.not=Diff.t-Diff.in)%>%
   separate(BgRatio, c("Bg.in", "Bg.t"), convert=T)%>%
   mutate(Bg.not=Bg.t-Bg.in)

cg2_mu <- odds.fun(cg2_mu)%>%mutate(gr=2)

df <- rbind(cg2_phi, cg2_mu)
rownames(df) <- NULL

col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

fig <- ggplot(df, aes(x=odds, y=Cluster))+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, color=MCls, alpha=factor(gr)),
                  size=0.5, height=0.2)+
   geom_point(size=2.5, color="grey40")+
   scale_colour_manual(values=col2)+
   scale_alpha_manual(values=c("1"=0.3, "2"=0.7),
      labels=c("1"="RNA varibility", "2"="RNA abundance"))+
   guides(colour="none")+
   ggtitle("cytokine-mediated signaling pathway")+ 
   scale_y_discrete("")+#, labels=Cluster)+
   scale_x_continuous("Odds ratio (log2 scale)")+
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5)) 

figfn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/Figure3.go.png"
png(figfn, width=700, height=600, res=120)
print(fig)
dev.off()
    

