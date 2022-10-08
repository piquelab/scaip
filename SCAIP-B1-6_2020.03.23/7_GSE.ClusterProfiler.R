###
library(tidyverse)
library(purrr)
library(furrr) 
library("BiocParallel")
library(qqman)
library(qvalue)
##
library(DESeq2)
library(annotables)
library(biobroom)
library(clusterProfiler)#,lib.loc="/wsu/home/groups/piquelab/apps/el7/R/4.0.0/lib64/R/library")
library(org.Hs.eg.db)
###
library(ggplot2)
## library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
 
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

theme_set(theme_grey())

rm(list=ls())

outdir <- "./7_GSE.ClusterProfiler_output/Filter2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)



#############################################################
### Enrichment analysis for differentially expressed gene ###
###     12-29-2020, last modified by Julong wei,           ###
#############################################################


##################################################################
### Example code for enrichment analysis using ClusterProfiler ###
##################################################################
#load("./6_DEG_CelltypeNew_output/YtX.comb.RData")
#geneBG <- rownames(YtX)
#geneBG.df <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

###example
#load("./6_DEG_CelltypeNew_output/Sigs.gene.DEG.RData")
#gene <- sigs[["LPS"]][["Monocyte"]]
#gene.df <- bitr(gene,fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),OrgDb=org.Hs.eg.db)
  
#ego2 <- enrichGO(gene=gene.df$ENTREZID,
#                 universe=geneBG.df$ENTREZID,
#                 #keyType="ENSEMBL",
#                 OrgDb=org.Hs.eg.db,
#                 ont="ALL",
#                 minGSSize=1,
#                 maxGSSize=nrow(gene.df),
#                 pvalueCutoff=1, qvalueCutoff=1, readable=T)
#ego3 <- gseGO(gene=gene.df$ENTREZID,
#              universe=geneBG.df$ENTREZID,
#              OrgDb=org.Hs.eg.db,
#              ont="ALL",
#              nPerm=1000,
#              minGSSize=1,
#              maxGSSize=nrow(gene.df),
#              pvalueCutoff=1, qvalueCutoff=1, verbose=FALSE)


#####################################################
### 0, parse gtt3 and obtain protein coding genes ###
#####################################################
#if(F){
#rm(list=ls())
#cat("0.", "Extract protein coding genes", "\n")
##/wsu/home/groups/piquelab/SCAIP/covariates/gencode.v34lift37.basic.annotation.gff3.gz
#fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
#con <- gzfile(fn,"rt")
#gff3 <- read.table(con,header=F)
#ID <- gsub(".*gene_id=", "", gff3$V9)
#ID <- gsub(";.*", "", ID)
#gene_type <- gsub(".*gene_type=", "", gff3$V9)
#gene_type <- gsub(";.*","",gene_type)
#
#anno <- gff3[,-9]%>%
#        mutate(ID=ID, gene_type=gene_type)%>%
#        filter( grepl("protein_coding",gene_type))%>%
#        distinct(ID, .keep_all=T)
#
#### 19,957 protein coding gene                 
#pList <- gsub("\\..*", "", anno$ID) 
#}

#######################################
### 1, gene set enrichment analysis ###
#######################################
#if (FALSE){
#### background gene list
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#geneBG <- gsub("\\..*", "", rownames(YtX))
#geneBG0 <- geneBG[geneBG%in%pList]  ##17,690
#geneBG.df <- bitr(geneBG0, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
#
###DEG gene list, Among 6413 DEG genes, 5930 protein gene
#load("./6_DEG_CelltypeNew_output/Sigs.gene.DEG.RData")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#geneList <- list()
#for(icont in contrast){
#   for (oneMCl in MCls){
#      gene <- sigs[[icont]][[oneMCl]]
#      gene0 <- gene[gene%in%pList]
#      gene.df <- bitr(gene0, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
#      comb <- paste(icont, oneMCl, sep="_")      
#      geneList[[comb]] <- gene.df$ENTREZID
#   }
#}
#   
#
#cg <- compareCluster(geneCluster=geneList, 
#                     universe=geneBG.df$ENTREZID,
#                     fun="enrichGO", 
#                     OrgDb="org.Hs.eg.db", 
#                     ont="ALL")
#                     
#fig0 <- dotplot(cg%>%filter(ONTOLOGY=="BP"),showCategory=5)+
#        theme(axis.text.x=element_text(angle=90, hjust=1,size=10),
#              axis.text.y=element_text(size=10))
#
#figfn <- "./7_GSE.ClusterProfiler_output/Figure1.1.BP.png"
#png(figfn, width=1200, height=1000, res=120)
#fig0
#dev.off()
#}


######################################################
### 1, Enrichment for up and down genes separately ###
######################################################

### background gene list
load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
rownames(YtX_sel) <- gsub("\\.[0-9].*", "", rownames(YtX_sel))

grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
anno <- data.frame(ensgene=rownames(YtX_sel))%>%mutate(rnz=rowSums(YtX_sel))%>%left_join(grch38_unq, by="ensgene")
geneBG <- anno%>%dplyr::pull(ensgene)##15,770
#geneBG0 <- geneBG[geneBG%in%pList]  ##17,690
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0,"Up", "Down"))


### GO enrichment 
# "GO analysis for direction (up and down) directly", "\n")
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

#cg0 <- as.data.frame(cg)
cg <- cg%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(cg,"./7_GSE.ClusterProfiler_output/Filter2/old/1_enrichGO.rds")

###
### combine direction of up and down
## "GO analysis for direction together" 
cg <- compareCluster(ENTREZID~contrast+MCls,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

write_rds(cg,"./7_GSE.ClusterProfiler_output/Filter2/1.2_enrichGO.combDirection.rds")

                                                 
###(3), KEGG enrichment analysis
cat("(3)", "KEGG analysis", "\n")
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=1,
                     qvalueCutoff=1,
                     minGSSize=0,
                     maxGSSize=1000)

#ck0 <- as.data.frame(ck)
ck <- ck%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(ck,"./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
              

######################################################
### 2. Show general figures of enrichment results  ###
######################################################

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
cluster2 <- setNames(cl2, cl)
lab2 <- setNames(gsub("-","+",cl),cl2)
        
                      
###(1), showing all the GO
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)]  
cg2 <- cg%>%mutate(ClusterNew=x,
                   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
   filter(maxGSSize>5&maxGSSize<500, p.adjust<0.1)                  
fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
              axis.text.y=element_text(size=10))

### pdf
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.GO.pdf"
pdf(figfn, width=15, height=10)
print(fig1)
dev.off()

### png
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.GO.png"
png(figfn, width=3500, height=2000, res=180)
print(fig1)
dev.off()


### (3). showing KEGG
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x,
                   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
  filter(maxGSSize>5&maxGSSize<500, p.adjust<0.1)

fig3 <- enrichplot::dotplot(ck2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
              axis.text.y=element_text(size=10))
        
### pdf
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure2.KEGG.pdf"
pdf(figfn, width=15, height=10)
print(fig3)
dev.off()

### png
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure3.KEGG.png"
png(figfn, width=3000, height=1500, res=150)
print(fig3)
dev.off()

### End


########################################
### response to reviewers, 6-3-2022, ###
########################################
rm(list=ls())
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cl2 <- 1:32

## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
cluster2 <- cl2
names(cluster2) <- cl

lab2 <- cl
names(lab2) <- cl2

###
###
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")    
    
###(1), showing all the GO
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg <- cg%>%    
    mutate(Cluster=gsub("-", "+", Cluster),
           contrast=gsub("-", "+", contrast),
           col.contrast=col1[contrast],
           col.MCls=col2[MCls],
           ClusterValue=as.numeric(cluster2[Cluster]),
           ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction}"),
## ClusterNew=paste(contrast, MCls, direction, sep="."),         
ClusterNew=fct_reorder(ClusterNew, ClusterValue),
           maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))

cg2 <- cg%>%filter(maxGSSize>5&maxGSSize<500, p.adjust<0.1)               
##
fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
        theme(axis.title=element_blank(),
              axis.text.x=element_markdown(angle=60, hjust=1, size=15),
              ## axis.text.x=element_text(angle=60, hjust=1, size=15),
              axis.text.y=element_text(size=12))


 
### pdf
## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.2_GO_reviews.pdf"
## pdf(figfn, width=18, height=12)
## print(fig1)
## dev.off()

###
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.2_GO_reviews.png"
png(figfn, width=3500, height=2200, res=180)
print(fig1)
dev.off()




### (3). showing KEGG
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
ck <- ck%>%    
    mutate(Cluster=gsub("-", "+", Cluster),
           contrast=gsub("-", "+", contrast),
           col.contrast=col1[contrast],
           col.MCls=col2[MCls],
           ClusterValue=as.numeric(cluster2[Cluster]),
           ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction}"),
## ClusterNew=paste(contrast, MCls, direction, sep="."),
ClusterNew=fct_reorder(ClusterNew, ClusterValue),
           maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))

ck2 <- ck%>%filter(maxGSSize>5&maxGSSize<500, p.adjust<0.1)               

fig3 <- enrichplot::dotplot(ck2, x="ClusterNew", showCategory=5)+
        theme(axis.title=element_blank(),
              axis.text.x=element_markdown(angle=60, hjust=1,size=14),
              ## axis.text.x=element_text(angle=60, hjust=1, size=14),
              axis.text.y=element_text(size=12))
        
### pdf
## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure2.2_KEGG_review.pdf"
## pdf(figfn, width=18, height=12)
## print(fig3)
## dev.off()


### png
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure2.2_KEGG_review.png"
png(figfn, width=2500, height=1500, res=150)
print(fig3)
dev.off()



### End




############################
### show combine results ###
############################

## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1.2_enrichGO.combDirection.rds")
## cg2 <- cg%>%mutate(ClusterNew=paste(MCls, contrast, sep="."),
##    maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
##    filter(maxGSSize<500, p.adjust<0.1)
## fig0 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=20)+
##         theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
##               axis.text.y=element_text(size=10))
## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.2.GO.png"
## png(figfn, width=2500, height=1500, res=120)
## print(fig0)
## dev.off()

##################################
### 3. show specific GO terms  ###
##################################


### background gene list
#load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
#rownames(YtX_sel) <- gsub("\\.[0-9].*", "", rownames(YtX_sel))
#geneBG <- rownames(YtX_sel)
#BG <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
#
#cg <- enrichGO(gene=BG$ENTREZID,
#                 universe=BG$ENTREZID,
#                 OrgDb=org.Hs.eg.db,
#                 ont="ALL",
#                 minGSSize=0,
#                 maxGSSize=nrow(BG),
#                 pvalueCutoff=1, qvalueCutoff=1)
#write_rds(cg, file="./7_GSE.ClusterProfiler_output/Filter2/4_GO.gene.rds")


#gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
#ens <- gene2%>%dplyr::pull(ENSEMBL)
## rm(list=ls())

## ExampleGOplot <- function(cg){
   
##    x <- str_split(cg$GeneRatio, "/", simplify=T)
##    GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
##    Drt2 <- c("Up"=1, "Down"=2) 
##    cg <- cg%>%mutate(Direction2=Drt2[as.character(direction)], 
##               contrast2=paste(contrast, Direction2, sep="."),
##               contrast2=gsub("-", "+", contrast2))
##    #             
##    cg$size <- rep(1,nrow(cg))
##    cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
##    cg$size[GeneRatio>=0.15] <- 3
##    #
##    ## cg$p2 <- cg$p.adjust
##    ## cg$p2[cg$p2>0.1] <- NA
   
##    fig0 <- ggplot(cg, aes(x=contrast2, y=MCls))+
##            geom_point(aes(size=factor(size), colour=p2))+
##            scale_x_discrete(labels=c("LPS.1"="LPS.Up", "LPS.2"="LPS.Down",
##                             "LPS+DEX.1"="LPS+DEX.Up", "LPS+DEX.2"="LPS+DEX.Down",
##                             "PHA.1"="PHA.Up", "PHA.2"="PHA.Down",
##                             "PHA+DEX.1"="PHA+DEX.Up", "PHA+DEX.2"="PHA+DEX.Down"))+
##            scale_colour_gradient(name="p.adjust",                           
##                                  low="blue", high="red", trans="reverse")+    #"#ffa500"
##            scale_size_manual(name="GeneRatio", 
##                              values=c("1"=1, "2"=2, "3"=3),
##                              labels=c("1"="<0.05", "2"="<0.15", "3"=">=0.15"))+
##            theme_bw()+
##            theme(axis.title=element_blank(),
##                  axis.text.x=element_text(angle=-90, size=9, hjust=0, vjust=0.5),
##                  axis.text.y=element_text(size=9),
##                  legend.background=element_blank(),
##                  legend.title=element_text(size=8),
##                  legend.text=element_text(size=6),
##                  legend.key.size=grid::unit(0.8,"lines"))
##    fig0
## }

## ### read data
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)

## ### type I interferon
## cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("Type I interferon signaling pathway")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.1_interferon.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 

## ### response to glucocorticoid ###
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(grepl("cellular response to glucocorticoid stimulus", Description))
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("cellular response to glucocorticoid stimulus")+
##         theme(plot.title=element_text(hjust=0.5, size=8))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.2_glucocorticoid.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ### inflammatory response
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="inflammatory response")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("inflammatory response")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.3_inflammatory.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ### cell activation involved in immune response 
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="cell activation involved in immune response")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("Cell activation involved in immune response")+
##         theme(plot.title=element_text(hjust=0.5, size=8))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.4_immuneresponse.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ### innate immune response
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="innate immune response")
## ##
## ## cg2$p2 <- cg2$p.adjust
## ## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("innate immune response")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.5_immuneresponse.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ### response to bacterium 

## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="response to bacterium")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("response to bacterium")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.6_bacterium.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ### cytokine-mediated signaling pathway
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="cytokine-mediated signaling pathway")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("cytokine-mediated signaling pathway")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.7_cytokine-mediated.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ### cytokine receptor binding
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="cytokine receptor binding")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("cytokine receptor binding")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.8_cytokine-receptor.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ###
## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
## cg0 <- as.data.frame(cg)
## cg2 <- cg0%>%filter(Description=="cytokine activity")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("cytokine actvity")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.9_cytokine-activity.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ###
## ### response to lipopolysaccharide
## cg2 <- cg0%>%filter(Description=="response to lipopolysaccharide")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("response to lipopolysaccharide")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.10_lipopoly.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 


## ###
## ### stress response to copper ion
## cg2 <- cg0%>%filter(Description=="stress response to copper ion")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("stress response to copper ion")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.11_copper.ion.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 

## ###
## ### stress response to metal ion
## cg2 <- cg0%>%filter(Description=="stress response to metal ion")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("stress response to metal ion")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.12_metal.ion.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ###
## ### 
## cg2 <- cg0%>%filter(Description=="myeloid leukocyte mediated immunity")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("myeloid leukocyte mediated immunity")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.13_myeloid.leukocyte.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ###
## ### 
## cg2 <- cg0%>%filter(Description=="granulocyte activation")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("granulocyte activation")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.14_granulocyte.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ###
## cg2 <- cg0%>%filter(Description=="myeloid cell activation involved in immune response")
## ##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(cg2)+
##         ggtitle("myeloid cell activation involved in immune response")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.15_myeloid.cell.activation.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


## ############
## ### KEGG ###
## ############

## ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds") 
## ck0 <- as.data.frame(ck)

## ### Influenza A
## ck2 <- ck0%>%filter(Description=="Influenza A")
## ck2$p2 <- ck2$p.adjust
## ck2$p2[ck2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(ck2)+
##         ggtitle("Influenza A")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.1_Influenza_A.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off() 

## ### cytokine-cytokine 
## ck2 <- ck0%>%filter(Description=="Cytokine-cytokine receptor interaction")
## ck2$p2 <- ck2$p.adjust
## ck2$p2[ck2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(ck2)+
##         ggtitle("Cytokine-cytokine receptor interaction")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.2_Cytokine.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()

## ### COVID-19
## ck2 <- ck0%>%filter(Description=="Coronavirus disease - COVID-19")
## ck2$p2 <- ck2$p.adjust
## ck2$p2[ck2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(ck2)+
##         ggtitle("Coronavirus disease-COVID-19")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.3_COVID-19.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()

## ### Asthma
## ck2 <- ck0%>%filter(Description=="Asthma")
## ck2$p2 <- ck2$p.adjust
## ck2$p2[ck2$p2>0.1] <- NA

## fig1 <- ExampleGOplot(ck2)+
##         ggtitle("Asthma")+
##         theme(plot.title=element_text(hjust=0.5, size=10))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.4_Asthma.png"
## png(figfn, width=500, height=400, res=120)
## print(fig1)
## dev.off()


###
###


###
### 
## rm(list=ls())
## load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
## geneBG <- gsub("\\..*", "", rownames(YtX_sel))
## df <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


## ck <- enrichKEGG(df$ENTREZID, organism='hsa',pvalueCutoff=1, qvalueCutoff=1)
## ck2 <- ck%>%filter(Description=="Asthma")
## geneID <- unlist(strsplit(ck2@result$geneID, split="/"))

## df2 <- df%>%filter(ENTREZID%in%geneID)


## cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds")
## cg2 <- cg%>%filter(grepl("type I interferon signaling pathway", Description))
## cg3 <- cg2%>%
##    filter(grepl("[(LPS)(PHA)].*Up", Cluster))
## geneList3 <- unique(unlist(str_split(cg3$geneID,"/")))
## cg4 <- cg2%>%
##    filter(grepl("DEX.*Down", Cluster))
## geneList4 <- unique(unlist(str_split(cg4$geneID,"/")))

## geneList <- unique(c(geneList3, geneList4))
## gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
## ens <- gene2%>%dplyr::pull(ENSEMBL)
