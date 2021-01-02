#
rm(list=ls())
source("./Bin/LibraryPackage.R")

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
#                 pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
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
#        distinct(ID, .keep_all=TRUE)
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
### 2, Enrichment for up and down genes separately ###
######################################################
if (FALSE){

cat("2.", "Enrichment analysis", "\n")
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


### (1) GO enrichment 
cat("2.1", "BP analysis", "\n")
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.1,
                     qvalueCutoff=0.2, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

#cg0 <- as.data.frame(cg)
cg <- cg%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(cg,"./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds")
                                                 
###(2), KEGG enrichment analysis
cat("2.2", "KEGG analysis", "\n")
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=0.1,
                     qvalueCutoff=0.2,
                     minGSSize=0,
                     maxGSSize=1000)

#ck0 <- as.data.frame(ck)
ck <- ck%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(ck,"./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
} ### 2, End  
              

#######################
### 3. Show figures ###
#######################

if (FALSE){
cat("3.", "Show figures dotplots", "\n")

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
cg2 <- cg%>%mutate(ClusterNew=x)                  
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
        
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.GO.png"
png(figfn, width=3500, height=2000, res=180)
print(fig1)
dev.off()

###(2), glucocorticoid
cg3 <- cg2%>%
       mutate(ii=grepl("glucocorticoid|corticosteroid|lipopolysaccharide", Description))%>%
       filter(ii,qvalue<0.01)
fig2 <- enrichplot::dotplot(cg3, x=~ClusterNew, showCategory=NULL)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
        
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.2.GO.png"
png(figfn, width=2000, height=1000, res=150)
print(fig2)
dev.off()       

###(3), type I interferon
cg3 <- cg2%>%
       filter(grepl("type I interferon", Description),qvalue<0.1)
       
fig3 <- enrichplot::dotplot(cg3, x=~ClusterNew, showCategory=NULL)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
        
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.3.GO.png"
png(figfn, width=2000, height=1000, res=150)
print(fig3)
dev.off() 


### (2), cluster the order
#cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
#cg0 <- as.data.frame(cg)
#
#x <- cg%>%group_by(Cluster)%>%dplyr::top_n(-10, wt=p.adjust)%>%as.data.frame()
#goID <- unique(x$ID)
#comb2 <- unique(cg0$Cluster)
#
#tmp <- map_dfc(comb2, function(ii){
#   x <- rep(0, length(goID))
#   names(x) <- goID
#   ID <- cg0%>%filter(Cluster==ii,p.adjust<0.05)%>%dplyr::pull(ID)
#   x[goID%in%ID] <- 1
#   x
#})
#tmp <- as.matrix(tmp)
#colnames(tmp) <- comb2
#rownames(tmp) <- goID
#hc <- hclust(dist(tmp))
##hcd <- as.dendrogram(hc)
##fig1 <- plot(hc) 
#dd <- cutree(hc, h=3.5)
#dd <- data.frame(bioID=names(dd),label=dd)
#opfn <- "./7_GSE.ClusterProfiler_output/Filter2/3_cutree.csv"
#write.csv(dd, opfn)
####dot plot

       
### (3). showing KEGG
ck <- read_rds("7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x)  
fig3 <- enrichplot::dotplot(ck2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
        
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure3.KEGG.png"
png(figfn, width=3000, height=1500, res=150)
print(fig3)
dev.off()
       
}



