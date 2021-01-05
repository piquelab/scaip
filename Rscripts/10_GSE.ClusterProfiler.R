rm(list=ls())

source("./Bin/LibraryPackage.R")
outdir <- "./10_RNA.Variance_output/GSE.ClusterProfiler/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

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


######################################################
### 2, Enrichment for up and down genes separately ###
######################################################
if (FALSE){

rm(list=ls())
cat("2.", "Enrichment analysis", "\n")
prefix <- "3_phiNew"

### (1) 
### background gene list
fn <- "./10_RNA.Variance_output/tmp7/3_phiNew.meta"
res <- read.table(fn, header=T)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


### (2)
### geneCluster for enrichment results
fn <- "./10_RNA.Variance_output/tmp7/3_phiNew.meta"
res <- read.table(fn, header=T)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0,"Up", "Down"))


### (3) 
### GO enrichment 
cat("2.1", "GO analysis", "\n")
cg <- compareCluster(ENTREZID~contrast+MCls, #+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.2,
                     qvalueCutoff=0.5, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- "./10_RNA.Variance_output/GSE.ClusterProfiler/3_phiNew.enrichGO.rds"
write_rds(cg, opfn1)
                                                 

### (4)
### KEGG enrichment analysis
cat("2.2", "KEGG analysis", "\n")
ck <- compareCluster(ENTREZID~contrast+MCls,#+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=0.2,
                     qvalueCutoff=0.5, 
                     minGSSize=0,
                     maxGSSize=1000)
ck <- ck%>%mutate(Cluster1=Cluster)
opfn2 <- "./10_RNA.Variance_output/GSE.ClusterProfiler/3_phiNew.enrichKEGG.rds"
write_rds(ck, opfn2)
} ### 2, End                


#######################
### 3. Show figures ###
#######################

if (TRUE){
Ni <- "3"
parm <- "phiNew"
cat("3.", "Show figures dotplots", "\n")

#cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
#            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
#cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")
#
#ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
#cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
#cluster2 <- setNames(cl2, cl) 
#lab2 <- setNames(gsub("-","+",cl),cl2)
         
CL <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep="_")
lab2 <- setNames(gsub("-","+",CL), CL)
                      
###(1)
fn <- "./10_RNA.Variance_output/GSE.ClusterProfiler/3_phiNew.enrichGO.rds"
cg <- read_rds(fn) 
#cg0 <- as.data.frame(cg)
#x <- cluster2[as.character(cg0$Cluster)] 
cg2 <- cg%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))  
#cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
#            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")               
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=15))
        
figfn <- "./10_RNA.Variance_output/GSE.ClusterProfiler/Figure1_phiNew.GO.png"
png(figfn, width=2500, height=2000, res=180)
print(fig1)
dev.off()


###(2)
#cg2 <- cg%>%
#       mutate(ii=grepl("glucocorticoid|corticosteroid|lipopolysaccharide", Description))%>%
#       filter(ii,qvalue<0.01)
#fig2 <- dotplot(cg0, x=~Cluster, showCategory=NULL)+
#        #facet_wrap(~direction,nrow=2)+
#        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
#              axis.text.y=element_text(size=10))
#        
#figfn <- "./10_RNA.Variance_output/GSE.ClusterProfiler/Figure1.1.BP.png"
#png(figfn, width=2000, height=1000, res=150)
#fig2
#dev.off()       
       
###(2)
fn <-  paste("./10_RNA.Variance_output/GSE.ClusterProfiler/3_phiNew.enrichKEGG.rds", sep="")
ck <- read_rds(fn)
#ck0 <- as.data.frame(ck)
#x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))     
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=12))
        
figfn <- "./10_RNA.Variance_output/GSE.ClusterProfiler/Figure2_phiNew.KEGG.png"
png(figfn, width=1500, height=1500, res=150)
print(fig2)
dev.off()
       
}



