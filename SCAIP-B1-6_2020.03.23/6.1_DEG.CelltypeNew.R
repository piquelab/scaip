###
rm(list=ls())
###
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
library(Rcpp)
library(reshape)
library(qqman)
library(qvalue)
##
library(Seurat)
library(DESeq2)
library(annotables)
library(biobroom)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(ggrastr)
##
theme_set(theme_grey())
outdir <- "./6_DEG.CelltypeNew_output/Filter2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#####################################################
v###    6. DEG analysis for each Cell type         ###
###    Batch one by one then meta analysis        ###
v###    12-25-2020, By Julong Wei, last modified   ###
#####################################################

#############################
### self-defined function ###
#############################

###my function (1)
myMeta <- function(dx){
   sub0 <-  !(is.na(dx[,"estimate"])|is.na(dx[,"stderror"])) 
   if (any(sub0)){
            baseMean <- dx[sub0,"baseMean"]
            b <- dx[sub0,"estimate"]
            std <- dx[sub0,"stderror"]
      ###
            vb_i <- 1/std^2
            w <- vb_i/sum(vb_i)
            baseMeanHat <- sum(w*baseMean)
            bhat <- sum(w*b)
            sdhat <- sqrt(1/sum(vb_i))
            z <- bhat/sdhat
            p <- 2*pnorm(abs(z),lower.tail=F)
            tmp <- c(baseMeanHat, bhat, sdhat, p)
   }else{
            tmp <- c(NA, NA, NA, NA)
   }
   tibble(baseMeanHat=tmp[1], beta=tmp[2],stderr=tmp[3],pval=tmp[4])
}

### my function (2)
myqval <- function(pval){
   qval <- pval
   ii0 <- !is.na(pval)
   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
   qval
}


##################################################
### 1, prepare data, pseudo-bulk seq data, YtX ###
##################################################

cat("1.", "Generate pseudo-bulk RNA seq data", "\n")
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
count <- sc@assays$RNA@counts         

##gene:64428, symbol:58149, ensgene:63697,
grch38_unq <- grch38%>%
              distinct(ensgene,.keep_all=T)%>%
              dplyr::select(ensgene, symbol, chr, start, end, biotype)
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn),
               ensgene2=gsub("[SU]-","",rn),
               uns=grepl("S-",rn), rnz=rowSums(count))%>%
        left_join(grch38_unq, by="ensgene")
        
#write.table(anno, file="./5_Identify_CelltypeOld_output/gene.infor", row.names=F, col.names=T, quote=F, sep="\t")        
#ii <- duplicated(grch38$ensgene)
#gene0 <- unique(grch38[ii,][["ensgene"]])
#nn <- sapply(gene0, function(one){
#x <- grch38[grch38$ensgene==one,]%>%select(ensgene, symbol,chr,start,end)
#comb <- unique(paste(x$symbol,x$chr,sep="_"))
#length(comb)
#})

autosome <- as.character(1:22)        
annoSel <- anno%>%filter(uns, rnz>0, chr%in%autosome)
### Out of 42,554 genes, 30,972 whole genes>20           
Ys <- count[annoSel$rn,]
rownames(Ys) <- annoSel$ensgene2 # rownames like ENSG00000187634.12

meta <- sc@meta.data
meta <- meta%>%
       mutate(bti=paste(MCls, treats, BEST.GUESS, BATCH, sep="_"))%>%
       dplyr::select(NEW_BARCODE, bti)

       
### summary the number of cells in 1536 combinations
dd <- meta %>%group_by(bti)%>%summarise(ncell=n())
save(dd, file="./6_DEG.CelltypeNew_output/0_ncell.RData")
#p <- ggplot(dd, aes(x=ncell))+
#     geom_histogram(alpha=0.5, color=NA, fill="blue")+xlab("No.cell")+
#     theme_bw()+
#     theme(plot.title=element_text(hjust=0.5))
#png("./6_DEG_CelltypeNew_output/Figure0.ncell.png", width=200, height=200, res=100)
#p
#dev.off()
### X, obs, bti
       
bti <- factor(meta$bti)       
X <- model.matrix(~0+bti)
YtX <- Ys %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))
opfn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
save(YtX, file=opfn)
###

### average 
#ncell <- colSums(X)
#YtXave <- sweep(YtX, 2, ncell, FUN="/")
#opfn <- "./6_DEG_CelltypeNew_output/YtX.ave.RData"
#save(YtXave, file=opfn)

### protein coding genes, rnz>20, combination with >20 cells
  
grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
anno <- data.frame(ensgene=gsub("\\..*", "", rownames(YtX)),ensgene2=rownames(YtX))%>%
        mutate(rnz=rowSums(YtX))%>%
        left_join(grch38_unq, by="ensgene")
geneSel <- anno%>%filter(rnz>20,grepl("protein_coding", biotype))%>%dplyr::pull(ensgene2)

load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)

YtX_sel <- YtX[geneSel, ddx$bti] ###rnz>20
opfn <- "./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
save(YtX_sel, file=opfn)


#######################################
### response to reviewers' cooments ###
#######################################
outdir <- "./6_DEG.CelltypeNew_output/Response_reviwers/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)
cvt <- str_split(ddx$bti, "_", simplify=T)

ddx <- ddx%>%mutate(MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3],
   treat2=gsub("-", "+", gsub("-EtOH", "", treats)))

MCls <- sort(unique(ddx$MCls))
fig_ls <- lapply(MCls, function(ii){
###
   dd2 <- ddx%>%filter(MCls==ii) 
   p1 <- ggplot(dd2, aes(x=treat2, y=ncell, fill=treat2))+
       geom_violin(aes(fill=treat2),width=1)+
       geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       ylab("Number of cells")+
       scale_fill_manual(values=c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
            "PHA"="#a6cee3", "PHA+DEX"="#1f78b4"))+
       ggtitle(ii)+
       theme_bw()+
       theme(legend.position="none",
             plot.title=element_text(hjust=0.5),
             axis.text=element_text(size=9),
             axis.title.x=element_blank())
    p1
})    
##

figfn <- paste(outdir, "Figure1_violin.png", sep="")
png(figfn, width=600, height=550, res=100)
plot_grid(plotlist=fig_ls, nrow=2, labels="AUTO", label_fontface="plain")
dev.off()






############################
## distinguish chemistry ###
############################
#meta <- sc@meta.data
#meta <- meta%>%
#       mutate(bti=paste(MCls, treats, BEST.GUESS, BATCH, chem, sep="_"))%>%
#       dplyr::select(NEW_BARCODE, bti)
       
#bti <- factor(meta$bti)       
#X <- model.matrix(~0+bti)
#YtX <- Ys %*% X
#YtX <- as.matrix(YtX)
#colnames(YtX) <- gsub("^bti", "", colnames(YtX))    

#ncell <- colSums(X)
#YtXave <- sweep(YtX, 2, ncell, FUN="/")
#opfn <- "./6_DEG_CelltypeNew_output/YtX.ave.2.RData"
#save(YtXave, file=opfn)   



##################################################
### 2, Differential gene expression analysis   ###
##################################################

cat("2.", "Differential gene expression analysis", "\n")
library("BiocParallel")  
library(DESeq2)
register(MulticoreParam(1))


### (1), prepare data
load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")  
rownames(YtX) <- gsub("\\.[0-9]*", "", rownames(YtX)) ##gene id, like ENSG00000187634 
cat("Genes,", nrow(YtX_sel), "Combination,", ncol(YtX_sel), "\n")
###
bti2 <- colnames(YtX_sel)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))

comb <- unique(cvt$comb)

### (2), one batch by batch in each cluster
MCls <- sort(unique(cvt$MCls))
Batch <- sort(unique(cvt$Batch2))
contrast.list <- list("LPS"=c("treats", "LPS-EtOH", "CTRL"),
                      "LPS-DEX"=c("treats", "LPS-DEX", "LPS-EtOH"),
                      "PHA"=c("treats", "PHA-EtOH", "CTRL"),
                      "PHA-DEX"=c("treats", "PHA-DEX", "PHA-EtOH"))
                     

##
### all the batch together
#cat("2.0", "run DESeq batch together", "\n")
### call DESeq
#for(oneMCl in MCls){ 
#   time0 <- Sys.time()
#   ##
#   cvt0 <- cvt %>% filter(MCls==oneMCl)
#   YtX0 <- YtX_sel[,cvt0$bti]
#   dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treats+Batch2)
#   dds <- DESeq(dds, parallel=T)
#   res <- contrast.list %>% map(~results(dds,contrast=.x, parallel=T))
#   res2 <- tibble(contrast=names(res),data=map(res,tidy)) %>% group_by(contrast) %>% unnest()  
   
#   opfn <- paste("./6_DEG_CelltypeNew_output/tmp2/Batch_", oneMCl, "_DESeq.results.txt", sep="")
#   write.table(res2, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
   ###
#   time1 <- Sys.time()
#   elapsed <- difftime(time1, time0, units="mins")
#   cat(oneMCl, elapsed, "Done\n")
#}

                                         
######################
### 2.1 call DESeq ###
######################
cat("2.1", "run DESeq batch separately", "\n")
   
res <- map_dfr(comb, function(oneX){   
   time0 <- Sys.time()
   cvt0 <- cvt%>% filter(comb==oneX) 
   oneMCl <- cvt0$MCls[1]
   oneBatch <- cvt0$Batch2[1]
   YtX0 <- YtX_sel[,cvt0$bti]
   dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treats)
   dds <- DESeq(dds, parallel=T)
   if (oneBatch %in% c("SCAIP2", "SCAIP3")){
      tmp <- results(dds, contrast=c("treats", "LPS-EtOH", "CTRL"), parallel=T)
      tmp2 <- tmp %>% tidy %>% mutate(contrast="LPS")
   }else{
      tmp <- contrast.list %>% map(~results(dds,contrast=.x, parallel=T))
      tmp2 <- tibble(contrast=names(tmp),data=map(tmp,tidy))%>%group_by(contrast)%>%unnest(c(data))  
   }
   tmp2 <- tmp2%>%mutate(Batch2=as.character(oneBatch),MCls=oneMCl)    
   
   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
   cat(oneX, elapsed, "Done\n")
   tmp2
})

opfn <- "./6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
write_rds(res,opfn)


##########################
### 2.2, meta analysis ###
##########################

#rm(list=ls())
cat("2.2", "meta analysis", "\n")

#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#res <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./6_DEG.CelltypeNew_output/tmp1/3_Batch1456_", oneMCl, "_meta.txt", sep="")
#   tmp <- read.table(file=fn, header=T)%>%mutate(MCls=oneMCl)
#   tmp
#})
#opfn <- "./6_DEG.CelltypeNew_output/tmp1/3_Batch1456.meta.rds"
#write_rds(res, opfn)
###

fn <- "./6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
res <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
       filter(Batch2%in%c("SCAIP1", "SCAIP4", "SCAIP5", "SCAIP6")) 
### (2) meta results
res2 <- res%>%group_by(rn)%>%
        nest()%>%
        mutate(outlist=mclapply(data, myMeta, mc.cores=1))%>%  #mutate(outlist=future_map(data,myMeta))%>%
        dplyr::select(-data)%>%unnest(c(outlist))%>%as.data.frame()           
cvt <- str_split(res2$rn, "_", simplify=T)
res2 <- res2%>%mutate(MCls=cvt[,1], contrast=cvt[,2], gene=cvt[,3])
           
### (3) add qvalue
res3 <- res2%>%group_by(MCls, contrast)%>%
        nest()%>%
        mutate(qval=map(data,~myqval((.x)$pval)))%>%
        unnest(c(data,qval))%>%as.data.frame()
#opfn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
opfn <- "./6_DEG.CelltypeNew_output/Filter2/3_Batch1456.meta.rds"
write_rds(res3, opfn)




############################################
### 3, Summary of results of DE analysis ###
############################################


##############################
### 3.1, figures, MA plots ###
##############################
## cat("3.1.", "Figure MA plots", "\n")
## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure1.1.MA.png"
## png(figfn, width=2000, height=2000, pointsize=12, res=300) 
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:16, 4, 4, byrow=T)
## layout(x)

## fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## res <- read_rds(fn)%>%mutate(color=ifelse(qval<0.1, T, F))%>%
##           dplyr::rename(baseMean=baseMeanHat)%>%drop_na(beta)

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## for (oneMCl in MCls){   
## ##1  
##    res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS")%>%dplyr::select(baseMean, beta, color, qval, pval)  
##    print(plotMA(res1[,1:3], colLine="NA", main="LPS", cex.main=1, cex.axis=0.8, cex.lab=1))

## ##2
##    res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")%>%dplyr::select(baseMean, beta, color, qval, pval)
##    print(plotMA(res2[,1:3], colLine="NA", main="LPS+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
## ##3
##    res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")%>%dplyr::select(baseMean, beta, color, qval, pval)
##    print(plotMA(res3[,1:3], colLine="NA", main="PHA", cex.main=1, cex.axis=0.8, cex.lab=1)) 
## ##4
##    res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")%>%dplyr::select(baseMean, beta, color, qval, pval)
##    print(plotMA(res4[,1:3], colLine="NA", main="PHA+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
   
##   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue") )        
## }
## dev.off()


## #####################
## ### 3.2, qq plots ###
## #####################
## cat("3.2.", "qq plots", "\n")
## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure1.2.qq.png"
## png(figfn, width=2000, height=2000, pointsize=12, res=300)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:16, 4, 4, byrow=T)
## layout(x)

## fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## res <- read_rds(fn)%>%
##        mutate(color=ifelse(qval<0.1, T, F))%>%
##        dplyr::rename(baseMean=baseMeanHat)%>%drop_na(beta)

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## for (oneMCl in MCls){   
## ##1  
##    res1 <- res%>%filter(MCls==oneMCl, contrast=="LPS")%>%dplyr::select(beta, color, qval, pval)  
##    print( qq(res1$pval, main="LPS", cex.main=1, cex.axis=0.8, cex.lab=1))

## ##2
##    res2 <- res%>%filter(MCls==oneMCl, contrast=="LPS-DEX")%>%dplyr::select(beta, color, qval, pval)
##    print( qq(res2$pval, main="LPS+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
## ##3
##    res3 <- res%>%filter(MCls==oneMCl, contrast=="PHA")%>%dplyr::select(beta, color, qval, pval)
##    print( qq(res3$pval, main="PHA", cex.main=1, cex.axis=0.8, cex.lab=1)) 
## ##4
##    res4 <- res%>%filter(MCls==oneMCl, contrast=="PHA-DEX")%>%dplyr::select(beta, color, qval, pval)
##    print( qq(res4$pval, main="PHA+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
##    print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))         
## }
## dev.off()

###
### 3.3, hist distribution of differential effects

rm(list=ls())

fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
dx <- read_rds(fn)%>%drop_na(beta)
dx <- dx%>%mutate(comb=paste(MCls, contrast, sep="_"),
                  gr=ifelse(qval<0.1, "sig", "not_sig"))
### for qq plots
comb <- unique(dx$comb)
dx <- map_dfr(comb, function(ii){
  di <- dx%>%filter(comb==ii)%>%arrange(pval)
  ngene <- nrow(di)
  di <- di%>%mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)))
  di
})  

###
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")

### qq plots
fig1 <- ggplot(dx, aes(x=expected,y=observed))+
   rasterise(geom_point(size=0.3, colour="grey30"),dpi=300)+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12)) 

### MA plots
## fig2 <- ggplot(dx, aes(x=baseMeanHat, y=beta))+
##    geom_point(aes(colour=factor(gr)), size=0.5)+
##    scale_colour_manual("", values=c("sig"="red", "not_sigs"="grey30"))+ 
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    geom_abline(colour="red")+
##    xlab("Mean expression")+
##    ylab(bquote(~log[2]~"fold change"))+
##    theme_bw()+
##    theme(legend.position="none", strip.text=element_text(size=12)) 

### hist
fig3 <- ggplot(dx, aes(x=beta))+
     geom_histogram(aes(y=..density..), fill="white", colour="grey30")+
     xlab(bquote(~log[2]~"fold change"))+
     ylab("Density")+
     facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
     theme_bw()+
     theme(strip.text=element_text(size=12))

## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure1.3.hist.png"
## png(filename=figfn, width=750, height=750, pointsize=12, res=120)  
## print(fig0)
## dev.off()
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure1.4_summryRes.pdf"
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig3, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()


#########################
### 3.3 boxplot show  ###
#########################
#if(FALSE){
#
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
#rn <- rownames(YtX) 
#rownames(YtX) <- gsub("\\..*", "", rn)
# 
#YtX0 <- YtX[sigs,]
#dx <- melt(YtX0)
#cvt <- str_split(dx$X2, "_", simplify=T)
#
#dx <- dx%>%
#      mutate(MCls=cvt[,1], treats=gsub("-EtOH","", cvt[,2]),sampleID=cvt[,3])%>%
#      dplyr::rename(ensgene=X1)
#
#dx2 <- dx%>%group_by(ensgene, MCls, treats)%>%summarise(y=mean(value))
#
#cols1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4") 
#cols2 <- c("Tcell"="#e41a1c", "NKcell"="#377eb8", 
#          "Bcell"="#4daf4a", "Monocyte"="#984ea3")
#          
#fig1 <- ggplot(dx2)+
#        geom_boxplot(aes(x=MCls, y=log2(y), color=MCls))+
#        scale_color_manual(values=cols2)+
#        ylab(bquote(log[2]~"(GE)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#fig2 <- ggplot(dx2)+
#        geom_boxplot(aes(x=treats, y=log2(y), color=treats))+
#        scale_color_manual(values=cols1)+
#        ylab(bquote(log[2]~"(GE)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#figfn <- "./6_DEG.CelltypeNew_output/Figure1.4.GE.boxplot.png"
#png(figfn, width=650, height=300, res=110)
#print(plot_grid(fig1, fig2, ncol=2, align="hv",
#                labels="AUTO", label_fontface="plain"))
#dev.off() 
#
#}
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
res <- res[,c(1,2,8,5,6,7,9)]
write.table(res, "./6_DEG.CelltypeNew_output/Differential_expression.txt", quote=F, sep="\t", row.names=F)

system("tar -czvf ./6_DEG.CelltypeNew_output/Differential_expression.txt.gz ./6_DEG.CelltypeNew_output/Differential_expression.txt")

########################################     
### 4, extract significant genes DEG ###
########################################

###4.1, Table of DEG
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1, abs(beta)>0.5)%>%drop_na(beta)
sigs <- unique(res$gene)
save(sigs, file="./6_DEG.CelltypeNew_output/Filter2/Sigs.gene.DEG.RData")

#######################################################################
### 4.2 show barplot DEG  from different contrast across cell types ###
#######################################################################

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

### (1). Show number of DEG with different contrast across cell type
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)

sigs <- res%>%
        group_by(contrast, MCls)%>%
        summarise(ngene=n())
x <- res%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))
x <- res%>%group_by(MCls)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

#cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
#          "NKcell"="#377eb8", "Tcell"="#e41a1c")
cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=MCls))+geom_bar(stat="identity", position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(name="", labels=lab2)+ylab("Nunber of DEG")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure2.1_DEG.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


## poster 
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure2.1.1_DEG.barplot.png"
png(filename=figfn, width=420, height=350, res=120)  
print(fig0+theme(legend.position="none"))
dev.off()


### barplots of DEG, up and down with light and deep colors, ***
###facet by MCls
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2,col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
        
## (3), barplot of DGE, facet by contrast and stacked up and down above x axis     
sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))

ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene))
#xpos <- c("Bcell"=0.5,"Monocyte"=1, "NKcell"=1.5, "Tcell"=2)
#ann2$xpos <- xpos[ann2$MCls]
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))

### 
fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col2comb,labels="")+ylab("DEG")+ylim(0,3500)+
        geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+100, fill=NULL), size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure2.3_DEG.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


###(4), barplot of DEG, facet by cell type, stacked up and down
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col1w <- colorspace::lighten(col1,0.3)
col1comb <- c(col1,col1w)
names(col1comb) <- paste(contrast, rep(c(1,2),each=4), sep="_") 

sig3 <- sigs%>%mutate(comb=paste(contrast, direction, sep="_"))
ann3 <- sig3%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene))
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## 
fig0 <- ggplot(sig3, aes(x=contrast, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col1comb, labels="")+ylab("DEG")+ylim(0,3500)+
        geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+100, fill=NULL), size=3)+
        scale_x_discrete(labels=lab2)+
        facet_grid(~MCls)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure2.4_DEG.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}

############################################################################ 
### (5) Final version, used for paper barplots of DEG, facet by contrast ###
###        and up above axis and down below axis                         ###
############################################################################   
Mybinom <- function(subdf){
   n1<- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
   n2 <- subdf%>%filter(direction==2)%>%dplyr::pull(ngene)
   ngene <- c(n1, n2)
   if(n1>n2){
      res <- binom.test(ngene, 0.5, alternative="greater")
   }else{
     res <- binom.test(ngene, 0.5, alternative="less") 
   }
   res$p.value
}
Mysymb <- function(pval){
  if(pval<0.001) symb <- "***"
  if(pval>=0.001 & pval<0.01) symb <- "**"
  if (pval>=0.01 & pval<0.05) symb <- "*"
  if(pval>0.05) symb <- ""
  symb
}
Mypos <- function(subdf){
 ny <- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
 ny
}
 
### facet by MCls
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2,col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_")

### read data
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
        
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2, -ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-1800,1800),5)
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
                          
###add star
anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
           mutate(pval=map_dbl(data, Mybinom), 
                  symb=map_chr(pval, Mysymb),
                  ypos=map_dbl(data, Mypos))%>%
           unnest(cols=c(contrast,MCls))                          

fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2))+
        geom_bar(aes(fill=comb),stat="identity")+
        scale_fill_manual(values=col2comb, labels="")+
        geom_hline(yintercept=0, color="grey60")+
        geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
                  vjust=ifelse(direction==2, 1.2, -0.2)), size=3)+        
        scale_y_continuous("", breaks=breaks_value, limits=c(-2000,2000),labels=abs(breaks_value))+         
        facet_grid(~contrast, labeller=facetlab)+   
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
fig0 <- fig0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure2.5_DEG.barplot5_final.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


##
#Mybinom <- function(subdf) {
#   n1<- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
#   n2 <- subdf%>%filter(direction==2)%>%dplyr::pull(ngene)
#   ngene <- c(n1, n2)
#   if(n1>n2){
#      res <- binom.test(ngene, 0.5, alternative="greater")
#   }else{
#     res <- binom.test(ngene, 0.5, alternative="less") 
#   }
#   res$p.value
#} 
#anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
#           mutate(pval=map_dbl(data,Mybinom))%>%
#           unnest(cols=c(contrast,MCls))



############################################
### 5,   heatmap for 16 conditions       ###
###  correlation between 16 conditions #####
##########################################

rm(list=ls())

load("./6_DEG.CelltypeNew_output/Filter2/Sigs.gene.DEG.RData")
DEGunq <- unique(sigs)

col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

######################
### (1) Using beta ###
######################  
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./6_DEG.CelltypeNew_output/Filter2/3_Batch1456.meta.rds"
res <- read_rds(fn)
TMP <- map_dfc(MCls, function(oneMCl){
   tmp <- map_dfc(Contrast, function(oneC){ 
      d0 <- res %>% filter(MCls==oneMCl, contrast==oneC)%>%mutate(z=beta/stderr)
      rownames(d0) <- d0$gene
      beta0 <- d0[DEGunq,"beta"]
      beta0
   })
   tmp 
})
rownames(TMP) <- DEGunq
TMP <- as.data.frame(TMP)
contrast2 <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
conditions <- paste(rep(MCls,each=4), rep(contrast2, times=4), sep="_")
names(TMP) <- conditions
ngene <- nrow(TMP)

###(1) heatmap
###mybreaks
ii <- rowSums(is.na(TMP))
TMP0 <- TMP[ii==0,]
y <- do.call(c, TMP0)
y0 <- y[abs(y)<2] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL

###colors
x <- str_split(conditions, "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- conditions
tmp_colors <- list(celltype=col2, treatment=col1) #brewer.pal(4,"Set1")

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)
fig1 <- pheatmap(TMP0, col=mycol, breaks=mybreaks, 
         scale="none",
         border_color="NA",
         cluster_rows=T, cluster_cols=T, 
         annotation_col=tmp_column,
         annotation_colors=tmp_colors, annotation_legend =T,
         show_colnames=T, show_rownames=F,
         na_col="white")

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure3.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
#              "Tcell_LPS+DEX", "Tcell_PHA+DEX", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
#              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#               "Tcell_LPS", "Tcell_PHA", "NKcell_LPS", "NKcell_PHA") 
Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
              "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
              "Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "NKcell_LPS+DEX", "NKcell_PHA+DEX", "Tcell_LPS+DEX", "Tcell_PHA+DEX")


corr <- cor(TMP0)[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)

x <- str_split(colnames(corr), "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)
tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F, na_col="white")

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure3.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()


###
library("grid")
library("ggplotify")
fig1 <- as.ggplot(fig1)
fig2 <- as.ggplot(fig2)
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure3.3_comb.heatmap.pdf"
pdf(figfn, width=10, height=5)
print(plot_grid(fig1, fig2, labels="AUTO", label_fontface="plain",
   label_x=0, label_y=0, hjust=-0.5, vjust=-0.5))
dev.off()

###overlap across conditions
### read data
#fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
#res <- read_rds(fn)%>%
#       filter(qval<0.1,abs(beta)>0.5)%>%
#       drop_na(beta)%>%
#       mutate(direction=ifelse(beta>0, "1", "2"), comb=paste(MCls, contrast, sep="_"))
#comb <- sort(unique(res$comb))
#ncomb <- length(comb)
#tmp <- NULL
#for ( i in 1:(ncomb-1)){
####
#  tmp0 <- map_dfr((i+1):ncomb, function(j){
#     ii <- comb[i]
#     jj <- comb[j]
#     xi <- res%>%filter(comb==ii)%>%dplyr::pull(gene)
#     xj <- res%>%filter(comb==jj)%>%dplyr::pull(gene)
#     n12 <- length(intersect(xi,xj))
#     dd <- data.frame(x=ii,y=jj, n12=n12)
#  })
#  tmp <- rbind(tmp, tmp0)
#}
#
####
#col2 <- rev(c("#f46d43", "#fdae61", "#fee08b"))
#mycol <- colorRampPalette(col2)(85)
#mycol2 <- c(rep("#FFFFC2",10), mycol, rep("#d53e4f",5))
# 
#fig1 <- ggplot(tmp, aes(x=x,y=y))+
#        geom_tile(aes(fill=n12))+
#        geom_text(aes(x=x, y=y, label=n12),size=3)+
#        scale_y_discrete(limits=rev)+
#        scale_fill_gradientn(colours=mycol)+
#        theme_void()+
#        theme(axis.title=element_blank(),
#              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
#              axis.text.y=element_text(hjust=1, vjust=0, size=10),
#              legend.title=element_blank(),
#              legend.position=c(0.8,0.8))
####
#figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure3.3.1_overlap.png"
#png(figfn, width=1000, height=1000,res=150)
#print(fig1)
#dev.off()



########################################################################################
### (5). scatterplots of beta(differential mean) between LPS/PHA and LPS-DEX/PHA-DEX ###
########################################################################################

rm(list=ls())

### label function        
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  eq 
}

##
xFun <- function(dx,a=0.5){
min1 <- min(dx$beta.x)
max2 <- max(dx$beta.x)
R <- max2-min1
xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
min1 <- min(dx$beta.y)
max2 <- max(dx$beta.y)
R <- max2-min1
ypos <- min1+a*R
}
  

    
### Read data

##load("./6_DEG.CelltypeNew_output/Filter2/Sigs.gene.DEG.RData")

fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%mutate(rn2=paste(MCls, gene, sep="_"))%>%drop_na(beta,qval)#%>%filter(gene%in%sigs)

             
### (1), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH  

dfa <- res%>%filter(contrast=="LPS")    
dfb <- res%>%filter(contrast=="LPS-DEX")%>%dplyr::select(rn2, beta, pval, qval)
       
df1 <- dfa%>%inner_join(dfb, by="rn2")

anno_df1 <- df1%>%group_by(MCls)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x,a=0.7)),
          ypos=map_dbl(data,~yFun(.x,a=1)))%>%
   dplyr::select(-data,-corr)
     
fig1 <- ggplot(df1, aes(x=beta.x, y=beta.y))+
   rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
   facet_wrap(~MCls, nrow=2, scales="free")+         
   scale_x_continuous("LPS effect on gene expression", expand=expansion(mult=0.1))+
   scale_y_continuous("LPS+DEX effect on gene expression", expand=expansion(mult=0.1))+
   theme_bw()+
   theme(strip.background=element_blank(),
         axis.title=element_text(size=12))
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure4.1_LPS.png"
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig1)
dev.off()

### (2), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
dfa <- res%>%filter(contrast=="PHA")    
dfb <- res%>%filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, pval, qval)
       
df2 <- dfa%>%inner_join(dfb,by="rn2")

anno_df2 <- df2%>%group_by(MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.7)),
          ypos=map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

fig2 <- ggplot(df2, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
    facet_wrap(~MCls, nrow=2, scales="free")+
    scale_x_continuous("PHA effect on gene expression", expand=expansion(mult=0.1))+
    scale_y_continuous("PHA+DEX effect on gene expression", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.background=element_blank(),
          axis.title=element_text(size=12))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure4.2_PHA.png"
png(filename=figfn, width=500, height=500, res=120)  
print(fig2)
dev.off()

###
###
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure4.5_comb.pdf"
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig2, nrow=1, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()

### 3, LPS and PHA
## dfa <- res%>%filter(contrast=="LPS")    
## dfb <- res%>%filter(contrast=="PHA")%>%dplyr::select(rn2, beta, pval, qval)
       
## df3 <- dfa%>%inner_join(dfb, by="rn2")

## anno_df3 <- df3%>%group_by(MCls)%>%
##            nest()%>%
##            mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##                   eq=map(corr,feq),
##                   r2=map_dbl(corr,~(.x)$estimate),
##                   xpos=map_dbl(data,~xFun(.x,a=0.7)),
##                   ypos=map_dbl(data,~yFun(.x,a=1)))%>%
##            dplyr::select(-data,-corr)
     
## fig3 <- ggplot(df3, aes(x=beta.x, y=beta.y))+
##         geom_point(size=0.3, color="grey50")+ 
##         geom_text(data=anno_df3, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
##         facet_wrap(~MCls, nrow=2, scales="free")+         
##         scale_x_continuous("LPS effect size on expression", expand=expansion(mult=0.1))+
##         scale_y_continuous("PHA effect size on expression", expand=expansion(mult=0.1))+
##         theme_bw()+
##         theme(strip.background=element_blank(),
##               axis.title=element_text(size=10))
## fig3 <- fig3+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure4.3_LPS-PHA.png"
## png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
## print(fig3)
## dev.off()


## ### 4, LPS-DEX and PHA-DEX
## dfa <- res%>%filter(contrast=="LPS-DEX")    
## dfb <- res%>%filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, pval, qval)
       
## df4 <- dfa%>%inner_join(dfb, by="rn2")

## anno_df4 <- df4%>%group_by(MCls)%>%
##            nest()%>%
##            mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##                   eq=map(corr,feq),
##                   r2=map_dbl(corr,~(.x)$estimate),
##                   xpos=map_dbl(data,~xFun(.x,a=0.7)),
##                   ypos=map_dbl(data,~yFun(.x,a=1)))%>%
##            dplyr::select(-data,-corr)
     
## fig4 <- ggplot(df4, aes(x=beta.x, y=beta.y))+
##         geom_point(size=0.3, color="grey50")+ 
##         geom_text(data=anno_df4, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
##         facet_wrap(~MCls, nrow=2, scales="free")+         
##         scale_x_continuous("LPS+DEX effect size on expression", expand=expansion(mult=0.1))+
##         scale_y_continuous("PHA+DEX effect size on expression", expand=expansion(mult=0.1))+
##         theme_bw()+
##         theme(strip.background=element_blank(),
##               axis.title=element_text(size=10))
## fig4 <- fig4+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure4.4_DEX.png"
## png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
## print(fig4)
## dev.off()



##
##




######################
### 6, upset plots ###
######################
## fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
## res <- res%>%
##         mutate(direction=ifelse(beta>0, "1", "2"),
##                comb=paste(MCls, contrast, sep="_"))
## comb <- sort(unique(res$comb))
## DEG <- map(comb, function(ii){
##    gene <- res%>%filter(comb==ii)%>%dplyr::pull(gene) 
##    gene
## })
## names(DEG) <- comb
## df2 <- list_to_matrix(DEG)

## df3 <- df2[,grepl("LPS$",comb)]
## mat <- make_comb_mat(df3)
## m2 <- mat[comb_degree(mat)>0]
## fig0 <- UpSet(m2, set_order=colnames(df3),comb_order=order(-comb_size(m2)), right_annotation=NULL, row_names_side="left")
## #              left_annotation=rowAnnotation("Set size"=anno_barplot(set_size(mat),
## #                                                       axis_param=list(direction="reverse"),
## #                                                       border=F,
## #                                                       width=unit(2,"cm"),
## #                                                       gp=gpar(fill="black"))),
## #              right_annotation=NULL,
## #              row_names_side=left)

## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure5.4_PHA-DEX.upset.png"
## png(figfn, width=1000, height=500, res=120)
## print(fig0)
## dev.off()


########################################
### 7, DEGs shared across cell types ###
########################################
library(ComplexHeatmap)
rm(list=ls())
### function
myData <- function(oneContrast, res){
   MCls <- sort(unique(res$MCls))
   DEG <- map(MCls,function(oneMCl){
      gene <- res%>%filter(MCls==oneMCl, contrast==oneContrast)%>%dplyr::pull(gene)
      gene
   })
   names(DEG) <- MCls
   df2 <- list_to_matrix(DEG)
   mat <- make_comb_mat(df2)
   dd <- data.frame(degree=comb_degree(mat), ngene=comb_size(mat))
   dd2 <- dd%>%
          group_by(degree)%>%
          summarise(ngene=sum(ngene), .groups="drop")%>%
          mutate(contrast=oneContrast)%>%filter(degree!=0)   
   dd2 <- dd2%>%mutate(prop=ngene/sum(ngene))
} 
##    
mypos <- function(d){ 
   ypos <- lapply(sort(unique(d$contrast)), function(ii){
      prop <- d%>%filter(contrast==ii)%>%dplyr::pull(prop)
      m1 <- 1-cumsum(prop)
      m2 <- c(1,m1[1:3])
      ytmp <- 0.5*(m1+m2)
      ytmp
   })
   ypos <- do.call(c,ypos)
}

### Read data
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
contrast <- sort(unique(res$contrast))
plotData <- map_dfr(contrast, myData, res=res)

plotData <- plotData%>%mutate(x=round(prop*100,digits=1))
plotData[4,"x"] <- 1.5
plotData[11,"x"] <- 8.3
plotData <- plotData%>%mutate(symbol=paste(x,"%",sep=""))
plotData$ypos <- mypos(plotData) 
##

contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
##
colw <- sapply(c(0.8, 0.6, 0.4, 0.2),function(x) colorspace::lighten("#800080", x))
names(colw) <- c("1", "2", "3", "4") 

fig0 <- ggplot(plotData, aes(x=factor(contrast), y=prop))+
        geom_bar(stat="identity", aes(fill=factor(degree)))+     
        geom_text(aes(label=symbol, y=ypos), size=3)+
        scale_fill_manual(values=colw)+
        scale_x_discrete(labels=lab2)+
        theme(legend.title=element_blank(),
              axis.title=element_blank(),
              axis.text.x=element_text(size=12, vjust=1),
              axis.text.y=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank(),
              panel.background=element_blank())

###              
## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure6.png"
## png(figfn, width=400, height=600, res=120)
figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure6.pdf"
pdf(figfn, width=4, height=6)
print(fig0)  
dev.off()              



############################################
### 6, Examples GSE enrichment analysis  ###
############################################
###
#lab1 <- c("CTRL"="CTRL", 
#          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
#          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
##col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
#
#### expression
#### fig 1
#if(FALSE){
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#count <- colSums(YtX)
#count_after <- median(count)
#count <- count/count_after
####
#bti2 <- colnames(YtX)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                  sampleID=cvt0[,3], Batch2=cvt0[,4])
#
##gene0 <- "ENSG00000126709"
##symbol0 <- "IFI6"
#gene0 <- "ENSG00000161980"
#symbol0 <- "POLR3K"                                                 
#cvt$y <- YtX[grepl(gene0,rownames(YtX)),]/count 
#
#d1 <- cvt%>%drop_na(y)%>%filter(y>0)
#fig1 <- ggplot(d1%>%filter(MCls=="Monocyte"), aes(x=treats, y=log2(y), colour=treats))+
#        geom_boxplot()+
#        ylab(bquote(log[2]~"(Expression)"))+
#        #scale_fill_manual("", values=col1, labels=lab1)+
#        scale_colour_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("", labels=lab1)+
#        ggtitle(bquote(~italic(.(symbol0))~" expressed in Monocyte cell"))+
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
#              axis.title.y=element_text(size=10),
#              plot.title=element_text(hjust=0.5, size=10),
#              legend.position="none")
#              
#figfn <- paste("./6_DEG.CelltypeNew_output/Figure4.2_", gene0, ".Mono.Boxplot.png", sep="")
#png(figfn, width=400, height=600, res=120)
#print(fig1)  
#dev.off()
#
#}        
#



###############################
### 7, Batch and chemistry ####
###############################
#if(FALSE){
#
#
#load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
#DEG <- sigs
# 
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX") 
###
#dx <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./6_DEG.CelltypeNew_output/1_Batch_", oneMCl, "_DESeq.results.txt", sep="")
#   res <- read.table(fn, header=T)%>%mutate(MCls=oneMCl)%>%filter(gene%in%DEG)
#   res
#})
#
#dx2 <- dx%>%
#      group_by(MCls, contrast, Batch2, sep="_")%>%
#      nest()%>%
#      mutate(rn=paste(MCls, contrast, Batch2, sep="_"),
#             beta=map(data,~(.x)$estimate))%>%
#      dplyr::select(-data)
#      
#TMP <- map_dfc(dx2$beta, ~(.x))
#names(TMP) <- dx2$rn
#comb <- names(TMP)
#
#ii <- rowSums(is.na(TMP))
#TMP0 <- TMP[ii==0,]
#y <- do.call(c, TMP0)
#y0 <- y[abs(y)<2] #99% percent quantile(abs(y),probs=0.99)
#mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
#names(mybreaks) <- NULL
#
####colors
#x <- str_split(names(TMP0), "_", simplify=T)
#tmp_column <- data.frame(MCls=x[,1], treats=x[,2], Batch=x[,3])
#rownames(tmp_column) <- comb
#tmp_colors <- list(MCls=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
#                          "NKcell"="#377eb8", "Tcell"="#e41a1c"), 
#                   treats=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#                            "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"),
#                   Batch=c("SCAIP1"="#e41a1c", "SCAIP2"="#377eb8",
#                           "SCAIP3"="#4daf4a", "SCAIP4"="#984ea3",
#                           "SCAIP5"="#ff7f00", "SCAIP6"="#ffff33")) #brewer.pal(4,"Set1")
#
#mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
#fig1 <- pheatmap(TMP0, col=mycol, breaks=mybreaks, 
#         scale="none",
#         border_color="NA",
#         cluster_rows=T, cluster_cols=T, 
#         annotation_col=tmp_column,
#         annotation_colors=tmp_colors,
#         show_colnames=T, show_rownames=F,
#         fontsize_row=6,
#         na_col="white")
#
#figfn <- "./6_DEG.CelltypeNew_output/Figure5.1_heatmap.beta.png"
#png(figfn, width=2000, height=1200,res=150)
#fig1
#dev.off()
#}  


### (2), Barplots show DGV, up and down separately, meanwhile with denotions of significance 
#if(FALSE){
#
#cat("(2).","DGV up and down", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#
#tmp <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./6_DEG.CelltypeNew_output/1_Batch_", oneMCl, "_meta.txt", sep="")
#   res <- read.table(file=fn, header=T)
#   res0 <- res%>%
#           mutate(MCls=oneMCl)%>%
#           filter(qval<0.1, abs(beta)>0.5, !is.na(qval))
#   res0 
#})
#
### up and down DGV
#sigs <- tmp%>%
#        mutate(direction=ifelse(beta>0, "1", "2"))%>%
#        group_by(contrast, MCls, direction)%>%
#        summarise(ngene=n())
#
###fun-1
#fmod <- function(subdf) {
#   res <- binom.test(subdf$ngene, 0.5, alternative="two.sided")
#   res$p.value
#}    
###fun-2       
#flab <- function(p){
#    if (p<=0.001){
#       mylab <- "***"
#    }else if(p<=0.01){
#       mylab <- "**"
#    }else if(p<=0.05){
#       mylab <- "*"
#    }else{
#       mylab <- "NS"
#    }
#    mylab
#}
#anno_df <- sigs%>%
#              group_by(contrast, MCls)%>%
#              nest()%>%
#              mutate(pval=map_dbl(data,fmod),
#                     y=map_dbl(data,~max((.x)$ngene)),
#                     label=map_chr(pval,flab))%>%
#              unnest()%>%
#              distinct(contrast, MCls,.keep_all=T)
#xpos <- c("LPS"=0.8, "LPS-DEX"=1.8, "PHA"=2.8, "PHA-DEX"=3.8)
#xmin <- xpos[anno_df$contrast]
#anno_df$xmin <- xmin
#anno_df <- anno_df%>%
#           mutate(xmax=xmin+0.4, y1=y+50)
#
#
#
#anno_df$group <- 1:nrow(anno_df)
#anno_df <- anno_df%>%filter(pval<0.05)
#
#cols <- c("1"="red","2"="blue")
#mylab <- c("1"="up","2"="down")
#               
#fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=direction))+geom_bar(stat="identity",position="dodge")+  #"stack"
#        scale_fill_manual(values=cols,labels=mylab)+
#        geom_signif(data=anno_df, 
#                    aes(xmin=xmin, xmax=xmax, annotations=label, y_position=y1, group=group),
#                    vjust=0.1, tip_length=0.05, manual=T)+   ##need to add group in aes
#        facet_wrap(~factor(MCls),ncol=2)+
#        xlab("")+ylab("No. DEG")+ylim(0,1800)+
#        theme_bw()+
#        theme(strip.background=element_blank(),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9),
#              strip.text.x=element_text(size=12))
####
#figfn <- "./6_DEG.CelltypeNew_output/Figure2.2_DEG.barplot2.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
#} ### End, (2)

##


  




