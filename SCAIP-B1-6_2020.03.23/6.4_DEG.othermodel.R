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
## library(Seurat)
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
outdir <- "./6_DEG.CelltypeNew_output/Response_reviews/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#############################################################################
###    6. DEG analysis for each Cell type                                 ###
###    Batch one by one then meta analysis, number of cells as covariates ###
###    5-25-2022, By Julong Wei, last modified                            ###
#############################################################################

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


## ##################################################
## ### 1, prepare data, pseudo-bulk seq data, YtX ###
## ##################################################

## cat("1.", "Generate pseudo-bulk RNA seq data", "\n")
## sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
## count <- sc@assays$RNA@counts         

## ##gene:64428, symbol:58149, ensgene:63697,
## grch38_unq <- grch38%>%
##               distinct(ensgene,.keep_all=T)%>%
##               dplyr::select(ensgene, symbol, chr, start, end, biotype)
## anno <- data.frame(rn=rownames(count))%>%
##         mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn),
##                ensgene2=gsub("[SU]-","",rn),
##                uns=grepl("S-",rn), rnz=rowSums(count))%>%
##         left_join(grch38_unq, by="ensgene")
        
## #write.table(anno, file="./5_Identify_CelltypeOld_output/gene.infor", row.names=F, col.names=T, quote=F, sep="\t")        
## #ii <- duplicated(grch38$ensgene)
## #gene0 <- unique(grch38[ii,][["ensgene"]])
## #nn <- sapply(gene0, function(one){
## #x <- grch38[grch38$ensgene==one,]%>%select(ensgene, symbol,chr,start,end)
## #comb <- unique(paste(x$symbol,x$chr,sep="_"))
## #length(comb)
## #})

## autosome <- as.character(1:22)        
## annoSel <- anno%>%filter(uns, rnz>0, chr%in%autosome)
## ### Out of 42,554 genes, 30,972 whole genes>20           
## Ys <- count[annoSel$rn,]
## rownames(Ys) <- annoSel$ensgene2 # rownames like ENSG00000187634.12

## meta <- sc@meta.data
## meta <- meta%>%
##        mutate(bti=paste(MCls, treats, BEST.GUESS, BATCH, sep="_"))%>%
##        dplyr::select(NEW_BARCODE, bti)

       
## ### summary the number of cells in 1536 combinations
## dd <- meta %>%group_by(bti)%>%summarise(ncell=n())
## save(dd, file="./6_DEG.CelltypeNew_output/0_ncell.RData")
## #p <- ggplot(dd, aes(x=ncell))+
## #     geom_histogram(alpha=0.5, color=NA, fill="blue")+xlab("No.cell")+
## #     theme_bw()+
## #     theme(plot.title=element_text(hjust=0.5))
## #png("./6_DEG_CelltypeNew_output/Figure0.ncell.png", width=200, height=200, res=100)
## #p
## #dev.off()
## ### X, obs, bti
       
## bti <- factor(meta$bti)       
## X <- model.matrix(~0+bti)
## YtX <- Ys %*% X
## YtX <- as.matrix(YtX)
## colnames(YtX) <- gsub("^bti", "", colnames(YtX))
## opfn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
## save(YtX, file=opfn)
## ###

### average 
#ncell <- colSums(X)
#YtXave <- sweep(YtX, 2, ncell, FUN="/")
#opfn <- "./6_DEG_CelltypeNew_output/YtX.ave.RData"
#save(YtXave, file=opfn)

### protein coding genes, rnz>20, combination with >20 cells
  
## grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
## anno <- data.frame(ensgene=gsub("\\..*", "", rownames(YtX)),ensgene2=rownames(YtX))%>%
##         mutate(rnz=rowSums(YtX))%>%
##         left_join(grch38_unq, by="ensgene")
## geneSel <- anno%>%filter(rnz>20,grepl("protein_coding", biotype))%>%dplyr::pull(ensgene2)

## load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
## ddx <- dd%>%filter(ncell>20)

## YtX_sel <- YtX[geneSel, ddx$bti] ###rnz>20
## opfn <- "./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
## save(YtX_sel, file=opfn)


##################################################
### 2, Differential gene expression analysis   ###
##################################################

library("BiocParallel")  
## library(DESeq2)
register(MulticoreParam(1))

load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)

### (1), prepare data
load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")  
rownames(YtX_sel) <- gsub("\\.[0-9]*", "", rownames(YtX_sel)) ##gene id, like ENSG00000187634 
###
bti2 <- colnames(YtX_sel)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))%>%left_join(ddx, by="bti")
cvt <- cvt%>%mutate(ncell=as.numeric(ncell))

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
   
res <- map_dfr(comb, function(oneX){   
   time0 <- Sys.time()
   cvt0 <- cvt%>% filter(comb==oneX) 
   oneMCl <- cvt0$MCls[1]
   oneBatch <- cvt0$Batch2[1]
   YtX0 <- YtX_sel[,cvt0$bti]
   dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treats+ncell)
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

opfn <- paste(outdir, "1_DESeq.results.rds", sep="")
write_rds(res,opfn)


##########################
### 2.2, meta analysis ###
##########################


#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#res <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./6_DEG.CelltypeNew_output/tmp1/3_Batch1456_", oneMCl, "_meta.txt", sep="")
#   tmp <- read.table(file=fn, header=T)%>%mutate(MCls=oneMCl)
#   tmp
#})
#opfn <- "./6_DEG.CelltypeNew_output/tmp1/3_Batch1456.meta.rds"
#write_rds(res, opfn)
###

fn <- paste(outdir, "1_DESeq.results.rds", sep="")
res <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))

## %>%
## filter(Batch2%in%c("SCAIP1", "SCAIP4", "SCAIP5", "SCAIP6")) 

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
opfn <- paste(outdir, "2_meta.rds", sep="")
write_rds(res3, opfn)



####
####
df1 <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
df1 <- df1%>%mutate(zscore=beta/stderr)%>%dplyr::select(MCls, contrast, rn, zscore)
###
df2 <- read_rds("./6_DEG.CelltypeNew_output/Response_reviews/2_meta.rds")
df2 <- df2%>%mutate(zscore=beta/stderr, rn=gsub("\\..*", "", rn))%>%dplyr::select(rn, zscore)

dfcomb <- df1%>%inner_join(df2, by="rn")

p <- ggplot(dfcomb, aes(x=zscore.x, y=zscore.y))+
   rasterise(geom_point(size=0.3, colour="grey30"), dpi=300)+
   geom_abline(colour="red")+
   facet_grid(MCls~contrast,
              labeller=labeller(contrast=c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                                           "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
   xlab("z score w/o number of cells")+
   ylab("z score w/. number of cells")+
   theme_bw()

figfn <- paste(outdir, "Figure2_scatter.plot.png", sep="")
png(figfn, width=600, height=600, res=120)
p
dev.off()



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



  




