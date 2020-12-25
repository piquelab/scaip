###
rm(list=ls())
source("./Bin/LibraryPackage.R")
outdir <- "./6_DEG.CelltypeNew_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#####################################################
###    6. DEG analysis for each Cell type         ###
###    Batch one by one then meta analysis        ###
###    11-26-2020, By Julong Wei, last modified   ###
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
if (FALSE){

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
annoSel <- anno%>%filter(uns, rnz>0, chr%in%autosome)#, grepl("protein_coding", biotype))
### Out of 17,697 protein coding gene, 15,770 with reads>20
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
opfn <- "./6_DEG.CelltypeNew_output/YtX.comb.RData"
save(YtX, file=opfn)
###

### average 
#ncell <- colSums(X)
#YtXave <- sweep(YtX, 2, ncell, FUN="/")
#opfn <- "./6_DEG_CelltypeNew_output/YtX.ave.RData"
#save(YtXave, file=opfn)


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

} ###1, End



##################################################
### 2, Differential gene expression analysis   ###
##################################################
if (FALSE){
cat("2.", "Differential gene expression analysis", "\n")
library("BiocParallel")  
library(DESeq2)
register(MulticoreParam(1))


### (1), prepare data
load("./6_DEG.CelltypeNew_output/YtX.comb.RData")  
rownames(YtX) <- gsub("\\.[0-9]*", "", rownames(YtX)) ##gene id, like ENSG00000187634 

grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
anno <- data.frame(ensgene=rownames(YtX))%>%mutate(rnz=rowSums(YtX))%>%left_join(grch38_unq, by="ensgene")
geneSel <- anno%>%filter(rnz>20,grepl("protein_coding", biotype))%>%dplyr::pull(ensgene)

load("./6_DEG.CelltypeNew_output/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)
rnz <- rowSums(YtX) 
YtX_sel <- YtX[geneSel, ddx$bti] ###rnz>20

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

opfn <- "./6_DEG.CelltypeNew_output/1_DESeq.results.rds"
write_rds(res,opfn)

} ###End, 2.1
### (3) plots

##########################
### 2.2, meta analysis ###
##########################
###
if(FALSE){
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

fn <- "./6_DEG.CelltypeNew_output/1_DESeq.results.rds"
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
opfn <- "./6_DEG.CelltypeNew_output/3_Batch1456.meta.rds"
write_rds(res3, opfn)

}   
###2.2, End



############################################
### 3, Summary of results of DE analysis ###
############################################
if (FALSE){
cat("3.", "Summary results", "\n")

##############################
### 3.1, figures, MA plots ###
##############################
cat("3.1.", "Figure MA plots", "\n")
figfn <- "./6_DEG.CelltypeNew_output/Figure1.1.MA.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300) 
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
res <- read_rds(fn)%>%mutate(color=ifelse(qval<0.1, T, F))%>%
          dplyr::rename(baseMean=baseMeanHat)%>%drop_na(beta)

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){   
##1  
   res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS")%>%dplyr::select(baseMean, beta, color, qval, pval)  
   print(plotMA(res1[,1:3], colLine="NA", main="LPS-ctrl vs CTRL", cex.main=1, cex.axis=0.8, cex.lab=1))

##2
   res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")%>%dplyr::select(baseMean, beta, color, qval, pval)
   print(plotMA(res2[,1:3], colLine="NA", main="LPS-DEX vs LPS-ctrl", cex.main=1, cex.axis=0.8, cex.lab=1))
##3
   res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")%>%dplyr::select(baseMean, beta, color, qval, pval)
   print(plotMA(res3[,1:3], colLine="NA", main="PHA-ctrl vs CTRL", cex.main=1, cex.axis=0.8, cex.lab=1)) 
##4
   res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")%>%dplyr::select(baseMean, beta, color, qval, pval)
   print(plotMA(res4[,1:3], colLine="NA", main="PHA-DEX vs PHA-ctrl", cex.main=1, cex.axis=0.8, cex.lab=1))
   
  print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue") )        
}
dev.off()


#####################
### 3.2, qq plots ###
#####################
cat("3.2.", "qq plots", "\n")
figfn <- "./6_DEG.CelltypeNew_output/Figure1.2.qq.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
res <- read_rds(fn)%>%
       mutate(color=ifelse(qval<0.1, T, F))%>%
       dplyr::rename(baseMean=baseMeanHat)%>%drop_na(beta)

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){   
##1  
   res1 <- res%>%filter(MCls==oneMCl, contrast=="LPS")%>%dplyr::select(beta, color, qval, pval)  
   print( qq(res1$pval, main="LPS-ctrl vs CTRL", cex.main=1, cex.axis=0.8, cex.lab=1))

##2
   res2 <- res%>%filter(MCls==oneMCl, contrast=="LPS-DEX")%>%dplyr::select(beta, color, qval, pval)
   print( qq(res2$pval, main="LPS-DEX vs LPS-ctrl", cex.main=1, cex.axis=0.8, cex.lab=1))
##3
   res3 <- res%>%filter(MCls==oneMCl, contrast=="PHA")%>%dplyr::select(beta, color, qval, pval)
   print( qq(res3$pval, main="PHA-ctrl vs CTRL", cex.main=1, cex.axis=0.8, cex.lab=1)) 
##4
   res4 <- res%>%filter(MCls==oneMCl, contrast=="PHA-DEX")%>%dplyr::select(beta, color, qval, pval)
   print( qq(res4$pval, main="PHA-DEX vs PHA-ctrl", cex.main=1, cex.axis=0.8, cex.lab=1))
   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))         
}
dev.off()

###
### 3.3, hist distribution of differential effects
fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
dx <- read_rds(fn)%>%drop_na(beta)
fig0 <- ggplot(dx, aes(x=beta))+
     geom_histogram(fill="grey70", color="grey20")+
     xlab(bquote("Effective size"~"("~beta~")"))+
     facet_grid(MCls~contrast, scales="free")+
     theme_bw()+
     theme(strip.background=element_blank())

figfn <- "./6_DEG.CelltypeNew_output/Figure1.3.hist.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig0)
dev.off()
} ###3, End


#########################
### 3.3 boxplot show  ###
#########################
if(FALSE){

load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
rn <- rownames(YtX) 
rownames(YtX) <- gsub("\\..*", "", rn)
 
YtX0 <- YtX[sigs,]
dx <- melt(YtX0)
cvt <- str_split(dx$X2, "_", simplify=T)

dx <- dx%>%
      mutate(MCls=cvt[,1], treats=gsub("-EtOH","", cvt[,2]),sampleID=cvt[,3])%>%
      dplyr::rename(ensgene=X1)

dx2 <- dx%>%group_by(ensgene, MCls, treats)%>%summarise(y=mean(value))

cols1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4") 
cols2 <- c("Tcell"="#e41a1c", "NKcell"="#377eb8", 
          "Bcell"="#4daf4a", "Monocyte"="#984ea3")
          
fig1 <- ggplot(dx2)+
        geom_boxplot(aes(x=MCls, y=log2(y), color=MCls))+
        scale_color_manual(values=cols2)+
        ylab(bquote(log[2]~"(GE)"))+ 
        theme_bw()+
        theme(axis.title.x=element_blank(),
              legend.position="none")
        
fig2 <- ggplot(dx2)+
        geom_boxplot(aes(x=treats, y=log2(y), color=treats))+
        scale_color_manual(values=cols1)+
        ylab(bquote(log[2]~"(GE)"))+ 
        theme_bw()+
        theme(axis.title.x=element_blank(),
              legend.position="none")
        
figfn <- "./6_DEG.CelltypeNew_output/Figure1.4.GE.boxplot.png"
png(figfn, width=650, height=300, res=110)
print(plot_grid(fig1, fig2, ncol=2, align="hv",
                labels="AUTO", label_fontface="plain"))
dev.off() 

}


########################################     
### 4, extract significant genes DEG ###
########################################
if (FALSE){

###4.1, Table of DEG

###
fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1, abs(beta)>0.5)%>%drop_na(beta)
sigs <- unique(res$gene)
save(sigs, file="./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")

} ###

#######################################################################
### 4.2 show barplot DEG  from different contrast across cell types ###
#######################################################################
if (FALSE){

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

### (1). Show number of DEG with different contrast across cell type
fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)

sigs <- res%>%
        group_by(contrast, MCls)%>%
        summarise(ngene=n())
x <- res%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))


#cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
#          "NKcell"="#377eb8", "Tcell"="#e41a1c")
cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=MCls))+geom_bar(stat="identity", position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(name="", labels=lab2)+ylab("No. DEG")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./6_DEG.CelltypeNew_output/Figure2.1_DEG.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
} ### End

### (3), barplots of DEG, up and down with light and deep colors, ***
if(FALSE){

cat("(3).","DGV up and down", "\n")

###facet by MCls
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2,col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)

## up and down DGE
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n())
###        
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

figfn <- "./6_DEG.CelltypeNew_output/Figure2.3_DEG.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### facet by contrast
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

figfn <- "./6_DEG.CelltypeNew_output/Figure2.4_DEG.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

### facet by MCls and up above axis and down below axis
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2,-ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-1800,1800),5)

facetlab <- as_labeller(lab2treat)
fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2, fill=comb))+
        geom_bar(stat="identity")+
        scale_fill_manual(values=col2comb, labels="")+
        geom_hline(yintercept=0, color="grey60")+
        geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
                  vjust=ifelse(direction==2, 1.2, -0.2)), size=3)+ #
        scale_y_continuous("DEG", breaks=breaks_value, limits=c(-2000,2000),labels=abs(breaks_value))+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

figfn <- "./6_DEG.CelltypeNew_output/Figure2.5_DEG.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

}



############################################
### 5,   heatmap for 16 conditions       ###
###  correlation between 16 conditions #####
##########################################
if (TRUE){

load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
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
fn <- "./6_DEG.CelltypeNew_output/3_Batch1456.meta.rds"
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

figfn <- "./6_DEG.CelltypeNew_output/Figure3.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "Tcell_LPS+DEX", "Tcell_PHA+DEX", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
               "Tcell_LPS", "Tcell_PHA", "NKcell_LPS", "NKcell_PHA") 

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

figfn <- "./6_DEG.CelltypeNew_output/Figure3.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()



#figfn <- "./6_DEG.CelltypeNew_output/tmp1/Figure4.2_corr.beta.png"
#png(figfn, width=1000, height=1000, res=180)
#print(corrplot(corr, method="color", order="hclust", hclust.method="complete", col=mycol,
#         tl.col="black", tl.cex=0.8, outline=F, diag=T))
#dev.off()

} ##End, 5



############################################
### 6, Examples GSE enrichment analysis  ###
############################################
###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
#          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

### expression
### fig 1
if(FALSE){
load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
count <- colSums(YtX)
count_after <- median(count)
count <- count/count_after
###
bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                  sampleID=cvt0[,3], Batch2=cvt0[,4])

#gene0 <- "ENSG00000126709"
#symbol0 <- "IFI6"
gene0 <- "ENSG00000161980"
symbol0 <- "POLR3K"                                                 
cvt$y <- YtX[grepl(gene0,rownames(YtX)),]/count 

d1 <- cvt%>%drop_na(y)%>%filter(y>0)
fig1 <- ggplot(d1%>%filter(MCls=="Monocyte"), aes(x=treats, y=log2(y), colour=treats))+
        geom_boxplot()+
        ylab(bquote(log[2]~"(Expression)"))+
        #scale_fill_manual("", values=col1, labels=lab1)+
        scale_colour_manual("", values=col1, labels=lab1)+
        scale_x_discrete("", labels=lab1)+
        ggtitle(bquote(~italic(.(symbol0))~" expressed in Monocyte cell"))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=10),
              legend.position="none")
              
figfn <- paste("./6_DEG.CelltypeNew_output/Figure4.2_", gene0, ".Mono.Boxplot.png", sep="")
png(figfn, width=400, height=600, res=120)
print(fig1)  
dev.off()

}        




###############################
### 7, Batch and chemistry ####
###############################
if(FALSE){


load("./6_DEG.CelltypeNew_output/Sigs.gene.DEG.RData")
DEG <- sigs
 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX") 
##
dx <- map_dfr(MCls, function(oneMCl){
   fn <- paste("./6_DEG.CelltypeNew_output/1_Batch_", oneMCl, "_DESeq.results.txt", sep="")
   res <- read.table(fn, header=T)%>%mutate(MCls=oneMCl)%>%filter(gene%in%DEG)
   res
})

dx2 <- dx%>%
      group_by(MCls, contrast, Batch2, sep="_")%>%
      nest()%>%
      mutate(rn=paste(MCls, contrast, Batch2, sep="_"),
             beta=map(data,~(.x)$estimate))%>%
      dplyr::select(-data)
      
TMP <- map_dfc(dx2$beta, ~(.x))
names(TMP) <- dx2$rn
comb <- names(TMP)

ii <- rowSums(is.na(TMP))
TMP0 <- TMP[ii==0,]
y <- do.call(c, TMP0)
y0 <- y[abs(y)<2] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL

###colors
x <- str_split(names(TMP0), "_", simplify=T)
tmp_column <- data.frame(MCls=x[,1], treats=x[,2], Batch=x[,3])
rownames(tmp_column) <- comb
tmp_colors <- list(MCls=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                          "NKcell"="#377eb8", "Tcell"="#e41a1c"), 
                   treats=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                            "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"),
                   Batch=c("SCAIP1"="#e41a1c", "SCAIP2"="#377eb8",
                           "SCAIP3"="#4daf4a", "SCAIP4"="#984ea3",
                           "SCAIP5"="#ff7f00", "SCAIP6"="#ffff33")) #brewer.pal(4,"Set1")

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
fig1 <- pheatmap(TMP0, col=mycol, breaks=mybreaks, 
         scale="none",
         border_color="NA",
         cluster_rows=T, cluster_cols=T, 
         annotation_col=tmp_column,
         annotation_colors=tmp_colors,
         show_colnames=T, show_rownames=F,
         fontsize_row=6,
         na_col="white")

figfn <- "./6_DEG.CelltypeNew_output/Figure5.1_heatmap.beta.png"
png(figfn, width=2000, height=1200,res=150)
fig1
dev.off()
}  


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


  




