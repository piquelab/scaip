##
library(tidyverse)
library(parallel)
library(purrr)
#library(furrr)
#library(reshape)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(annotables)
library(biobroom)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

##############################
### Diffferential function ###
##############################

###Fun-1, 

### myDE, batch separately
myDE <- function(y, X, gene, nna, threshold=8){

   con.ls <- list("LPS"=c("CTRL","LPS-EtOH"),
                  "LPS-DEX"=c("LPS-EtOH","LPS-DEX"),
                  "PHA"=c("CTRL", "PHA-EtOH"),
                  "PHA-DEX"=c("PHA-EtOH","PHA-DEX"))
   Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

   #y <- try(qqnorm(y,plot.it=FALSE)$x,silent=T) 
   y <- try(log2(y),silent=T)
   x1 <- X[,1]
   lm0 <- try(lm(y~0+x1),silent=T)
   
   if ( (class(y)!="try-error") & (class(lm0)!="try-error") ){
      b <- coef(lm0)
      nn <- gsub("x[12]", "", names(b))
      names(b) <- nn
      vb <- diag(vcov(lm0))
      names(vb) <- nn
         
      ## Contrast       
      dd <- lapply(Contrast,function(one){
         con0 <- con.ls[[one]] 
         if ( (all(con0%in%nn))&(all(nna[con0]>=threshold))  ){
            bhat <- b[con0[2]]-b[con0[1]]
            sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
            z <- bhat/sdhat
            p <- 2*pnorm(-abs(z))
            dd1 <- data.frame(gene=gene, beta=bhat, stderr=sdhat, pval=p, contrast=one)            
         }else{
            dd1 <- NA
         }           
         dd1           
      })
      dd <- dd[!is.na(dd)]
      dd <- do.call(rbind,dd)             
   }else{
      dd <- NA
   } ## End if try-error
   
} ###

###Batch together
myDE2 <- function(y, X, gene, nna, threshold=8){

   con.ls <- list("LPS"=c("CTRL","LPS-EtOH"),
                  "LPS-DEX"=c("LPS-EtOH","LPS-DEX"),
                  "PHA"=c("CTRL", "PHA-EtOH"),
                  "PHA-DEX"=c("PHA-EtOH","PHA-DEX"))
   Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

   #y <- try(qqnorm(y,plot.it=FALSE)$x,silent=T) 
   y <- try(log2(y),silent=T)
   x1 <- X[,1]       
   x2 <- X[,2]   
   lm0 <- try(lm(y~0+x1+x2),silent=T)
   
   if ( (class(y)!="try-error") & (class(lm0)!="try-error") ){
      b <- coef(lm0)
      nn <- gsub("x[12]", "", names(b))
      names(b) <- nn
      vb <- diag(vcov(lm0))
      names(vb) <- nn
         
      ## Contrast       
      dd <- lapply(Contrast,function(one){
         con0 <- con.ls[[one]] 
         nna1 <- nna[grepl(con0[1],names(nna))]
         nna2 <- nna[grepl(con0[2],names(nna))]
         if ( (all(con0%in%nn))&(any(nna1>=threshold))&(any(nna2>=threshold)) ){
            bhat <- b[con0[2]]-b[con0[1]]
            sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
            z <- bhat/sdhat
            p <- 2*pnorm(-abs(z))
            dd1 <- data.frame(gene=gene, beta=bhat, stderr=sdhat, pval=p, contrast=one)            
         }else{
            dd1 <- NA
         }           
         dd1           
      })
      dd <- dd[!is.na(dd)]
      dd <- do.call(rbind,dd)             
   }else{
      dd <- NA
   } ## End if try-error
   
} ###

### Fun-2, myMeta
myMeta <- function(dx){
   sub0 <-  !(is.na(dx[,"beta"])|is.na(dx[,"stderr"])) 
   if (any(sub0)){
      b <- dx[sub0,"beta"]
      std <- dx[sub0,"stderr"]
      ###
      vb_i <- 1/std^2
      w <- vb_i/sum(vb_i)
      bhat <- sum(w*b)
      sdhat <- sqrt(1/sum(vb_i))
      z <- bhat/sdhat
      p <- 2*pnorm(abs(z),lower.tail=F)
      tmp <- c(bhat, sdhat, p)
   }else{
      tmp <- c(NA, NA, NA)
   }
   tibble(beta=tmp[1],stderr=tmp[2],pval=tmp[3])
}

### Fun-3, myqval
myqval <- function(pval){
   qval <- pval
   ii0 <- !is.na(pval)
   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
   qval
}

### Fun-4, countNA
countNA <- function(X,y){
   cvt <- data.frame(x1=X, y=y)
   d2 <- cvt%>%group_by(x1)%>%summarise(n=sum(!is.na(y)), .groups="drop")
   n <- d2$n
   names(n) <- d2$x1
   n
}

countNA2 <- function(X,y){
   cvt <- data.frame(x1=X$x1, x2=X$x2, y=y)%>%mutate(rn=paste(x1, x2, sep="_"))
   d2 <- cvt%>%group_by(rn)%>%summarise(n=sum(!is.na(y)), .groups="drop")
   n <- d2$n
   names(n) <- d2$rn
   n
}


### 4, differential of gene mean  ###
#####################################

if(FALSE){
load("./10_RNA.Variance_output/tmp9/1.2_Sel.Bx.RData")
rn <- rownames(Bx)
rownames(Bx) <- gsub("\\.[0-9]*", "", rn)
###
bti2 <- colnames(Bx)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
comb <- unique(cvt$comb)

### estimate differetial results
cat("Differential analysis by batch, gene,", nrow(Bx), "\n")
res <- map_dfr(comb, function(oneX){
   cat(oneX,"\n")
   cvti <- cvt %>%filter(comb==oneX)
   oneComb <- unlist(strsplit(oneX, "_"))
   oneMCl <- oneComb[1]
   oneBatch <- oneComb[2]
   Bi <- Bx[,cvti$bti]
   X <- data.frame(x1=cvti$treats)
   rn <- rownames(Bi)
   ## Start loop By gene 
   TMP <- mclapply(rn, function(ii){
      y <- Bi[ii,]
      nna <- countNA(X$x1,y)
      dd <- myDE(y, X, ii, nna, threshold=3)
      dd
   }, mc.cores=1) ### End loop by gene
   ###  
   TMP <- TMP[!is.na(TMP)]
   TMP <- as.data.frame(do.call(rbind, TMP))%>%mutate(MCls=oneMCl, batch=oneBatch)
   TMP       
})

opfn <- "./10_RNA.Variance_output/tmp9/4_mu.results"
write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
} ###

###"Meta analysis"
if(FALSE){
cat("Meta analysis", "\n")
###Read data
fn <- "./10_RNA.Variance_output/tmp9/4_mu.results"
res <- read.table(fn,header=T)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
       filter(batch%in%c("SCAIP1","SCAIP4", "SCAIP5", "SCAIP6"))

### (2) Meta
res2 <- res%>%group_by(rn)%>%
        nest()%>%
        mutate(outlist=mclapply(data, myMeta, mc.cores=1))%>%
        dplyr::select(-data)%>%unnest(c(outlist))%>%as.data.frame()
cvt <- str_split(res2$rn, "_", simplify=T)
res2 <- res2%>%mutate(MCls=cvt[,1], contrast=cvt[,2], gene=cvt[,3])
          
### (3) add qvalue
res3 <- res2%>%group_by(MCls, contrast)%>%
        nest()%>%
        mutate(qval=map(data, ~myqval((.x)$pval)))%>%
        unnest(c(data,qval))%>%as.data.frame()
        
#opfn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
opfn <- "./10_RNA.Variance_output/tmp9/4_mu.meta2"
write.table(res3, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")

} ###End, Differential analysis

###
#### 5.3, estimate differetial results all the batch together
#if(FALSE){
#cat("5.3", "Read data", "\n") 
#
#load("./10_RNA.Variance_output/tmp7/1_RNA.Bx.RData")
#rn <- rownames(Bx)
#rownames(Bx) <- gsub("\\.[0-9]*", "", rn)
####
#bti2 <- colnames(Bx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
#cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
#
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#res <- map_dfr(MCls, function(oneMCl){
#   cat(oneMCl,"\n")   
#   cvti <- cvt %>% filter(MCls==oneMCl)
#   Bi <- Bx[,cvti$bti]
#   X <- data.frame(x1=cvti$treats, x2=cvti$Batch2)
#   rn <- rownames(Bi)
#   ## Start loop By gene 
#   TMP <- mclapply(rn, function(ii){
#      y <- Bi[ii,]
#      nna <- countNA2(X, y)
#      dd <- myDE2(y, X, ii, nna)
#      dd
#   }, mc.cores=1) ### End loop by gene
#   ###  
#   TMP <- TMP[!is.na(TMP)]
#   TMP <- as.data.frame(do.call(rbind, TMP))%>%mutate(MCls=oneMCl)
#
#   TMP2 <- TMP%>%group_by(MCls, contrast)%>%
#           nest()%>%
#           mutate(qval=map(data, ~myqval((.x)$pval)))%>%
#           unnest(c(data,qval))%>%as.data.frame()   
#   TMP2
#})      
#opfn <- "./10_RNA.Variance_output/tmp7/4.1_mu"
#write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
#} ###


#######################
### Summary results ###
#######################


#####################
### (1). qq plots ###
#####################

## figfn <- "./10_RNA.Variance_output/tmp9/Figure4.1.qq.png"
## png(figfn, width=2000, height=2000, pointsize=12, res=300)
## par(mfrow=c(4,8),mar=c(4,4,1.5,2),mgp=c(2,1,0))
## x <- matrix(1:16, 4, 4, byrow=T)
## layout(x)

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
## res <- read.table(fn,header=T)%>%drop_na(pval)
## for (oneMCl in MCls){
## ##1  
##    res1 <- res %>%filter(MCls==oneMCl, contrast=="LPS") 
##    print(qq(res1$pval, main="LPS", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))

## ##2
##    res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")
##    print(qq(res2$pval, main="LPS-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
## ##3
##    res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")
##    print(qq(res3$pval, main="PHA", cex.main=0.8, cex.axis=0.8, cex.lab=0.8)) 
   
## ##4
##    res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")
##    print(qq(res4$pval, main="PHA-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
##    print(mtext(oneMCl, side=4, line=0.5, cex=0.8, col="blue") )
## }
## dev.off() 

## Sys.sleep(5)

## ############################
## ### (2), histogram plots ###
## ############################
## rm(list=ls())
## cat("hist plots for differential effects size", "\n")
## lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
## dx <- read.table(fn,header=T)%>%drop_na(beta) 
## fig0 <- ggplot(dx, aes(x=beta))+
##      geom_histogram(fill="grey70", color="grey20")+
##      xlab(bquote("Effective size"~"("~beta~")"))+
##      facet_grid(MCls~contrast, scales="free",labeller=labeller(contrast=lab2))+
##      theme_bw()+
##      theme(strip.background=element_blank())

## figfn <- "./10_RNA.Variance_output/tmp9/Figure4.2.hist.png"
## png(filename=figfn, width=800, height=800, pointsize=12, res=130)  
## print(fig0)
## dev.off()

## Sys.sleep(5)

rm(list=ls())

fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
dx <- read.table(fn, header=TRUE)%>%drop_na(beta)
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

### hist
fig2 <- ggplot(dx, aes(x=beta))+
   geom_histogram(aes(y=..density..), fill="white", colour="grey30")+
    xlab(bquote(~log[2]~"fold change"))+
    ylab("Density")+
    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
    theme_bw()+
    theme(strip.text=element_text(size=12))


figfn <- "./10_RNA.Variance_output/tmp9/Figure4.2_summryRes.pdf"
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig2, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()


############################
### (3), Barplots of DGM ###
############################

if(FALSE){

fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- unique(res$gene)
save(sigs, file="./10_RNA.Variance_output/tmp9/Sig4.DMG.RData")

}


if(FALSE){
###Barplots show No. DGE
cat("DMG Barplots", "\n")
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- res2%>%group_by(contrast, MCls)%>%summarise(ngene=n(),.groups="drop")

x <- res2%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
                    
fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=MCls))+
        geom_bar(stat="identity",position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(labels=lab2)+ylab("No. DMG")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.1_DMG.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}

### (3), barplots of DGM, up and down with light and deep colors, ***
if(FALSE){

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
###read data
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)

## up and down DGM
sigs <- res2%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
        

###Figure4.3.3, facet by contrast, and up and down together        
sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))
ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))

### 
fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col2comb, labels="")+ylab("DMG")+ylim(0,1600)+
        geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+50, fill=NULL), size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.3_DMG.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

### Figure4.3.4, facet by  MCls, up and down togther (stack)
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col1w <- colorspace::lighten(col1,0.3)
col1comb <- c(col1, col1w)
names(col1comb) <- paste(contrast, rep(c(1,2),each=4), sep="_") 

sig3 <- sigs%>%mutate(comb=paste(contrast, direction, sep="_"))
sig3$facet_fill_color <- col2[sig3$MCls] 
ann3 <- sig3%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
### 
fig0 <- ggplot(sig3, aes(x=contrast, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col1comb, labels="")+ylab("DMG")+ylim(0,1600)+
        geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+50, fill=NULL), size=3)+
        scale_x_discrete(labels=lab2)+
        facet_grid(~MCls)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.4_DMG.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

} ###

#############################################
### Final used for publication            ###
### Figure4.3.5, facet by contrast        ###
### and up above axis and down below axis ###
#############################################
rm(list=ls())
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


MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2, 0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
###read data
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)
x <- res%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))
x <- res%>%group_by(MCls)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

## up and down DGV
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2,-ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-800,800),5)
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
                          
###add star
anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
           mutate(pval=map_dbl(data, Mybinom), 
                  symb=map_chr(pval, Mysymb),
                  ypos=map_dbl(data, Mypos))%>%
           unnest(cols=c(contrast,MCls))   
                                  
###When provide new data frame into geom_text, use geom_bar(aes(fill))                           
fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2))+
   geom_bar(aes(fill=comb),stat="identity")+
   scale_fill_manual(values=col2comb, labels="")+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
                 vjust=ifelse(direction==2, 1.1, -0.2)), size=3)+ #
   scale_y_continuous(breaks=breaks_value, limits=c(-800,800),labels=abs(breaks_value))+
   ylab("Differentially expressed genes")+ 
   facet_grid(~contrast, labeller=facetlab)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
         strip.text=element_text(size=12),
         plot.margin=unit(c(5.5, 15, 5.5, 5.5), "points"))
              
fig0 <- fig0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)        

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.5_DMG.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
###

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.5_DMG.barplot5.pdf"
pdf(figfn, width=8, height=4)  
print(fig0)
grid.text("upregulated", x=unit(0.98,"npc"), y=unit(0.7,"npc"),
           rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
grid.text("downregulated", x=unit(0.98,"npc"), y=unit(0.38,"npc"),
           rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
dev.off()


### binomial test between up and down regulated genes
if(FALSE){
Mybinom <- function(subdf) {
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
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)

## up and down DGV
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
           mutate(pval=map_dbl(data,Mybinom))%>%
           unnest(cols=c(contrast,MCls))
           
}


##################################################################
### (4), correlation of differential effects across conditions ###
##################################################################


load("./10_RNA.Variance_output/tmp9/Sig4.DMG.RData")
Geneunq <- sigs

col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
          
######################
### (1) Using beta ###
######################  
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta2"
res <- read.table(fn,header=T)
TMP <- map_dfc(MCls, function(oneMCl){
   tmp <- map_dfc(Contrast, function(oneC){ 
      d0 <- res %>% filter(MCls==oneMCl, contrast==oneC)
      rownames(d0) <- d0$gene
      beta0 <- d0[Geneunq,"beta"]
      beta0
   })
   tmp 
})
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
y0 <- y[abs(y)<1.47] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL

###colors
x <- str_split(conditions, "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- conditions
tmp_colors <- list(celltype=col2, 
                   treatment=col1) #brewer.pal(4,"Set1")

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)
fig1 <- pheatmap(TMP0, col=mycol, breaks=mybreaks, 
         scale="none",
         border_color="NA",
         cluster_rows=T, cluster_cols=T, 
         annotation_col=tmp_column,
         annotation_colors=tmp_colors,
         show_colnames=T, show_rownames=F,
         na_col="white")

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.4.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "NKcell_LPS+DEX", "Tcell_LPS+DEX", 
#              "NKcell_PHA+DEX", "Tcell_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
#               "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#              "NKcell_PHA", "Tcell_PHA", "NKcell_LPS", "Tcell_LPS") 
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

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.4.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()
#corr <- cor(TMP0)
##mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)
#
#figfn <- "./10_RNA.Variance_output/tmp6/Figure4.4.2_corr.beta.png"
#png(figfn, width=1000, height=1000, res=180)
#print(corrplot(corr, method="color", order="hclust", hclust.method="complete", col=mycol,
#         tl.col="black", tl.cex=0.8, outline=F, diag=T))
#dev.off()



########################################################################################
### (5). scatterplots of beta(differential mean) between LPS/PHA and LPS-DEX/PHA-DEX ###
########################################################################################

if(FALSE){

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
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(file=fn,header=T)%>%mutate(rn2=paste(MCls, gene, sep="_"))
             
### (1), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH  
cat("(1)", "compare beta(differetial dispersion) between LPS and LPS-DEX", "\n")
dfa <- res%>%filter(contrast=="LPS")    
dfb <- res%>%filter(contrast=="LPS-DEX")%>%dplyr::select(rn2, beta, pval, qval)
       
df1 <- dfa%>%inner_join(dfb, by="rn2")

anno_df1 <- df1%>%group_by(MCls)%>%
           nest()%>%
           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
                  eq=map(corr,feq),
                  r2=map_dbl(corr, ~(.x)$estimate),
                  xpos=map_dbl(data,~xFun(.x,a=0.7)),
                  ypos=map_dbl(data,~yFun(.x,a=1)))%>%
           dplyr::select(-data,-corr)
     
fig1 <- ggplot(df1, aes(x=beta.x,y=beta.y))+
        geom_point(size=0.3, color="grey50")+ 
        geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
        facet_wrap(~MCls, nrow=2, scales="free")+
        scale_x_continuous("LPS effect size on mean", expand=expansion(mult=0.1))+
        scale_y_continuous("LPS+DEX effect size on mean", expand=expansion(mult=0.1))+
        theme_bw()+
        theme(strip.background=element_blank(),
              axis.title=element_text(size=10))
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- "./10_RNA.Variance_output/tmp9/Figure4.5.1_LPS.png"
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig1)
dev.off()

### (2), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
cat("(2)", "compare beta(differetial dispersion) between PHA and PHA-DEX", "\n")
dfa <- res%>%filter(contrast=="PHA")    
dfb <- res%>%filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, pval, qval)
       
df2 <- dfa%>%inner_join(dfb,by="rn2")

anno_df2 <- df2%>%group_by(MCls)%>%
           nest()%>%
           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
                  eq=map(corr,feq),
                  r2=map_dbl(corr, ~(.x)$estimate),
                  xpos=map_dbl(data,~xFun(.x,a=0.7)),
                  ypos=map_dbl(data,~yFun(.x,a=1)))%>%
           dplyr::select(-data,-corr)
     
fig2 <- ggplot(df2, aes(x=beta.x,y=beta.y))+
        geom_point(size=0.3, color="grey50")+ 
        geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
        facet_wrap(~MCls, nrow=2, scales="free")+
        scale_x_continuous("PHA effect size on mean", expand=expansion(mult=0.1))+
        scale_y_continuous("PHA+DEX effect size on mean", expand=expansion(mult=0.1))+
        theme_bw()+
        theme(strip.background=element_blank(),
              axis.title=element_text(size=10))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- "./10_RNA.Variance_output/tmp9/Figure4.5.2_PHA.png"
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig2)
dev.off()
} ###



##############################
### 6, compare DVG and DVG ###
##############################

###############################
### 6.1, Table show overlap ###
###############################
if(FALSE){

rm(list=ls())
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
x <- read.table("./10_RNA.Variance_output/tmp7/4_mu.meta", header=T)
gene2 <- unique(x$gene)

### (1). pseudo-bulk differential 
fn <- "./6_DEG.CelltypeNew_output/2_meta.rds"
df1 <- read_rds(fn)%>%
       mutate(zscore=beta/stderr)%>%
       drop_na(beta,qval)%>%filter(qval<0.1, abs(beta)>0.5)#filter(gene%in%gene2,qval<0.1,abs(beta)>0.5)
sig1 <- df1%>%
        group_by(contrast, MCls)%>%
        nest()%>%
        mutate(ngene=map_dbl(data,nrow),
               gene=map(data,~(.x)$gene),
               rn=paste(MCls,contrast,sep="_"))%>%
        dplyr::select(-data)


### (2). variance
fn <- "./10_RNA.Variance_output/tmp7/2.1_va2"
df2 <- read.table(file=fn,header=T)%>%
       mutate(zscore=beta/stderr)%>%
       drop_na(beta,qval)%>%filter(qval<0.1, abs(beta)>0.5)  
sig2<- df2%>%
        group_by(contrast, MCls)%>%
        nest()%>%
        mutate(ngene=map_dbl(data,nrow),
               gene=map(data,~(.x)$gene),
               rn=paste(MCls,contrast,sep="_"))%>%
        dplyr::select(-data)


### (3). dispersion
fn <- "./10_RNA.Variance_output/tmp7/3_phi.meta"
df3 <- read.table(fn,header=T)%>%
       mutate(zscore=beta/stderr)%>%
       drop_na(beta,qval)%>%filter(qval<0.1, abs(beta)>0.5)
sig3 <- df3%>%
        group_by(contrast, MCls)%>%
        nest()%>%
        mutate(ngene=map_dbl(data,nrow),
               gene=map(data,~(.x)$gene),
               rn=paste(MCls,contrast,sep="_"))%>%
        dplyr::select(-data)


### (4). residual dispersion
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew3.meta"
df4 <- read.table(fn,header=T)%>%mutate(zscore=beta/stderr)%>%drop_na(beta,qval)%>%filter(qval<0.1, abs(beta)>0.5)
sig4 <- df4%>%
        group_by(contrast, MCls)%>%
        nest()%>%
        mutate(ngene=map_dbl(data,nrow),
               gene=map(data,~(.x)$gene),
               rn=paste(MCls,contrast,sep="_"))%>%
        dplyr::select(-data)

### (5). mean expression
fn <- "./10_RNA.Variance_output/tmp7/4_mu.meta"
df5 <- read.table(fn, header=T)%>%
       mutate(zscore=beta/stderr)%>%
       drop_na(beta,qval)%>%
       filter(qval<0.1,abs(beta)>0.5)

sig5 <- df5%>%
        group_by(contrast, MCls)%>%
        nest()%>%
        mutate(ngene=map_dbl(data,nrow),
               gene=map(data,~(.x)$gene),
               rn=paste(MCls,contrast,sep="_"))%>%
        dplyr::select(-data)


##overlap functions
Overlap <- function(dx, dy){

   flen <- function(x,y) length(intersect(x,y))
   flen2 <- function(x,y) length(union(x,y))
            
   res0 <- dx%>%left_join(dy, by="rn")%>%
        mutate(n12=map2_dbl(gene.x, gene.y, flen),nt=map2_dbl(gene.x, gene.y, flen2))%>%
        dplyr::select(-gene.x, -gene.y, -contrast.y, -MCls.y)%>%
        dplyr::rename(contrast=contrast.x, MCls=MCls.x, nx=ngene.x, ny=ngene.y)%>%
        mutate(Jindex=n12/nt)
   res0
}

x1 <- Overlap(sig1, sig2) ## pseudo-bulk vs variance
x2 <- Overlap(sig1, sig3) ## pseudo-bulk vs dispersion
x3 <- Overlap(sig1, sig4) ## pseudo-bulk vs corrected dispersion
x4 <- Overlap(sig1, sig5) ## pseudo-bulk vs mean      
##
} ###6.1, End 


#########################################################
#### 6.2, scatter plot of beta between va, phi and mu ###
#########################################################

rm(list=ls())
### my fun 1, generate data frame used for plots 
myDFxy <- function(dfx, dfy){
###
   dfx <- dfx%>%dplyr::select(zscore, beta, qval, rn, contrast, MCls, gene)
   dfy <- dfy%>%dplyr::select(zscore, beta, qval, rn)       
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

### label expression function       
feq <- function(x){
  #r <- format(as.numeric(x$estimate),digits=1)
  r <- round(as.numeric(x$estimate), digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  eq
  #r 
}

#################
### Read data ###
#################

### df1, pseudo-bulk differential 
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
df1 <- read_rds(fn)%>%drop_na(beta, qval)%>%mutate(zscore=beta/stderr)

### df2, residual dispersion
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(file=fn,header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 
       
### df3, mean expression
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df3 <- read.table(file=fn, header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 

###
mycol <- c("1"="grey", "2"="#118C4F", "3"="blue", "4"="#9400D3")
## mycol <- c("1"="grey20", "2"="red", "3"="blue", "4"="#9400D3")
## mycol <- c("1"="#bababa", "2"="#de2d26", "3"="#6baed6", "4"="#756bb1")
Newcon2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")


#### (1). mean vs residual dispersion
cat("(1).", "mean vs residual dispersion", "\n")
dfxy1 <- myDFxy(df3, df2)
mylabel <- c("1"="NS", "2"="DEG(only)", "3"="DVG(only)", "4"="Both")

anno_df1 <- dfxy1%>%
            group_by(contrast, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
                   eq=map(corr,feq))%>%
           dplyr::select(-data,-corr)
        
fig1 <- ggplot(dfxy1, aes(x=zscore.x, y=zscore.y))+
        rasterise(geom_point(aes(colour=factor(gr)), size=0.3),dpi=300)+
        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
        geom_hline(yintercept=0, color="grey60")+
        geom_vline(xintercept=0, color="grey60")+
        geom_text(data=anno_df1, x=0, y=22, aes(label=eq), size=3, parse=T)+
        xlab("Z score of gene mean expression")+ylab("Z score of gene variability")+
        facet_grid(MCls~contrast, labeller=labeller(contrast=Newcon2))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.text=element_text(size=9),
              #legend.key.size=unit(0.4,units="cm"),
              axis.text=element_text(size=9))
###              
## figfn <- "./10_RNA.Variance_output/tmp9/Figure5.1_DMGvsDVG.png"
## png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
## print(fig1)
## dev.off()

figfn <- "./10_RNA.Variance_output/tmp9/Figure5.1_DMGvsDVG.pdf"
pdf(figfn, width=8, height=7)  
print(fig1)
dev.off()


###sumamary number
x <- dfxy1%>%
     filter(gr==4)%>%
     group_by(MCls, contrast)%>%nest()%>%
     mutate(ngene=map_dbl(data,nrow))


### (2). mean and gene expression
dfxy2 <- myDFxy(df3, df1)
mylabel <- c("1"="NS", "2"="DEG(only, mean)", "3"="DEG(only, pseudo-bulk)", "4"="Both")
anno_df2 <- dfxy2%>%
   group_by(contrast, MCls)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
          eq=map(corr,feq))%>%
   dplyr::select(-data,-corr)
        
fig2 <- ggplot(dfxy2, aes(x=zscore.x, y=zscore.y))+
   rasterise(geom_point(aes(colour=factor(gr)), size=0.3))+
   scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
   geom_text(data=anno_df2, x=0, y=50, aes(label=eq), size=3, parse=T)+
   geom_hline(yintercept=0, color="grey60")+
   geom_vline(xintercept=0, color="grey60")+
   xlab("Z score of gene mean")+ylab("Z score of gene expression")+
   facet_grid(MCls~contrast, labeller=labeller(contrast=Newcon2))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         #legend.key.size=unit(0.4,units="cm"),
         axis.text=element_text(size=9))

###              
figfn <- "./10_RNA.Variance_output/tmp9/Figure5.2_DMGvsDEG.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig2)
dev.off()

###
figfn <- "./10_RNA.Variance_output/tmp9/Figure5.2_DMGvsDEG.pdf"
pdf(figfn, width=8, height=7)  
print(fig2)
dev.off()


##################################
### barplots, show DMG and DVG ###
##################################

##rm(list=ls())

### df2, residual dispersion
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(file=fn,header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 
       
### df3, mean expression
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df3 <- read.table(file=fn, header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 

dfcomb <- myDFxy(df3, df2)
#mylabel <- c("1"="NS", "2"="DMG(only)", "3"="DVG(only)", "4"="Both")

dfcomb <- dfcomb%>%filter(gr!=1)
dfcomb$gr[dfcomb$gr==2] <- 1
dfcomb$gr[dfcomb$gr==4] <- 2
##1="DEG", 2="Both", 3="DVG"  

sigs <- dfcomb%>%group_by(MCls, contrast, gr)%>%summarise(ngene=n(),.groups="drop")
sigs <- sigs%>%mutate(comb=paste(MCls, gr, sep="_"))%>%as.data.frame()

anno2 <- dfcomb%>%
    group_by(MCls, contrast)%>%summarise(ngene=n(),.groups="drop")%>%as.data.frame()

###
sigs%>%filter(MCls=="Monocyte")%>%
    group_by(MCls,contrast)%>%
    mutate(ntotal=sum(ngene),prop=ngene/ntotal)%>%filter(gr==2)


lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

## colw <- lapply(MCls, function(ii){
##          x1 <- colorspace::lighten(col2[ii], 0)
##          x2 <- colorspace::lighten(col2[ii], 0.6)
##          x3 <- colorspace::lighten(col2[ii], 0.3)
##          xx <- c(x1, x2, x3)
##          names(xx) <- paste(ii, 1:3, sep="_")
##          xx
##          })
## colw <- unlist(colw)                                             


### (1)
fig0 <- ggplot(sigs)+
   geom_bar(stat="identity", position=position_stack(reverse=T),
            aes(x=MCls, y=ngene, fill=MCls, alpha=factor(gr)))+
   scale_fill_manual(values=col2)+
   scale_alpha_manual(values=c("1"=1, "2"=0.3, "3"=0.6))+
   geom_text(data=anno2, aes(x=MCls, label=ngene, y=ngene+50, fill=NULL), size=3)+
   scale_y_continuous(breaks=seq(0,1500,300), limits=c(0,1600))+ 
   facet_grid(~contrast, labeller=facetlab)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
legend2 <- get_legend(
   ggplot(sigs%>%filter(MCls=="Monocyte"),aes(x=contrast,y=ngene))+
   geom_bar(stat="identity", position=position_stack(reverse=T),
      fill=col2[["Monocyte"]], aes(alpha=factor(gr)))+
   #aes(fill=comb))+
   scale_alpha_manual(values=c("1"=1, "2"=0.3, "3"=0.6), 
      labels=c("1"="DEG", "2"="Both", "3"="DVG"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.text=element_text(size=8),
         strip.text=element_text(size=12),
         legend.key.size=grid::unit(1,"lines")))

###
figfn <- "./10_RNA.Variance_output/tmp9/Figure6.1_bar.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)
print(plot_grid(fig0, legend2, rel_widths=c(4,0.5)))
dev.off()

##pdf
figfn <- "./10_RNA.Variance_output/tmp9/Figure6.1_bar.pdf"
pdf(figfn, width=8, height=4)
print(plot_grid(fig0, legend2, rel_widths=c(4,0.5)))
dev.off()


###
### (2), equal height
fig0 <- ggplot(sigs)+
   geom_bar(stat="identity", position=position_fill(reverse=T),
            aes(x=MCls, y=ngene, fill=MCls, alpha=factor(gr)))+
   scale_fill_manual(values=col2)+
   scale_alpha_manual(values=c("1"=1, "2"=0.3, "3"=0.6))+
   geom_text(data=anno2, aes(x=MCls, y=1.2, label=ngene), vjust=-0.2, position="fill", size=3)+
   scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+ 
   facet_grid(~contrast, labeller=facetlab)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
legend2 <- get_legend(
   ggplot(sigs%>%filter(MCls=="Monocyte"),aes(x=contrast,y=ngene))+
   geom_bar(stat="identity", position=position_fill(reverse=T),
      fill=col2[["Monocyte"]], aes(alpha=factor(gr)))+
   #aes(fill=comb))+
   scale_alpha_manual(values=c("1"=1, "2"=0.3, "3"=0.6), 
      labels=c("1"="DEG", "2"="Both", "3"="DVG"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(1,"lines")))

###
figfn <- "./10_RNA.Variance_output/tmp9/Figure6.2_bar.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)
print(plot_grid(fig0, legend2, rel_widths=c(4,0.5)))
dev.off()  
      
