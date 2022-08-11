###
library(tidyverse)
library(parallel)
library(purrr)
library(reshape)
library(qqman)
library(qvalue)
##
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
library(RColorBrewer)
library(viridis)

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


#################################################################
### 3 (2) Differential dispersion after removing mean effects ###
#################################################################


####################################
### Differential procedure ###
####################################

## Differential residual dispertion
load("./10_RNA.Variance_output/tmp9/1.2_Sel.PhxNew.RData")
rn <- rownames(PhxNew2)
rownames(PhxNew2) <- gsub("\\.[0-9]*", "", rn)

###
bti2 <- colnames(PhxNew2)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
comb <- unique(cvt$comb)


### 4.1, estimate differetial results by batch
cat("Differential analysis by batch, gene,", nrow(PhxNew2), "\n")
res <- map_dfr(comb, function(oneX){
   cat(oneX,"\n")
   cvti <- cvt %>%filter(comb==oneX)
   oneComb <- unlist(strsplit(oneX, "_"))
   oneMCl <- oneComb[1]
   oneBatch <- oneComb[2]
   Phi <- PhxNew2[,cvti$bti]
   X <- data.frame(x1=cvti$treats)
   rn <- rownames(Phi)
   ## Start loop By gene 
   TMP <- mclapply(rn, function(ii){
      y <- Phi[ii,]
      nna <- countNA(X$x1,y)
      dd <- myDE(y, X, ii, nna, threshold=3)
      dd
   }, mc.cores=1) ### End loop by gene
   ###  
   TMP <- TMP[!is.na(TMP)]
   TMP <- as.data.frame(do.call(rbind, TMP))%>%mutate(MCls=oneMCl, batch=oneBatch)
   TMP       
})

opfn <- "./10_RNA.Variance_output/tmp9/3_phiNew.results"
write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")



### meta analysis
### "Meta analysis"
### Read data
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.results"
res <- read.table(fn,header=T)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
       filter(batch%in%c("SCAIP1","SCAIP4", "SCAIP5", "SCAIP6"))

### Meta
res2 <- res%>%group_by(rn)%>%
        nest()%>%
        mutate(outlist=mclapply(data, myMeta, mc.cores=1))%>%
        dplyr::select(-data)%>%unnest(c(outlist))%>%as.data.frame()
cvt <- str_split(res2$rn, "_", simplify=T)
res2 <- res2%>%mutate(MCls=cvt[,1], contrast=cvt[,2], gene=cvt[,3])
          
### add qvalue
res3 <- res2%>%group_by(MCls, contrast)%>%
        nest()%>%
        mutate(qval=map(data, ~myqval((.x)$pval)))%>%
        unnest(c(data,qval))%>%as.data.frame()
        
#opfn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
opfn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta2" ##remove batch 2 and 3
write.table(res3, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")

} ###End, Differential analysis


#######################
### Summary results ###
#######################

#####################
### (1). qq plots ###
#####################

## rm(list=ls())
## cat("Show qq plots", "\n")
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

## figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.1.qq.png"
## png(figfn, width=2000, height=2000, pointsize=12, res=300)
## par(mfrow=c(4,8),mar=c(4,4,1.5,2),mgp=c(2,1,0))
## x <- matrix(1:16, 4, 4, byrow=T)
## layout(x)

## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T)%>%drop_na(pval) 
## ###
## for (oneMCl in MCls){   
## ##1  
##    res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS") 
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
## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## dx <- read.table(fn,header=T)%>%drop_na(beta) 
## lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
##           "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## fig0 <- ggplot(dx, aes(x=beta))+
##      geom_histogram(fill="grey70", color="grey20")+
##      xlab(bquote("Effective size"~"("~beta~")"))+
##      facet_grid(MCls~contrast,scales="free",labeller=labeller(contrast=lab2))+
##      theme_bw()+
##      theme(strip.background=element_blank())

## figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.2.hist.png"
## png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
## print(fig0)
## dev.off()

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=TRUE)%>%drop_na(beta)

res <- res[,c(1,2,7,4:6,8)]
write.table(res, "./10_RNA.Variance_output/Differential_variability.txt",
   quote=F, sep="\t", row.names=F)

system("tar -czvf ./10_RNA.Variance_output/Differential_variability.txt.gz ./10_RNA.Variance_output/Differential_variability.txt")

## Sys.sleep(5)
rm(list=ls())

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
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
   ggrastr::rasterise(geom_point(size=0.3, colour="grey30"), dpi=300)+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12)) 

### hist
fig2 <- ggplot(dx, aes(x=beta))+
   geom_histogram(aes(y=..density..), fill="white", colour="grey30")+
   geom_vline(xintercept=c(-0.5,0.5), color="red", linetype="dashed", alpha=0.5, size=0.5)+ 
    xlab(bquote(~log[2]~"fold change"))+
    ylab("Density")+
    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
    theme_bw()+
    theme(strip.text=element_text(size=12))


figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.2_summryRes_reviews.pdf"
pdf(figfn, width=12, height=6)
print(plot_grid(fig1, fig2, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()



############################
### (3), Barplots of DGP ###
############################


fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- unique(res$gene)
save(sigs, file="./10_RNA.Variance_output/tmp9/Sig3x.DGP.RData") 


#### Barplots show NO.DGV together(Up and down)
## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(file=fn,header=T)
## res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
## sigs <- res2%>%group_by(contrast, MCls)%>%summarise(ngene=n(), .groups="drop")

## x <- res2%>%group_by(contrast)%>%summarise(ngene=n_distinct(gene), .groups="drop")
## x2 <- res2%>%group_by(MCls)%>%summarise(ngene=n_distinct(gene), .groups="drop")

## cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##           "NKcell"="#aa4b56", "Tcell"="#ffaa00")
## lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
##           "PHA"="PHA", "PHA-DEX"="PHA+DEX")          
## fig0 <- ggplot(sigs,aes(x=contrast,y=ngene,fill=MCls))+
##         geom_bar(stat="identity",position=position_dodge())+
##         scale_fill_manual(values=cols)+
##         scale_x_discrete(labels=lab2)+ylab("No. DVG")+
##         theme_bw()+
##         theme(legend.title=element_blank(),
##               axis.title.x=element_blank())
## ###
## figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.1_DGP.barplot.png"
## png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
## print(fig0)
## dev.off()


## ### (3), barplots of DGP, up and down with light and deep colors, ***
## if(FALSE){

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

## ### colors
## col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##           "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
## col2w <- colorspace::lighten(col2,0.3)
## col2comb <- c(col2, col2w)
## names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
## ###read data
## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)

## ## up and down DGV
## sigs <- res%>%
##         mutate(direction=ifelse(beta>0, "1", "2"))%>%
##         group_by(contrast, MCls, direction)%>%
##         summarise(ngene=n(),.groups="drop")
        
        
## ### Figure3x.3.3, facet by contrast and up and down together(stack)       
## sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))
## ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
## #xpos <- c("Bcell"=0.5,"Monocyte"=1, "NKcell"=1.5, "Tcell"=2)
## #ann2$xpos <- xpos[ann2$MCls]
## facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
##                           "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
 
## fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
##         geom_bar(stat="identity", position="stack")+ 
##         scale_fill_manual(values=col2comb, labels="")+ylab("DVG")+ylim(0,600)+
##         geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+20, fill=NULL), size=3)+
##         facet_grid(~contrast, labeller=facetlab)+
##         theme_bw()+
##         theme(legend.position="none",
##               axis.title.x=element_blank(),
##               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

## figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.3_DGP.barplot3.png"
## png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
## print(fig0)
## dev.off()


## ### Figure3x.3.4, facet by cell type, up and down together(stack)
## contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
## col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## col1w <- colorspace::lighten(col1,0.3)
## col1comb <- c(col1, col1w)
## names(col1comb) <- paste(contrast, rep(c(1,2),each=4), sep="_") 

## sig3 <- sigs%>%mutate(comb=paste(contrast, direction, sep="_"))
## sig3$facet_fill_color <- col2[sig3$MCls]
## lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
##           "PHA"="PHA", "PHA-DEX"="PHA+DEX")
                
## ann3 <- sig3%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
## ### 
## fig0 <- ggplot(sig3, aes(x=contrast, y=ngene, fill=comb))+
##         geom_bar(stat="identity", position="stack")+ 
##         scale_fill_manual(values=col1comb, labels="")+ylab("DVG")+ylim(0,600)+
##         geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+20, fill=NULL), size=3)+
##         scale_x_discrete(labels=lab2)+
##         facet_grid(~MCls)+
##         theme_bw()+
##         theme(legend.position="none",
##               axis.title.x=element_blank(),
##               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

## figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.4_DGP.barplot4.png"
## png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
## print(fig0)
## dev.off()

#################################################
### Final Barplot used for paper Figure3x.3.5 ###
### barplots of DVG, facet by contrast        ###
### and up above axis and down below axis *** ###
#################################################
### used for paper   
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
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
###read data
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
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
breaks_value <- pretty(c(-400,300),5)
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
                  vjust=ifelse(direction==2, 1.1, -0.2)), size=3)+ #
        scale_y_continuous("", breaks=breaks_value, limits=c(-450,300),labels=abs(breaks_value))+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
## fig0 <- fig0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3) 

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.5_DGP.barplot5_reviews.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()



### binomial test between up and down regulated genes
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
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)

## up and down DGV
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
           mutate(pval=map_dbl(data,Mybinom))%>%
           unnest(cols=c(contrast,MCls))
           


##################################################################
### (4), correlation of differential effects across conditions ###
##################################################################

load("./10_RNA.Variance_output/tmp9/Sig3x.DGP.RData")
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
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta2"
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
rownames(TMP) <- Geneunq
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
y0 <- y[abs(y)<6.18] #99% percent quantile(abs(y),probs=0.99)
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
         annotation_colors=tmp_colors,
         show_colnames=T, show_rownames=F,
         na_col="white")

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.4.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
#Neworder <- c("Bcell_LPS+DEX", "Bcell_LPS", "Bcell_PHA+DEX",  "Bcell_PHA",
#              "Monocyte_LPS+DEX", "Monocyte_LPS", "Monocyte_PHA+DEX", "Monocyte_PHA", 
#              "Tcell_LPS+DEX","Tcell_LPS", "Tcell_PHA+DEX", "Tcell_PHA", 
#               "NKcell_LPS+DEX", "NKcell_LPS", "NKcell_PHA+DEX", "NKcell_PHA") 
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "NKcell_LPS+DEX", "Tcell_LPS+DEX",
#              "NKcell_PHA+DEX", "Tcell_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",               
#              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#              "NKcell_PHA", "Tcell_PHA","NKcell_LPS", "Tcell_LPS") 
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

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.4.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()


#corr <- cor(TMP0)
##mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)
#
#figfn <- "./10_RNA.Variance_output/tmp6/Figure3x.4.2_corr.beta.png"
#png(figfn, width=1000, height=1000, res=180)
#print(corrplot(corr, method="color", order="hclust", hclust.method="complete", col=mycol,
#         tl.col="black", tl.cex=0.8, outline=F, diag=T))
#dev.off()



############################################################
### (5), compare beta (differential residual dispersion) ###
###         from LPS/PHA and LPS-DEX/PHA-DEX             ###
############################################################
#
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
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(file=fn,header=T)%>%mutate(rn2=paste(MCls, gene, sep="_"))
             
### (1), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH  
cat("(1)", "compare beta(differential residual dispersion) between LPS and LPS-DEX", "\n")
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
     
fig1 <- ggplot(df1, aes(x=beta.x,y=beta.y))+
   geom_point(size=0.3,color="grey50")+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
   facet_wrap(~MCls, nrow=2, scales="free")+
   scale_x_continuous("LPS effect on gene variability", expand=expansion(mult=0.1))+
   scale_y_continuous("LPS+DEX effect on gene variability", expand=expansion(mult=0.1))+
   theme_bw()+
   theme(strip.background=element_blank(),
         axis.title=element_text(size=12))
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.5.1_LPS.png"
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
   scale_x_continuous("PHA effect on gene variability", expand=expansion(mult=0.1))+
   scale_y_continuous("PHA+DEX on gene variability", expand=expansion(mult=0.1))+
   theme_bw()+
   theme(strip.background=element_blank(),
         axis.title=element_text(size=12))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)        
                           
figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.5.2_PHA.png"
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig2)
dev.off()


###
###
figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.5.3_comb.pdf"
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig2, nrow=1, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()
    
    
