#module load R 
#setwd("/nfs/rprdata/julong/SCAIP/analyses/TWAS/UKBB2.summary")
##
#library(org.Hs.eg.db)
library(annotables) 
library(data.table)
library(tidyverse) 
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra) 
library(ggExtra)
library(RColorBrewer)
theme_set(theme_grey())
##

### blood
Child <- c("LINGO4", "PSMD4", "S100A12", "TDRKH", "CCL20", "TLR10",
           "MLEC", "SPPL3", "GSDMA", "GSDMB", "IKZF3", "LINC00672",
           "MED24", "MSL1", "ORMDL3", "STARD3", "STK19P")
Share <- c("AFF4", "IL4", "PDLIM4", "SLC22A5", "AC007248.6", "IL18R1", "IL18RAP",
           "D2HGDH", "SLC25A46","AGER", "ATF6B", "BAG6", "C4A", "C4B", "C6orf48",
           "CLIC1", "CYP21A2", "DHFRP2", "FKBPL", "FLOT1",
           "HLA-C", "HLA-DOB", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
           "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DRB9", "HSPA1L", "LY6G5B",
           "MICA", "MICB", "NOTCH4", "PSMB9", "TAP2", "TCF19", "ZBTB12", "OR2H2", "RP11-672A2.4", "RP11-619A14.3", "RPS26", "SUOX", "SMAD3", "IL4R")
res <- list(Child=Child, Share=Share)
opfn <- "./results/Old_blood.rds"
write_rds(res, opfn)

### Lung
Adult <- c("VPS52")
Child <- c("FLG", "LINGO4", "PSMD4", "S100A12", "TDRKH", "FAM114A1", "STEAP1B", "AP5B1", "SPPL3", "CSF3", "GSDMA", "GSDMB", "LINC00672", "MED24", "MSL1", "ORMDL3",
        "PGAP3", "PNMT", "SMARCE1", "HLA-DPA1", "STK19P", "RP11-672A2.4")
Share <- c("IL4", "KIF3A", "RAD50", "SLC22A5", "IL18R1", "IL18RAP", "AC093642.3", "D2HGDH", "GAL3ST2", "SLC25A46", "TMEM232", "ABHD16A", "AGER", "APOM", "ATF6B", "BAG6", "C4A", "C4B", "C6orf48", "CYP21A2", "FLOT1", "GPANK1",
  "HLA-C", "HLA-DOB", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
  "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DRB9", "LY6G5B", "MDC1", 
   "MICA", "MICB", "NOTCH4", "PPP1R18", "PSMB9", "TAP2", "TCF19", "TNXA",
   "RPS26", "SUOX", "SMAD3", "IL4R", "CLEC16A")
res <- list(Adult=Adult, Child=Child, Share=Share)
opfn <- "./results/Old_lung.rds"
write_rds(res, opfn)




###
### whole blood

################
### qq plots ###
################  

res <- fread("./results/Blood.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
                aes(x = 0, xend = expected, y = 0, yend = expected),
                size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
   labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
        y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
   theme_bw()
  
figfn <- "./figures/Figure1.1_blood.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


#######################
### Manhattan plots ###
#######################

res <- fread("./results/Blood.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%
       distinct(ensgene,chr,.keep_all=T)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure1.2_blood.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


###
###
res <- fread("./results/Blood.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       distinct(ensgene, chr, .keep_all=T)%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%arrange(chr)

write.csv(x, file="./results/Blood.sign.gene.csv", row.names=F)

x1 <- read.csv("./results/Blood.sign.gene.csv")
x2 <- read_rds("./results/Old_blood.rds")

############
### skin ###
############

res <- fread("./results/Skin_Not_Sun.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
  theme_bw()
  
figfn <- "./figures/Figure2.1_skin_Notsun.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


#######################
### Manhattan plots ###
#######################

res <- fread("./results/Skin_Not_Sun.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure2.2._skin_Notsun.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


###
### 
res <- fread("./results/Skin_Not_Sun.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%
     arrange(chr)%>%
     distinct(gene,gene_name,chr,.keep_all=TRUE)

write.csv(x, file="./results/Skin_Not_Sun.sign.gene.csv", row.names=F)



########################
### skin sun exposed ###
########################

res <- fread("./results/Skin_Sun.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
  theme_bw()
  
figfn <- "./figures/Figure3.1_skin_sun.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


#######################
### Manhattan plots ###
#######################

res <- fread("./results/Skin_Sun.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure3.2._skin_sun.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


###
### 
res <- fread("./results/Skin_Sun.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%
     arrange(chr)%>%
     distinct(gene,gene_name,chr,.keep_all=TRUE)

write.csv(x, file="./results/Skin_Sun.sign.gene.csv", row.names=F)




############
### Lung ###
############

### qq plots ###

res <- fread("./results/Lung.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
  theme_bw()
  
figfn <- "./figures/Figure4.1_lung.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


### Manhattan plots ###

res <- fread("./results/Lung.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure4.2._lung.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


### 
res <- fread("./results/Lung.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%arrange(chr)%>%
     distinct(gene, gene_name, chr, .keep_all=TRUE)

write.csv(x, file="./results/Lung.sign.gene.csv", row.names=F)

### overlap 
x1 <- read.csv("./results/Lung.sign.gene.csv")
gene1 <- unique(x1$gene_name) 
x2 <- read_rds("./results/Old_lung.rds")

ii <- unique(intersect(gene1, x2$Share))

#######################
### small intestine ###
#######################


### qq plots ###

res <- fread("./results/Small_intestine.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
  theme_bw()
  
figfn <- "./figures/Figure5.1_small_intestine.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


### Manhattan plots ###

res <- fread("./results/Small_intestine.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure5.2_small_intestine.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


### 
res <- fread("./results/Small_intestine.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%arrange(chr)%>%
     distinct(gene, gene_name, chr, .keep_all=TRUE)

write.csv(x, file="./results/Small_intestine.sign.gene.csv", row.names=F)



#######################
### Spleen ###
#######################


### qq plots ###

res <- fread("./results/Spleen.csv", header=T)
ngene <- nrow(res)               
plotdata <- data.frame(observed=-log10(sort(res$pvalue)),
                       expected=-log10(ppoints(ngene)),
                       clower=-log10(qbeta(p=(1-0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))),
                       cupper=-log10(qbeta(p=(1+0.95)/2,shape1=seq(ngene),shape2=rev(seq(ngene)))))
                       
p <- ggplot(plotdata, aes(x=expected,y=observed))+
   geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey30", alpha=0.5) +
   geom_step(color="red", size = 1.1, direction = "vh") +
   geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
  theme_bw()
  
figfn <- "./figures/Figure6.1_spleen.qq.png"
png(figfn, width=400, height=400, pointsize=12, res=120)
print(p)
dev.off()


### Manhattan plots ###

res <- fread("./results/Spleen.csv", header=T)

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       drop_na(chr)%>%as.data.frame()
       
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:max(res$chr)){
   npi <- max(res[res$chr==i,"start"])
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- s+npi
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)
   
### plot   
sig <- -log10(0.05/ngene)
p2 <- ggplot(res, aes(x=BPcum, y=-log10(pvalue), 
                      color=factor(chr), size=-log10(pvalue)))+
      geom_point(alpha=0.75)+
      geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
      scale_x_continuous("", label=axis.set$chr, breaks=axis.set$midpoint,limits=c(min(res$BPcum),max(res$BPcum)))+
      scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"))+
      scale_color_manual(values=rep(c("grey0", "grey70"), nchr))+
      scale_size_continuous(range=c(0.2,2)) +
      theme_bw()+
      theme(legend.position="none")

###
figfn <- "./figures/Figure6.2_spleen.manhattan.png"
png(figfn, width=800, height=400, pointsize=12, res=120)
print(p2)
dev.off()


### 
res <- fread("./results/Spleen.csv", header=T)

ngene <- nrow(res)
sigs <- 0.05/ngene

anno <- grch38%>%dplyr::select(ensgene, chr, start, end)
res <- res%>%
       mutate(ensgene=gsub("\\..*", "", gene))%>%
       left_join(anno, by="ensgene")%>%
       mutate(chr=as.numeric(chr))%>%
       arrange(desc(zscore))%>%as.data.frame()
       
x <- res%>%filter(pvalue<sigs)%>%arrange(chr)%>%
     distinct(gene, gene_name, chr, .keep_all=TRUE)

write.csv(x, file="./results/Spleen.sign.gene.csv", row.names=F)



