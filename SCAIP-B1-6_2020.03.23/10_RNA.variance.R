#
rm(list=ls())
source("./Bin/LibraryPackage.R")

outdir <- "./10_RNA.Variance_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)
outdir2 <- "./10_RNA.Variance_output/tmp9/"
if(!file.exists(outdir2)) dir.create(outdir2, showWarnings=F)
source("./Bin/Funs.R")

#############################################################
###    gene expression variance and Differetial analysis  ###
###           Last modified by Julong wei, 2020-12-31     ###
#############################################################


########################################################
### 1. Calculate gene variance, dispersion and mean  ###
###       based on negative binomial distribution    ###
########################################################


######################
### (1). Read data ###
######################

### 1.1, read data
if(FALSE){
cat("1.","Read data from seurat object", "\n")
cat("Generate counts data(X) spliced genes in the autochromosome", "\n")
cat("Generate combinations(bti)", "\n")
 
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
count <- sc@assays$RNA@counts         
meta <- sc@meta.data

grch38_unq <- grch38%>%
              distinct(ensgene,.keep_all=T)%>%
              dplyr::select(ensgene, symbol, chr, biotype)
vars <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn),
               ensgene2=gsub("[SU]-", "", rn), 
               uns=grepl("S-",rn), rnz=rowSums(count))%>%
        left_join(grch38_unq, by="ensgene")

autosome <- as.character(1:22)        
varsSel <- vars%>%filter(uns, rnz>20, chr%in%autosome)#, grepl("protein_coding", biotype))
### Out of 17,697 protein coding gene, 15,770 with reads>20
### Out of 42,554 genes, 30,972 whole genes>20       
X <- count[varsSel$rn,]
rownames(X) <- varsSel$ensgene2

### obs, cells information
obs <- meta%>%#filter(MCls=="Monocyte")%>% 
       mutate(bti=paste(MCls, treats, BEST.GUESS, BATCH, sep="_"))%>%
       dplyr::select(NEW_BARCODE, bti, MCls)

dd <- obs%>%group_by(bti)%>%summarise(ncell=n())
#save(dd, file="./10_RNA.Variance_output/tmp9/0_ncell.RData")
bti <- dd%>%filter(ncell>20)%>%dplyr::select(bti)%>%unlist() ##15 at 5%; 25 at 10%, default 20

} ### End, 1

### 1.2, Distribution of number of cells per combination(MCl+individual+treatment, 1536) 
### 5%, 15
if(FALSE){
cat("0.", "Distribution of number of cells per combination", "\n")
load("./10_RNA.Variance_output/tmp9/0_ncell.RData")
cvt <- str_split(dd$bti, "_", simplify=T)
dd0 <- data.frame(MCls=cvt[,1], treats=cvt[,2], ncell=dd$ncell)
p <- ggplot(dd0, aes(x=ncell))+
     geom_histogram(fill="grey70", color="grey40")+xlab("No.cell")+
     facet_wrap(~MCls,nrow=2,scales="free")+
     theme_bw()+
     theme(plot.title=element_text(hjust=0.5))

png("./10_RNA.Variance_output/tmp9/Figure0.0.ncell.png", height=600, width=700, res=120)
print(p)
dev.off()
}


##############################################
### (2), Estimate gene expression variance ###
##############################################

### 2.1, Estimate gene expression variance based NB model 
if(FALSE){  
cat("2.1", "Estimate gene expression variance based NB model", "\n")
## calcualte for each combinations, from the same individual, the same cell type and treats
size_after <-  median(colSums(X))  #tmp9
#size_after <- 1e+06                #tmp10
TMP <- lapply(1:length(bti), function(i){   
### 
   time0 <- Sys.time()
   bti0 <- as.character(bti[i])
   celli <- as.character(obs[obs$bti==bti[i], "NEW_BARCODE"])
   
   gene <- rownames(X)
   v <- rep(NA,nrow(X))
   names(v) <- gene
   b <- rep(NA,nrow(X))
   names(b) <- gene
   phi <- rep(NA,nrow(X))
   names(phi) <- gene
   ##
   se.mu <- rep(NA,nrow(X))
   names(se.mu) <- gene
   se.phi <- rep(NA,nrow(X))
   names(se.phi) <- gene
###
   Xi <- X[,celli]
   rnz <- rowSums(Xi)
   size <- colSums(Xi)/size_after  ###0.001        
   
   ncell <- ncol(Xi)
   
   Xe <- Xi>0
   nnz <- rowSums(Xe) ## number of cells with non-zero reads 
   #Xe <- Xi==1
   #nn1 <- rowSums(Xe)      
   gene0 <- gene[(rnz>15)&(nnz>15)]

   ## for each gene
   tmp <- mclapply(gene0,function(k){
      xk <- as.numeric(Xi[k,])
      ## given initial value by mean, var 
      mu0 <- mean(xk/size)
      va0 <- var(xk/size)
      phi0 <- va0/mu0^2
      theta <- c(log(mu0),log(1/phi0))
       
      parm <- try(optim(par=theta, fn=nb_llik, x=xk, size=size, 
                        method="L-BFGS-B", hessian=T,
                        lower=c(-Inf,-Inf),upper=c(Inf, Inf)), silent=T)
      
      ### if-1, "try-error"                 
      if ( class(parm)!="try-error"){
      
         ## give initial value by mean and expectation value
         #th1 <- parm$par[1]
         #th2hat <- a0+a1*th1
         #if ( abs(th2hat+parm$par[2])>8 ){
           #theta <- c(log(mu0), -(a0+log(mu0)*a1))
           #parm <- try(optim(par=theta, fn=nb_llik, x=xk, size=size, 
           #                  method="L-BFGS-B", hessian=T, 
           #                  lower=c(-Inf,-Inf), upper=c(Inf, Inf)), silent=T)
         #}
    
         ### if-3
         #if (class(parm)!="try-error"){
         
         ### if-4, it is convergence
         if ( parm$convergence==0){
            parr <- parm$par
            mu <- exp(parr[1])
            phk <- 1/exp(parr[2])
            bk <- mu
            vk <- mu^2*phk
            sek <- sqrt(diag(pseudoinverse(parm$hessian)))
            res <- c(vk, bk, phk, sek[1], sek[2])
            if (parm$value<=0) res <- c(NA,NA,NA,NA,NA)
         }else{
            res <- c(NA,NA,NA,NA,NA)
         } ### End, if-4
         
         #}else{
         #   res <- c(NA,NA,NA,NA,NA)
         #} ###End, if-3
         
      ###   
      }else{
         res <- c(NA,NA,NA,NA,NA)
      } ### End, if-1
      
      return(res)
   },mc.cores=1)
   tmp <- do.call(rbind, tmp)
  
   ##variance   
   v1 <- tmp[,1]
   v[gene0] <- v1
   ##mean value
   b1 <- tmp[,2]
   b[gene0] <- b1   
   ##dispersion
   phi1 <- tmp[,3]
   phi[gene0] <- phi1
   ##standard error
   se.mu[gene0] <- tmp[,4]
   se.phi[gene0] <- tmp[,5]

   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
   cat(i, bti[i], length(gene0), elapsed, "\n")
   return(list(v, b, phi, se.mu, se.phi, length(gene0)))
})

v <- lapply(TMP,function(ii) ii[[1]])
Vx <- do.call(cbind, v)
colnames(Vx) <- as.character(bti)
save(Vx, file="./10_RNA.Variance_output/tmp9/1_RNA.Vx.RData")

##
b <- lapply(TMP, function(ii) ii[[2]])
Bx <- do.call(cbind, b)
colnames(Bx) <- as.character(bti)
save(Bx, file="./10_RNA.Variance_output/tmp9/1_RNA.Bx.RData")

##
phx <- lapply(TMP, function(ii) ii[[3]])
Phx <- do.call(cbind, phx)
colnames(Phx) <- as.character(bti)
save(Phx, file="./10_RNA.Variance_output/tmp9/1_RNA.Phx.RData")

##
Sx.mu <- lapply(TMP, function(ii) ii[[4]])
Sx.mu <- do.call(cbind, Sx.mu)
colnames(Sx.mu) <- as.character(bti)
save(Sx.mu, file="./10_RNA.Variance_output/tmp9/1_RNA.Sx.mu.RData")

##
Sx.phi <- lapply(TMP, function(ii) ii[[5]])
Sx.phi <- do.call(cbind, Sx.phi)
colnames(Sx.phi) <- as.character(bti)
save(Sx.phi, file="./10_RNA.Variance_output/tmp9/1_RNA.Sx.phi.RData")

###
ngene <- lapply(TMP, function(ii) ii[[6]])
ngene <- do.call(c,ngene)
names(ngene) <- as.character(bti)
save(ngene, file="./10_RNA.Variance_output/tmp9/1_RNA.ngene.RData")

} ###2, End


### 2.2. Calcualte residual dispersion ###
if(FALSE){
rm(list=ls())
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
load("./10_RNA.Variance_output/tmp10/1_RNA.Vx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Bx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Phx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Sx.mu.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Sx.phi.RData")

### (1), Regression
d1 <- melt(Vx)
d2 <- melt(Bx)
d3 <- melt(Phx)
d4 <- melt(Sx.mu)
d5 <- melt(Sx.phi)

ddx <- data.frame(X1=d1$X1, X2=d1$X2,
                 va=d1$value, mu=d2$value, phi=d3$value,
                 se.mu=d4$value, se.phi=d5$value)%>%
       drop_na(va, mu, phi, se.mu, se.phi)#%>%
       #filter(se.mu<3, se.phi<21.81, phi>1.6e-06)

cvt0 <- str_split(ddx$X2, "_", simplify=T)

dd3 <- ddx%>%
       mutate(MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]))%>% 
       mutate(x=log2(mu),y=log2(phi))
       #dplyr::rename(ensgene=X1)%>%
       #group_by(ensgene, MCls, treats)%>%
       #summarise(va=mean(va), mu=mean(mu), phi=mean(phi))%>%
       
       
lmx <- dd3%>%group_by(MCls)%>%
       nest()%>%
       mutate(lmr=map(data, ~lm(y~x,data=.x)))
lmr <- lmx$lmr
names(lmr) <- lmx$MCls
#names(lmr) <- c("Bcell", "Monocyte", "NKcell", "Tcell")
       

### (2), Residual Phx
#Phx[Sx.phi>546.143] <- NA
#Phx[Sx.mu>0.43] <- NA

#Bx[Sx.phi>546.143] <- NA
#Bx[Sx.mu>0.43] <- NA

tmp <- str_split(colnames(Phx), "_", simplify=T)
cvt <- data.frame(rn=colnames(Phx), MCls=tmp[,1])#,Cluster=tmp[,2])

Phx <- log2(Phx) 
Bx <- log2(Bx)
PhxNew <- Phx 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Cluster <- c("3","6") 
for (ii in MCls){       
cvt0 <- cvt%>%filter(MCls==ii)
rn0 <- cvt0[["rn"]]
Phx0 <- Phx[,rn0]  
Bx0 <- Bx[,rn0]
##
lm0 <- coef(lmr[[ii]])
a <- lm0[1]
b1 <- lm0[2]
#b2 <- lm0[3]
#b3 <- lm0[4]
PhxNew[,rn0] <- Phx0-(a+b1*Bx0)#b2*Bx0^2+b3*Bx0^3)
}
PhxNew2 <- 2^PhxNew

save(PhxNew2, file="./10_RNA.Variance_output/tmp10/1_RNA.PhxNew.RData")
}


##############################
### (3) summary parameters ###
##############################

### 3.1, density plot for variance, mean and dispersion
if (FALSE){
load("./10_RNA.Variance_output/tmp9/1_RNA.Vx.RData")
dd1 <- melt(Vx)
load("./10_RNA.Variance_output/tmp9/1_RNA.Bx.RData")
dd2 <- melt(Bx)
load("./10_RNA.Variance_output/tmp9/1_RNA.Phx.RData")
dd3 <- melt(Phx)

load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.mu.RData")
dd4 <- melt(Sx.mu)
load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.phi.RData")
dd5 <- melt(Sx.phi)

dd <- data.frame(X1=dd1$X1,X2=dd1$X2,
                 va=dd1$value, mu=dd2$value, phi=dd3$value,
                 se.mu=dd4$value, se.phi=dd5$value)
ddx <- dd%>%
      drop_na(va, mu, phi,se.mu,se.phi)
      
##variance, 10%, 0.00048
##dispersion, 10%, 0.035
##variance
## at 10% percentage
p1 <- ggplot(ddx, aes(x=log10(va+1e-04)))+
     geom_density()+xlim(-5,5)+
     theme_bw()+ggtitle("Variance")+theme(plot.title=element_text(hjust=0.5))

##mean value
p2 <- ggplot(ddx, aes(x=log10(mu)))+xlim(-5,5)+
     geom_density()+
     theme_bw()+ggtitle("Mean")+theme(plot.title=element_text(hjust=0.5))
        
###dispersion
## 10%  
p3 <- ggplot(ddx, aes(x=log10(phi+1e-02)))+
     geom_density()+xlim(-5,5)+
     theme_bw()+ggtitle("Dispersion")+theme(plot.title=element_text(hjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure1.1.density.png"
png(figfn, width=1000, height=500, res=120)
print(plot_grid(p1, p2, p3, ncol=3))
dev.off()
} ##2.2, End


### 3.2, Distribution of se
if(FALSE){
cat("(2).", "Distribution of se.mu and se.phi", "\n")
load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.mu.RData")
dd1 <- melt(Sx.mu)
load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.phi.RData")
dd2 <- melt(Sx.phi)

cvt0 <- str_split(dd1$X2, "_", simplify=T)
dd <- data.frame(X1=dd1$X1,X2=dd1$X2, 
                 MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                 sampleID=cvt0[,3],  Batch2=cvt0[,4], 
                 se.mu=dd1$value, se.phi=dd2$value)
dd <- dd%>%drop_na(se.mu,se.phi)

#tmp <- dd%>%group_by(MCls)%>%nest()%>%
#       mutate(qq1=map_dbl(data,~quantile((.x)$se.mu, probs=0.9999)),
#              qq2=map_dbl(data,~quantile((.x)$se.phi, probs=0.99)))

## se of mu
fig1 <- ggplot(dd,aes(x=se.mu))+
        geom_histogram(fill="grey70", color="grey40", bins=50)+
        facet_wrap(~MCls, nrow=2, scales="free")+
        theme_bw()
        
figfn <- "./10_RNA.Variance_output/tmp9/Figure1.2.1.seMu.png"
png(figfn, width=600, height=500, res=120)
print(fig1)
dev.off()

## se of phi
fig2 <- ggplot(dd, aes(x=se.phi))+
        geom_histogram(fill="grey70", color="grey40", bins=50)+
        facet_wrap(~MCls, nrow=2, scales="free")+
        theme_bw()
figfn <- "./10_RNA.Variance_output/tmp9/Figure1.2.2.sePhi1.png"
png(figfn, width=600, height=500, res=120)
print(fig2)
dev.off()

## se of phi
## 2, 99%, 117.97
## 3, 95%, 31.79
## 4, 90%, 6.45 
fig2 <- ggplot(dd%>%filter(se.phi<6.45), aes(x=se.phi))+
        geom_histogram(fill="grey70", color="grey40", bins=50)+
        facet_wrap(~MCls, nrow=2, scales="free")+
        theme_bw()
        
figfn <- "./10_RNA.Variance_output/tmp9/Figure1.2.2.sePhi4.png"
png(figfn, width=600, height=500, res=120)
print(fig2)
dev.off()
}


### 3.3, scatter plots, showing relations between mean variance, mean dispersion ###
if (FALSE){
load("./10_RNA.Variance_output/tmp9/1_RNA.Vx.RData")
load("./10_RNA.Variance_output/tmp9/1_RNA.Bx.RData")
load("./10_RNA.Variance_output/tmp9/1_RNA.Phx.RData")
load("./10_RNA.Variance_output/tmp9/1_RNA.PhxNew.RData")
load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.mu.RData")
load("./10_RNA.Variance_output/tmp9/1_RNA.Sx.phi.RData")

d1 <- melt(Vx)
d2 <- melt(Bx)
d3a <- melt(Phx)
d3b <- melt(PhxNew2)
d4 <- melt(Sx.mu)
d5 <- melt(Sx.phi)

dd <- data.frame(X1=d1$X1, X2=d1$X2,
                 va=d1$value, mu=d2$value, phi=d3a$value, phiNew=d3b$value,
                 se.mu=d4$value, se.phi=d5$value)
ddx <- dd%>%
      drop_na(va, mu, phi, se.mu, se.phi)##%>%
      ##filter(se.mu<0.3, se.phi<21.81, phi>1.6e-06)

cvt0 <- str_split(ddx$X2, "_", simplify=T)

dd2 <- ddx%>%
       mutate(MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]))%>% #Cluster=cvt0[,2],
       dplyr::rename(ensgene=X1)

###
fmod <- function(df){
   lm0 <- lm(y~x, data=df)
   r2 <- round(summary(lm0)$r.squared, digits=3)
   #summary(lm0)$r.squared
   eq <- bquote(italic(R)^2==.(r2))
   as.character(as.expression(eq))
}

### label function       
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  as.character(as.expression(eq)) 
} 
#
feq2 <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  r
  #eq <- bquote(italic(R)==.(r)~","~.(symb))
  #as.character(as.expression(eq)) 
} 
             
##figure 1 
dd3 <- dd2%>%mutate(x=log10(mu), y=log10(va))  

anno_df1 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr,feq),
                   r=map_dbl(corr, feq2))
x <- anno_df1%>%dplyr::select(-data, -corr, -corr2)

fig1 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df1, aes(x=-0.5, y=4, label=corr2), parse=T, size=2.5)+  ##CPM 2,9
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        ylab(bquote(log[10]~"(va)"))+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        facet_grid(treats~MCls)+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()+
        theme(legend.title=element_blank())
                 
#fig1 <- ggplot(dd3, aes(x=x, y=y))+
#        geom_point(size=0.05,alpha=0.01)+
        #scale_color_manual("",values=cols1, guide=guide_legend(override.aes=list(size=3)))+
#        xlab(bquote(log["10"]~"("~mu~")"))+
#        ylab(bquote(log["10"]~"(va"~+~1e-08~")"))+          
#        facet_wrap(~MCls,nrow=2)+    
#        geom_smooth(method="loess",formula=y~poly(x,2),span=0.3)+

figfn <- "./10_RNA.Variance_output/tmp9/Figure1.3.1.va_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
#png(figfn, width=600,height=700,res=150)
print(fig1)
dev.off() 

### figure 2, scatter plots show the correlation between dispersion and  mu
dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phi))  ##phi+0.01

anno_df2 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr,feq),
                   r=map_dbl(corr, feq2))
x <- anno_df2%>%dplyr::select(-corr,-corr2)                   
fig2 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df2, aes(x=0.9, y=-2, label=corr2), parse=T, size=2.5)+ ##CPM, 4, -1.5
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.2))+
        facet_grid(treats~MCls)+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()
        
#fig2 <- ggplot(dd2%>%filter(mu.mean<2.6e+06),aes(x=log10(mu.mean), y=log10(phi.mean)))+
#        geom_point(size=0.05)+
#        #scale_color_manual("",values=cols1)+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"("~phi~+~0.001~")"))+
#        geom_smooth(method="lm",formula=y~x)+
#        facet_wrap(~MCls,nrow=2)+
#        theme_bw()+theme(legend.position="none")
figfn <- "./10_RNA.Variance_output/tmp9/Figure1.3.2.phi_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
#png(figfn, width=600, height=700, res=150)
print(fig2)
dev.off()

#### figure 3, scatter plots show the correlation between new dispersion and  mu
dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phiNew))  ##phi+0.01

anno_df3 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr, feq),
                   r=map_dbl(corr, feq2))
                   
x <- anno_df3%>%dplyr::select(-corr,-corr2)
  
fig3 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df3, aes(x=0.5, y=1.5, label=corr2), parse=T, size=2.5)+   #CPM, 3, 1
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.2))+
        facet_grid(treats~MCls)+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()
#        
##fig2 <- ggplot(dd2%>%filter(mu.mean<2.6e+06),aes(x=log10(mu.mean), y=log10(phi.mean)))+
##        geom_point(size=0.05)+
##        #scale_color_manual("",values=cols1)+
##        xlab(bquote(log[10]~"("~mu~")"))+
##        ylab(bquote(log[10]~"("~phi~+~0.001~")"))+
##        geom_smooth(method="lm",formula=y~x)+
##        facet_wrap(~MCls,nrow=2)+
##        theme_bw()+theme(legend.position="none")
figfn <- "./10_RNA.Variance_output/tmp9/Figure1.3.3.phiNew_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
#png(figfn, width=600, height=700, res=150)
print(fig3)
dev.off()
}###

#
####
####new 
#ddx <- dd%>%
#      drop_na(va, mu, phi, phiNew, se.mu, se.phi)%>%
#      filter(se.mu<0.43, se.phi<30)
#
#cvt0 <- str_split(ddx$X2, "_", simplify=T)
#
#dd2 <- ddx%>%
#       mutate(MCls=cvt0[,1],treats=gsub("-EtOH", "", cvt0[,2]))%>%
#       dplyr::rename(ensgene=X1)
#dd2 <- dd2%>%
#       group_by(ensgene, MCls, treats)%>%
#       summarise(va=mean(va), mu=mean(mu), phi=mean(phi), phiNew=mean(phiNew))
#       
#dd3 <- dd2%>%mutate(x=log2(mu),y=log2(phiNew))
#corr <- cor.test(dd3$x, dd3$y)
#       
#fig0 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        geom_text(x=0, y=2.2, label=feq(corr), parse=T, size=3)+
#        xlim(-8,1)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log[2]~"("~mu~")"))+
#        ylab(bquote(log[2]~"("~phi~")"))+          
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()+
#        theme(legend.title=element_blank())
#figfn <- "./10_RNA.Variance_output/tmp6/Figure1.5.4_phi2.mu.png"
#png(figfn, width=550, height=500, res=120)
#print(fig0)
#dev.off()     

### 1.6 plots of average value
#if (FALSE){
#dd2 <- dd2%>%
#       group_by(ensgene, MCls, treats)%>%
#       summarise(va=mean(va), mu=mean(mu), phi=mean(phi),.groups="drop")#,phiNew=mean(phiNew))
#       
#dd3 <- dd2%>%mutate(x=log10(mu), y=log10(va))%>%filter(y>-5)
#
#anno_df1 <- dd3%>%
#            group_by(treats, MCls)%>%
#            nest()%>%
#            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
#                   corr2=map(corr,feq),
#                   r=map_dbl(corr, feq2))
#                   
#fig1 <- ggplot(dd3, aes(x=x, y=y))+
#        geom_point(size=0.1,alpha=0.1)+
#        geom_text(data=anno_df1, aes(x=0.5, y=2.5, label=corr2), parse=T, size=2.5)+       #-9
#        xlab(bquote(log["10"]~"("~mu~")"))+
#        ylab(bquote(log["10"]~"(va)"))+          
#        facet_grid(treats~MCls)+ 
#        geom_smooth(method="lm",formula=y~x, size=0.5)+
#        theme_bw()   
#
#figfn <- "./10_RNA.Variance_output/tmp7/Figure1.7.1.va_mu.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig1)
#dev.off() 
#
#dd3 <- dd2%>%mutate(x=log10(mu), y=log10(phi))%>%filter(y>-2)
#anno_df2 <- dd3%>%
#            group_by(treats, MCls)%>%
#            nest()%>%
#            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
#                   corr2=map(corr,feq),
#                   r=map_dbl(corr, feq2))
#fig2 <- ggplot(dd3,aes(x=x, y=y))+
#        geom_point(size=0.1, alpha=0.1)+
#        geom_text(data=anno_df2, aes(x=0.5, y=1.5, label=corr2), parse=T, size=2.5)+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"("~phi~")"))+
#        geom_smooth(method="lm", formula=y~x, size=0.5)+
#        facet_grid(treats~MCls)+
#        theme_bw()
#figfn <- "./10_RNA.Variance_output/tmp9/Figure1.7.2.phi_mu.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig2)
#dev.off()
#}
 
        
### 3.4, scatter plots, showing relations between mean variance, mean dispersion, log10(CPM) ###
if (FALSE){
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.Vx.RData")
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.Bx.RData")
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.Phx.RData")
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.PhxNew.RData")
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.Sx.mu.RData")
load("./10_RNA.Variance_output/trash/tmp10/1_RNA.Sx.phi.RData")

d1 <- melt(Vx)
d2 <- melt(Bx)
d3a <- melt(Phx)
d3b <- melt(PhxNew2)
d4 <- melt(Sx.mu)
d5 <- melt(Sx.phi)

dd <- data.frame(X1=d1$X1, X2=d1$X2,
                 va=d1$value, mu=d2$value, phi=d3a$value, phiNew=d3b$value,
                 se.mu=d4$value, se.phi=d5$value)
ddx <- dd%>%
      drop_na(va, mu, phi, se.mu, se.phi)##%>%
      ##filter(se.mu<0.3, se.phi<21.81, phi>1.6e-06)

cvt0 <- str_split(ddx$X2, "_", simplify=T)

dd2 <- ddx%>%
       mutate(MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]))%>% #Cluster=cvt0[,2],
       dplyr::rename(ensgene=X1)

###
fmod <- function(df){
   lm0 <- lm(y~x, data=df)
   r2 <- round(summary(lm0)$r.squared, digits=3)
   #summary(lm0)$r.squared
   eq <- bquote(italic(R)^2==.(r2))
   as.character(as.expression(eq))
}

### label function       
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r))#~","~.(symb))
  as.character(as.expression(eq)) 
} 
#
feq2 <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  r
  #eq <- bquote(italic(R)==.(r)~","~.(symb))
  #as.character(as.expression(eq)) 
} 
             
##figure 1 
dd3 <- dd2%>%mutate(x=log10(mu), y=log10(va))  

anno_df1 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr,feq),
                   r=map_dbl(corr, feq2))
x <- anno_df1%>%dplyr::select(-data, -corr, -corr2)

fig1 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df1, aes(x=3, y=9, label=corr2), parse=T, size=2.5)+  ##CPM 2,9
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        ylab(bquote(log[10]~"(va)"))+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        facet_grid(treats~MCls)+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()+
        theme(legend.title=element_blank())
                 
figfn <- "./10_RNA.Variance_output/tmp10/Figure1.3.1.va_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
#png(figfn, width=600,height=700,res=150)
print(fig1)
dev.off() 

### figure 2, scatter plots show the correlation between dispersion and  mu
dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phi))  ##phi+0.01

anno_df2 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr,feq),
                   r=map_dbl(corr, feq2))
x <- anno_df2%>%dplyr::select(-corr,-corr2)                   
fig2 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df2, aes(x=3, y=1.5, label=corr2), parse=T, size=2.5)+ ##CPM, 4, -1.5
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.2))+
        facet_grid(treats~MCls)+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()

figfn <- "./10_RNA.Variance_output/tmp10/Figure1.3.2.phi_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
#png(figfn, width=600, height=700, res=150)
print(fig2)
dev.off()

#### figure 3, scatter plots show the correlation between new dispersion and  mu
dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phiNew))  ##phi+0.01

anno_df3 <- dd3%>%
            group_by(treats, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
                   corr2=map(corr, feq),
                   r=map_dbl(corr, feq2))
                   
x <- anno_df3%>%dplyr::select(-corr,-corr2)

lab1 <- c("CTRL"="CTRL","LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
  
fig3 <- ggplot(dd3, aes(x=x,y=y))+
        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
        geom_text(data=anno_df3, aes(x=3, y=2, label=corr2), parse=T, size=2.5)+   #CPM, 3, 1
        scale_fill_viridis_c()+
        xlab(bquote(log[10]~"("~mu~")"))+
        scale_y_continuous(bquote(log[10]~"("~phi~")"), expand=expansion(mult=0.5))+
        facet_grid(treats~MCls, labeller=labeller(treats=lab1))+          
        #facet_wrap(~MCls, nrow=2, scales="free")+
        geom_smooth(method="lm",formula=y~x, size=0.5)+
        theme_bw()
#        
figfn <- "./10_RNA.Variance_output/tmp10/Figure1.3.3.phiNew_mu.scatter.png"
png(figfn, width=900, height=800, res=150)
print(fig3)
dev.off()
}###       
       


#########################################
### (4) extract protein coding genes  ###
#########################################
if (FALSE){ 
load("./10_RNA.Variance_output/tmp10/1_RNA.Vx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Bx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.Phx.RData")
load("./10_RNA.Variance_output/tmp10/1_RNA.PhxNew.RData")

grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
anno <- data.frame(ensgene=gsub("\\..*", "", rownames(Vx)), ensgene2=rownames(Vx))%>%
        left_join(grch38_unq, by="ensgene")
geneSel <- anno%>%filter(grepl("protein_coding", biotype))%>%dplyr::pull(ensgene2)

Vx <- Vx[geneSel,]
opfn1 <- "./10_RNA.Variance_output/tmp10/1.2_Sel.Vx.RData"
save(Vx, file=opfn1)
Bx <- Bx[geneSel,]
opfn2 <- "./10_RNA.Variance_output/tmp10/1.2_Sel.Bx.RData"
save(Bx, file=opfn2)
Phx <- Phx[geneSel,]
opfn3 <- "./10_RNA.Variance_output/tmp10/1.2_Sel.Phx.RData"
save(Phx, file=opfn3)
PhxNew2 <- PhxNew2[geneSel,]
opfn4 <- "./10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData"
save(PhxNew2, file=opfn4)

}


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


########################################
### 2, differential of gene variance ###
#######################################

if (FALSE){   

#### 2.0, Read data
load("./10_RNA.Variance_output/tmp9/1.2_Sel.Vx.RData")
rn <- rownames(Vx)
rownames(Vx) <- gsub("\\.[0-9]*", "", rn)

###
bti2 <- colnames(Vx)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
comb <- unique(cvt$comb)


### 2.1, estimate differetial results by batch
res <- map_dfr(comb, function(oneX){
   cat(oneX,"\n")
   cvti <- cvt %>%filter(comb==oneX)
   oneComb <- unlist(strsplit(oneX, "_"))
   oneMCl <- oneComb[1]
   oneBatch <- oneComb[2]
   Vi <- Vx[,cvti$bti]
   X <- data.frame(x1=cvti$treats)
   rn <- rownames(Vi)
   ## Start loop By gene 
   TMP <- mclapply(rn, function(ii){
      y <- Vi[ii,]
      nna <- countNA(X$x1,y)
      dd <- myDE(y, X, ii, nna, threshold=3)
      dd
   }, mc.cores=1) ### End loop by gene
   ###  
   TMP <- TMP[!is.na(TMP)]
   TMP <- as.data.frame(do.call(rbind, TMP))%>%mutate(MCls=oneMCl, batch=oneBatch)
   TMP       
})

opfn <- "./10_RNA.Variance_output/tmp9/2_va.results"
write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
}

### 2.2, "Meta analysis"
if(FALSE){
fn <- "./10_RNA.Variance_output/tmp9/2_va.results"
res <- read.table(fn,header=T)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
       filter(batch%in%c("SCAIP1","SCAIP4", "SCAIP5", "SCAIP6"))

## meta
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
        
#opfn <- "./10_RNA.Variance_output/tmp9/2_va.meta"
opfn <- "./10_RNA.Variance_output/tmp9/2_va.meta2" #remove batch 2 and 3
write.table(res3, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")

} ###End, 3.1-3.2, Differential analysis


####
#### 3.3, estimate differetial results all the batch together
#if(FALSE){
#cat("3.0", "Read data", "\n") 
#
#load("./10_RNA.Variance_output/tmp7/1_RNA.Vx.RData")
#rn <- rownames(Vx)
#rownames(Vx) <- gsub("\\.[0-9]*", "", rn)
####
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
#cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
#
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#res <- map_dfr(MCls, function(oneMCl){
#   cat(oneMCl,"\n")   
#   cvti <- cvt %>% filter(MCls==oneMCl)
#   Vi <- Vx[,cvti$bti]
#   X <- data.frame(x1=cvti$treats, x2=cvti$Batch2)
#   rn <- rownames(Vi)
#   ## Start loop By gene 
#   TMP <- mclapply(rn, function(ii){
#      y <- Vi[ii,]
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
#opfn <- "./10_RNA.Variance_output/tmp7/2.1_va2"
#write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
#
#} ###


#######################
### Summary results ###
#######################


if(FALSE){

cat("(1).", "Show qq plots", "\n")
### (1) qq plots
figfn <- "./10_RNA.Variance_output/tmp9/Figure2.1.qq.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300) 
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

res <- read.table("./10_RNA.Variance_output/tmp9/2_va.meta", header=T)
res <- res%>%drop_na(pval) 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){  
##1  
   res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS") 
   print(qq(res1$pval, main="LPS", cex.main=1, cex.axis=0.8, cex.lab=1))
   
##2
   res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")
   print(qq(res2$pval, main="LPS+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
   
##3
   res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")
   print(qq(res3$pval, main="PHA", cex.main=1, cex.axis=0.8, cex.lab=1))
    
##4
   res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")
   print(qq(res4$pval, main="PHA+DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
   
   print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))
}
dev.off() 


### (2), hist distribution of effect size
#cat("(2).", "hist distribution of effect size", "\n")
lab1 <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
fn <- "./10_RNA.Variance_output/tmp9/2_va.meta"
dx <- read.table(fn, header=T)%>%drop_na(beta)
fig0 <- ggplot(dx, aes(x=beta))+
     geom_histogram(fill="grey70", colour="grey20")+
     xlab(bquote("Effective size"~"("~beta~")"))+
     facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
     theme_bw()+
     theme(strip.background=element_blank())

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.2.hist.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig0)
dev.off()

### (3), zscore
#lab1 <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#dx <- read.table("./10_RNA.Variance_output/tmp9/2_va.meta", header=T)%>%
#      mutate(zscore=beta/stderr)%>%
#      drop_na(pval)
#z0 <- quantile(abs(dx$zscore),probs=0.99)
#fig0 <- ggplot(dx%>%filter(abs(zscore)<z0), aes(x=zscore))+
#     geom_histogram(fill="grey70", colour="grey40")+
#     xlab("z score")+
#     facet_grid(MCls~contrast, scales="free_y", labeller=labeller(contrast=lab1))+
#     theme_bw()+
#     theme(strip.background=element_blank())
#
#figfn <- "./10_RNA.Variance_output/tmp7/Figure2.3.hist.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
} ###3.4, End


###################################################################
### (3) show barplot of signifcantly differential gene variance ###
###################################################################

if(FALSE){  
fn <- "./10_RNA.Variance_output/tmp9/2_va.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)   
sigs <- unique(res2$gene)
save(sigs, file="./10_RNA.Variance_output/tmp9/Sig2.DGV.RData")
}


if(FALSE){
#### Barplots show NO.DGV together(Up and down)
fn <- "./10_RNA.Variance_output/tmp9/2_va.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- res2%>%group_by(contrast, MCls)%>%summarise(ngene=n(), .groups="drop")

x <- res2%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")          
fig0 <- ggplot(sigs,aes(x=contrast,y=ngene,fill=MCls))+
        geom_bar(stat="identity",position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(labels=lab2)+ylab("No. DGV")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./10_RNA.Variance_output/tmp9/Figure2.3.1_DGV.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}

### (3), barplots of DEG, up and down with light and deep colors, ***
if(FALSE){

### colors
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
###read data
fn <- "./10_RNA.Variance_output/tmp9/2_va.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)

## up and down DGV
sigs <- res2%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(), .groups="drop")        

###figure2.3.3, facet by contrast, up and down together (stack)  
sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))
ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
### 
fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col2comb, labels="")+ylab("DGV")+ylim(0,1200)+
        geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+50, fill=NULL), size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.3.3_DGV.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### figure2.3.4, facet by cell type,
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col1w <- colorspace::lighten(col1,0.3)
col1comb <- c(col1, col1w)
names(col1comb) <- paste(contrast, rep(c(1,2),each=4), sep="_") 

sig3 <- sigs%>%mutate(comb=paste(contrast, direction, sep="_"))
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
               "PHA"="PHA", "PHA-DEX"="PHA+DEX") 
ann3 <- sig3%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")

### 
fig0 <- ggplot(sig3, aes(x=contrast, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col1comb, labels="")+ylab("DGV")+ylim(0,1200)+
        geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+50, fill=NULL), size=3)+
        scale_x_discrete(labels=lab2)+
        facet_grid(~MCls)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.3.4_DGV.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### figure2.3.5, facet by treatment and Up above x axis and down below x axis
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2,-ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-600,600),5)
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2, fill=comb))+
        geom_bar(stat="identity")+
        scale_fill_manual(values=col2comb, labels="")+
        geom_hline(yintercept=0, color="grey60")+
        geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
                  vjust=ifelse(direction==2, 1, -0.2)), size=3)+ #
        scale_y_continuous("", breaks=breaks_value, limits=c(-600,600),labels=abs(breaks_value))+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.3.5_DGV.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

}###


##################################################################
### (4), correlation of differential effects across conditions ###
##################################################################

if (FALSE){

load("./10_RNA.Variance_output/tmp9/Sig2.DGV.RData")
Geneunq <- unique(sigs)

col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./10_RNA.Variance_output/tmp9/2_va.meta2"
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
y0 <- y[abs(y)<6.49] #99% percent quantile(abs(y),probs=0.99)
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

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.4.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "Tcell_LPS+DEX", "Tcell_PHA+DEX", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
              "Monocyte_LPS", "Monocyte_PHA", "NKcell_LPS", "NKcell_PHA",
              "Bcell_LPS", "Bcell_PHA", "Tcell_LPS", "Tcell_PHA") 

corr <- cor(TMP0)#[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)

x <- str_split(colnames(corr), "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)
tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
                 cluster_rows=T, cluster_cols=T,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F, na_col="white")

figfn <- "./10_RNA.Variance_output/tmp9/Figure2.4.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()

#mycol <- viridisLite::viridis(100)

#figfn <- "./10_RNA.Variance_output/tmp6/Figure2.4.2_corr.beta.png"
#png(figfn, width=1000, height=1000, res=180)
#print(corrplot(corr, method="color", order="hclust", hclust.method="complete", col=mycol,
#         tl.col="black", tl.cex=0.8, outline=F, diag=T))
#dev.off()

} ##End, 5




##################################
### 3, differential dispersion ###
##################################

###############################
### Differential dispersion ###
###############################
if (FALSE){   

cat("Differential dispersion", "\n")
#### Read data
load("./10_RNA.Variance_output/tmp9/1.2_Sel.Phx.RData")
rn <- rownames(Phx)
rownames(Phx) <- gsub("\\.[0-9]*", "", rn)

###
bti2 <- colnames(Phx)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
cvt <- cvt%>%mutate(comb=paste(MCls, Batch2, sep="_"))
comb <- unique(cvt$comb)


### 4.1, estimate differetial results by batch
#cat("4.1", "Differential analysis by batch", "\n")
res <- map_dfr(comb, function(oneX){
   cat(oneX,"\n")
   cvti <- cvt %>%filter(comb==oneX)
   oneComb <- unlist(strsplit(oneX, "_"))
   oneMCl <- oneComb[1]
   oneBatch <- oneComb[2]
   Phi <- Phx[,cvti$bti]
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
opfn <- "./10_RNA.Variance_output/tmp9/3_phi.results"
write.table(res, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")
}

### 4.2, "Meta analysis"
#cat("4.2", "Meta analysis", "\n")
if(FALSE){
fn <- "./10_RNA.Variance_output/tmp9/3_phi.results"
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
        
#opfn <- "./10_RNA.Variance_output/tmp9/3_phi.meta"
opfn <- "./10_RNA.Variance_output/tmp9/3_phi.meta2" ##remove batch 2 and 3
write.table(res3, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")

} ###End, Differential analysis


#######################
### Summary results ###
#######################

#####################
### (1), qq plots ###
#####################

if(FALSE){

cat("(1).", "Show qq plots", "\n")
figfn <- "./10_RNA.Variance_output/tmp9/Figure3.1.qq.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300)
par(mfrow=c(4,8),mar=c(4,4,1.5,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

res <- read.table("./10_RNA.Variance_output/tmp9/3_phi.meta", header=T)%>%drop_na(pval)
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")       
for (oneMCl in MCls){  
##1  
   res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS") 
   print(qq(res1$pval, main="LPS", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))

##2
   res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")
   print(qq(res2$pval, main="LPS-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
##3
   res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")
   print(qq(res3$pval, main="PHA", cex.main=0.8, cex.axis=0.8, cex.lab=0.8)) 
   
##4
   res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")
   print(qq(res4$pval, main="PHA-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
   print(mtext(oneMCl, side=4, line=0.5, cex=0.8, col="blue") )
}
dev.off() 

Sys.sleep(5)


############################
### (2), histogram plots ###
############################

rm(list=ls())
cat("(2).", "hist plots for differential effects size", "\n")

dx <- read.table("./10_RNA.Variance_output/tmp9/3_phi.meta",header=T)%>%drop_na(beta)
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fig0 <- ggplot(dx, aes(x=beta))+
     geom_histogram(fill="grey70", color="grey20")+
     xlab(bquote("Effective size"~"("~beta~")"))+
     facet_grid(MCls~contrast,scales="free_y", labeller=labeller(contrast=lab2))+
     theme_bw()+
     theme(strip.background=element_blank())

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.2.hist.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig0)
dev.off()

Sys.sleep(5)

} ##4.3,4.4, End


#######################
### (3), Barplots   ###
#######################

if(FALSE){
res <- read.table("./10_RNA.Variance_output/tmp9/3_phi.meta",header=T)%>%
       drop_na(qval)%>%
       filter(qval<0.1, abs(beta)>0.5)
sigs <- unique(res$gene)
save(sigs, file="./10_RNA.Variance_output/tmp9/Sig3.DGP.RData")
}

if(FALSE){
#### Barplots show NO.DGV together(Up and down)
fn <- "./10_RNA.Variance_output/tmp9/3_phi.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- res2%>%group_by(contrast, MCls)%>%summarise(ngene=n(), .groups="drop")

x <- res2%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")          
fig0 <- ggplot(sigs,aes(x=contrast,y=ngene,fill=MCls))+
        geom_bar(stat="identity",position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(labels=lab2)+ylab("No. DGP")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./10_RNA.Variance_output/tmp9/Figure3.3.1_DGP.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}


### (3), barplots of DGP, up and down with light and deep colors, ***
if(FALSE){

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
                   
###read data
fn <- "./10_RNA.Variance_output/tmp9/3_phi.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)

## up and down DGV
sigs <- res2%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(), .groups="drop")             


###Figure3.3.3,  facet by contrast and stacked up and down above x axis        
sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))
ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))

### 
fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col2comb, labels="")+ylab("DGP")+ylim(0,600)+
        geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+20, fill=NULL), size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.3.3_DGP.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### Figure3.3.4, facet by cell type, stacked up and down
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
        scale_fill_manual(values=col1comb, labels="")+ylab("DGP")+ylim(0,600)+
        geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+20, fill=NULL), size=3)+
        facet_grid(~MCls)+
        scale_x_discrete(labels=lab2)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.3.4_DGP.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### Figure3.3.5, facet by contrast, and up above axis and down below axis
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2, -ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-400,300),5)

facetlab <- as_labeller(lab2)
fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2, fill=comb))+
        geom_bar(stat="identity")+
        scale_fill_manual(values=col2comb, labels="")+
        geom_hline(yintercept=0, color="grey60")+
        geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
                  vjust=ifelse(direction==2, 1.1, -0.2)), size=3)+ #
        scale_y_continuous("", breaks=breaks_value, limits=c(-400,300),labels=abs(breaks_value))+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.3.5_DGP.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


}###


##################################################################
### (4), correlation of differential effects across conditions ###
##################################################################

if (FALSE){

load("./10_RNA.Variance_output/tmp9/Sig3.DGP.RData")
Geneunq <- sigs
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
          
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./10_RNA.Variance_output/tmp9/3_phi.meta2"
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
y0 <- y[abs(y)<6.09] #99% percent quantile(abs(y),probs=0.99)
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

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.4.1_heatmap.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig1)
dev.off()                     


### (2) correlation heatmap ###
#Neworder <- c("Bcell_LPS+DEX", "Bcell_LPS", "Bcell_PHA+DEX",  "Bcell_PHA",
#              "Monocyte_LPS+DEX", "Monocyte_LPS", "Monocyte_PHA+DEX", "Monocyte_PHA", 
#              "Tcell_LPS+DEX","Tcell_LPS", "Tcell_PHA+DEX", "Tcell_PHA", 
#               "NKcell_LPS+DEX", "NKcell_LPS", "NKcell_PHA+DEX", "NKcell_PHA") 
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
#              "Tcell_LPS+DEX", "Tcell_PHA+DEX", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
#              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#               "Tcell_LPS", "Tcell_PHA", "NKcell_LPS", "NKcell_PHA") 
Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", 
              "NKcell_LPS+DEX", "Tcell_LPS+DEX",
              "NKcell_PHA+DEX", "Tcell_PHA+DEX", 
              "Bcell_LPS+DEX", "Bcell_PHA+DEX",               
              "Monocyte_LPS", "Monocyte_PHA", 
              "Bcell_LPS", "Bcell_PHA", 
              "NKcell_PHA", "Tcell_PHA",
              "NKcell_LPS", "Tcell_LPS") 

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

figfn <- "./10_RNA.Variance_output/tmp9/Figure3.4.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()          
                   
} ##End, 5


#################################################################
### 3 (2) Differential dispersion after removing mean effects ###
#################################################################


####################################
### Differential procedure ###
####################################

if(FALSE){

cat("Differential residual dispertion","\n")
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
} ###


### meta analysis
if (FALSE){
### "Meta analysis"
cat("Meta analysis", "\n")
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
if(FALSE){
rm(list=ls())
cat("Show qq plots", "\n")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.1.qq.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300)
par(mfrow=c(4,8),mar=c(4,4,1.5,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%drop_na(pval) 
for (oneMCl in MCls){   
##1  
   res1 <- res %>% filter(MCls==oneMCl, contrast=="LPS") 
   print(qq(res1$pval, main="LPS", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))

##2
   res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")
   print(qq(res2$pval, main="LPS-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
##3
   res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")
   print(qq(res3$pval, main="PHA", cex.main=0.8, cex.axis=0.8, cex.lab=0.8)) 
   
##4
   res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")
   print(qq(res4$pval, main="PHA-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
   print(mtext(oneMCl, side=4, line=0.5, cex=0.8, col="blue") )
}
dev.off() 

Sys.sleep(5)

############################
### (2), histogram plots ###
############################

rm(list=ls())
cat("hist plots for differential effects size", "\n")
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
dx <- read.table(fn,header=T)%>%drop_na(beta) 
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fig0 <- ggplot(dx, aes(x=beta))+
     geom_histogram(fill="grey70", color="grey20")+
     xlab(bquote("Effective size"~"("~beta~")"))+
     facet_grid(MCls~contrast,scales="free",labeller=labeller(contrast=lab2))+
     theme_bw()+
     theme(strip.background=element_blank())

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.2.hist.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig0)
dev.off()

Sys.sleep(5)

} ##4.3,4.4, End


############################
### (3), Barplots of DGP ###
############################

if(FALSE){
##
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- unique(res$gene)
save(sigs, file="./10_RNA.Variance_output/tmp9/Sig3x.DGP.RData")
} 


if(FALSE){
#### Barplots show NO.DGV together(Up and down)
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(file=fn,header=T)
res2 <- res%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
sigs <- res2%>%group_by(contrast, MCls)%>%summarise(ngene=n(), .groups="drop")

x <- res2%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")          
fig0 <- ggplot(sigs,aes(x=contrast,y=ngene,fill=MCls))+
        geom_bar(stat="identity",position=position_dodge())+
        scale_fill_manual(values=cols)+
        scale_x_discrete(labels=lab2)+ylab("No. DVG")+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.title.x=element_blank())
###
figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.1_DGP.barplot.png"
png(filename=figfn, width=600, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}


### (3), barplots of DGP, up and down with light and deep colors, ***
if(FALSE){

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

## up and down DGV
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
        
        
### Figure3x.3.3, facet by contrast and up and down together(stack)       
sig2 <- sigs%>%mutate(comb=paste(MCls, direction, sep="_"))
ann2 <- sig2%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
#xpos <- c("Bcell"=0.5,"Monocyte"=1, "NKcell"=1.5, "Tcell"=2)
#ann2$xpos <- xpos[ann2$MCls]
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
 
fig0 <- ggplot(sig2, aes(x=MCls, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col2comb, labels="")+ylab("DVG")+ylim(0,600)+
        geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+20, fill=NULL), size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.3_DGP.barplot3.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()


### Figure3x.3.4, facet by cell type, up and down together(stack)
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col1w <- colorspace::lighten(col1,0.3)
col1comb <- c(col1, col1w)
names(col1comb) <- paste(contrast, rep(c(1,2),each=4), sep="_") 

sig3 <- sigs%>%mutate(comb=paste(contrast, direction, sep="_"))
sig3$facet_fill_color <- col2[sig3$MCls]
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
                
ann3 <- sig3%>%group_by(MCls, contrast)%>%summarise(ngene=sum(ngene),.groups="drop")
### 
fig0 <- ggplot(sig3, aes(x=contrast, y=ngene, fill=comb))+
        geom_bar(stat="identity", position="stack")+ 
        scale_fill_manual(values=col1comb, labels="")+ylab("DVG")+ylim(0,600)+
        geom_text(data=ann3, aes(x=contrast, label=ngene, y=ngene+20, fill=NULL), size=3)+
        scale_x_discrete(labels=lab2)+
        facet_grid(~MCls)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.4_DGP.barplot4.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
} ###


###Figure3x.3.5,  barplots of DVG, facet by contrast, and up above axis and down below axis ***
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

if(FALSE){
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
              
fig0 <- fig0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)        

figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.3.5_DGP.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()

}###


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
           
}


##################################################################
### (4), correlation of differential effects across conditions ###
##################################################################

if (FALSE){

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

} ##End, 5


############################################################
### (5), compare beta (differential residual dispersion) ###
###         from LPS/PHA and LPS-DEX/PHA-DEX             ###
############################################################
#
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
        scale_x_continuous("LPS effect size on dispersion", expand=expansion(mult=0.1))+
        scale_y_continuous("LPS+DEX effect size on dispersion", expand=expansion(mult=0.1))+
        theme_bw()+
        theme(strip.background=element_blank(),
              axis.title=element_text(size=10))
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
        scale_x_continuous("PHA effect size on dispersion", expand=expansion(mult=0.1))+
        scale_y_continuous("PHA+DEX effect size on dispersion", expand=expansion(mult=0.1))+
        theme_bw()+
        theme(strip.background=element_blank(),
              axis.title=element_text(size=10))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)        
                           
figfn <- "./10_RNA.Variance_output/tmp9/Figure3x.5.2_PHA.png"
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig2)
dev.off()
} ###


#####################################
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
if(FALSE){

cat("Show qq plots", "\n")

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.1.qq.png"
png(figfn, width=2000, height=2000, pointsize=12, res=300)
par(mfrow=c(4,8),mar=c(4,4,1.5,2),mgp=c(2,1,0))
x <- matrix(1:16, 4, 4, byrow=T)
layout(x)

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn,header=T)%>%drop_na(pval)
for (oneMCl in MCls){
##1  
   res1 <- res %>%filter(MCls==oneMCl, contrast=="LPS") 
   print(qq(res1$pval, main="LPS", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))

##2
   res2 <- res %>% filter(MCls==oneMCl, contrast=="LPS-DEX")
   print(qq(res2$pval, main="LPS-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
##3
   res3 <- res %>% filter(MCls==oneMCl, contrast=="PHA")
   print(qq(res3$pval, main="PHA", cex.main=0.8, cex.axis=0.8, cex.lab=0.8)) 
   
##4
   res4 <- res %>% filter(MCls==oneMCl, contrast=="PHA-DEX")
   print(qq(res4$pval, main="PHA-DEX", cex.main=0.8, cex.axis=0.8, cex.lab=0.8))
   
   print(mtext(oneMCl, side=4, line=0.5, cex=0.8, col="blue") )
}
dev.off() 

Sys.sleep(5)

############################
### (2), histogram plots ###
############################
rm(list=ls())
cat("hist plots for differential effects size", "\n")
lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
dx <- read.table(fn,header=T)%>%drop_na(beta) 
fig0 <- ggplot(dx, aes(x=beta))+
     geom_histogram(fill="grey70", color="grey20")+
     xlab(bquote("Effective size"~"("~beta~")"))+
     facet_grid(MCls~contrast, scales="free",labeller=labeller(contrast=lab2))+
     theme_bw()+
     theme(strip.background=element_blank())

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.2.hist.png"
png(filename=figfn, width=800, height=800, pointsize=12, res=130)  
print(fig0)
dev.off()

Sys.sleep(5)

} ##5.3,5.4, End



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

### Figure4.3.5, facet by contrast, and up above axis and down below axis
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

if(FALSE){
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
        scale_y_continuous("", breaks=breaks_value, limits=c(-800,800),labels=abs(breaks_value))+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
fig0 <- fig0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)        

figfn <- "./10_RNA.Variance_output/tmp9/Figure4.3.5_DMG.barplot5.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)  
print(fig0)
dev.off()
}###


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

if (FALSE){

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

} ##End, 5


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



##################
### 6, overlap ###
##################

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

if (FALSE){
### Read data
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
mycol <- c("1"="grey20", "2"="red", "3"="blue", "4"="#9400D3")
#mycol <- c("1"="#bababa", "2"="#de2d26", "3"="#6baed6", "4"="#756bb1")
Newcon2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")

#### (1). mean vs residual dispersion
cat("(1).", "mean vs residual dispersion", "\n")
dfxy1 <- myDFxy(df3, df2)
mylabel <- c("1"="NS", "2"="DMG(only)", "3"="DVG(only)", "4"="Both")

anno_df1 <- dfxy1%>%
            group_by(contrast, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
                   eq=map(corr,feq))%>%
           dplyr::select(-data,-corr)
        
fig1 <- ggplot(dfxy1, aes(x=zscore.x, y=zscore.y))+
        geom_point(aes(colour=factor(gr)), size=0.3)+
        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
        geom_hline(yintercept=0, color="grey60")+
        geom_vline(xintercept=0, color="grey60")+
        geom_text(data=anno_df1, x=0, y=22, aes(label=eq), size=3, parse=T)+
        xlab("zscore of gene mean expression")+ylab("zscore of gene variability")+
        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.text=element_text(size=9),
              #legend.key.size=unit(0.4,units="cm"),
              axis.text=element_text(size=9))
###              
figfn <- "./10_RNA.Variance_output/tmp9/Figure5.1_DMGvsDVG.png"
png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
print(fig1)
dev.off() 
###sumamary number
x <- dfxy1%>%
     filter(gr==4)%>%
     group_by(MCls, contrast)%>%nest()%>%
     mutate(ngene=map_dbl(data,nrow))

### (2). mean and gene expression
cat("(2).", "mean vs gene expression", "\n")
dfxy2 <- myDFxy(df3, df1)
mylabel <- c("1"="NS", "2"="DMG(only)", "3"="DEG(only)", "4"="Both")
anno_df2 <- dfxy2%>%
            group_by(contrast, MCls)%>%
            nest()%>%
            mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
                   eq=map(corr,feq))%>%
            dplyr::select(-data,-corr)
        
fig2 <- ggplot(dfxy2, aes(x=zscore.x, y=zscore.y))+
        geom_point(aes(colour=factor(gr)), size=0.3)+
        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
        geom_text(data=anno_df2, x=0, y=50, aes(label=eq), size=3, parse=T)+
        geom_hline(yintercept=0, color="grey60")+
        geom_vline(xintercept=0, color="grey60")+
        xlab("zscore of gene mean")+ylab("zscore of gene expression")+
        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
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

} ### 6.2, End


if(TRUE){
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
sigs <- sigs%>%mutate(comb=paste(MCls, gr, sep="_"))
#sig2 <- sigs%>%group_by(MCls, contrast)%>%mutate(prop=ngene/sum(ngene))

lab2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
colw <- lapply(MCls, function(ii){
         x1 <- colorspace::lighten(col2[ii], 0)
         x2 <- colorspace::lighten(col2[ii], 0.6)
         x3 <- colorspace::lighten(col2[ii], 0.3)
         xx <- c(x1, x2, x3)
         names(xx) <- paste(ii, 1:3, sep="_")
         xx
         })
colw <- unlist(colw)                                             

          
fig0 <- ggplot(sigs,aes(x=MCls, y=ngene))+
        geom_bar(stat="identity", position=position_stack(reverse=T), aes(fill=comb))+
        scale_fill_manual(values=colw)+
        #geom_text(aes(label=ngene), position="stack", hjust=0.5, vjust=3, size=3)+
        facet_grid(~contrast, labeller=facetlab)+
        theme_bw()+
        theme(legend.position="none",
              axis.title=element_blank(),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
legend2 <- get_legend(
        ggplot(sigs%>%filter(MCls=="Monocyte"),aes(x=contrast,y=ngene))+
        geom_bar(stat="identity", position=position_stack(), aes(fill=comb))+
        scale_fill_manual(values=colw[grepl("Monocyte",names(colw))], 
                          labels=c("Monocyte_1"="DEG", "Monocyte_2"="Both", "Monocyte_3"="DVG"))+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.background=element_blank(),
              legend.text=element_text(size=8),
              legend.key.size=grid::unit(1,"lines")))

###
figfn <- "./10_RNA.Variance_output/tmp9/Figure5.3_bar.png"
png(filename=figfn, width=800, height=400, pointsize=12, res=120)
print(plot_grid(fig0, legend2, rel_widths=c(4,0.5)))
dev.off()  
}
      


#########################################
### 3.7, show plots for specific gene ###
#########################################
#if(FALSE){
#
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Phx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.PhxNew.RData")
#
#count <- colSums(YtX)
#count_after <- median(count)
#count <- count/count_after
#
#
#bti2 <- colnames(YtX)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                  sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#bti2 <- colnames(PhxNew2)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                   sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#
#
#lab1 <- c("CTRL"="CTRL", 
#          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
#          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
##col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
#
#
#
##Naii <- is.na(Vx)
##Num <- rowSums(Naii)
##Vtmp <- Vx[Num==0,]
##range(Vtmp)
##cols <- pal_npg("nrc",alpha=0.6)(5)
##cols1 <- c("CTRL"="#828282", 
##           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
##           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")   
##
#### fig2, hist plot for variance
##gene0 <- "ENSG00000169442"
##cvt$v <- Vx[gene0,]
##dd1 <- cvt %>% filter(!is.na(v))
##fig1 <- ggplot(dd1,aes(x=log(v)))+
##        geom_histogram(color="#e9ecef", alpha=0.6, position="identity")+
##        xlab("log(variance)")+facet_grid(MCls~treats)+theme_bw()
#####
##figfn <- paste("./10_RNA.Variance_output/Figure2.1_va.", gene0, ".hist.log.pdf", sep="")
##pdf(figfn)
##fig1
##dev.off() 
##"ENSG00000187608"
##"ENSG00000126709"
##"ENSG00000089127"
##"ENSG00000111331"
##"ENSG00000147434"
#
#gene0 <- "ENSG00000157601"
#symbol <- "MX1"
###
#cvt$y <- YtX[grepl(gene0,rownames(YtX)),]/count 
#d1 <- cvt%>%drop_na(y)%>%filter(MCls=="Tcell")
#fig1 <- ggplot(d1, aes(x=treats, y=log2(y), colour=treats))+
#        geom_boxplot()+ylab("Expression")+
#        scale_colour_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("",labels=lab1)+
#        ggtitle(symbol)+
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0),
#              plot.title=element_text(hjust=0.5,face="italic", size=8), 
#              legend.position="none")
#        
#        
#
#### fig2. box plot for variance
#cvt2$y <- PhxNew2[grepl(gene0,rownames(PhxNew2)),]
#d2 <- cvt2%>%drop_na(y)%>%filter(MCls=="Tcell")
#fig2 <- ggplot(d2, aes(x=treats, y=log2(y+1e-02), colour=treats))+
#        geom_boxplot()+ylab(bquote(log[2]~"("~italic(phi)~")"))+
#        scale_colour_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("",labels=lab1)+
#        ggtitle(symbol)+ 
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0, size=6),
#              axis.title.y=element_text(size=8),
#              plot.title=element_text(hjust=0.5, face="italic", size=8),
#              legend.position="none")
#              #plot.title = element_text(hjust = 0.5))
#              
#
#figfn <- paste("./10_RNA.Variance_output/tmp6/Example/Figure6.", gene0, ".box.png", sep="")
##png(figfn, width=600, height=400, res=120)
##print(plot_grid(fig1, fig2, ncol=2)) 
#png(figfn, width=200, height=200, res=100)
#print(fig2) 
#dev.off()
#
#}



#########################
### (5). gene example ###
#########################
#if (FALSE){
#
#lab1 <- c("CTRL"="CTRL", 
#          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
#          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#           
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                   sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#gene0 <- "ENSG00000140043"
#symbol0 <- "PTGR2"
#
#cvt2$y <- Vx[grepl(gene0,rownames(Vx)),]
#d2 <- cvt2%>%drop_na(y)
#
#fig2 <- ggplot(d2%>%filter(MCls=="Tcell"), aes(x=treats, y=log2(y), colour=treats))+
#        geom_boxplot()+ylab(bquote(log[2]~"(variability)"))+
#        scale_colour_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("", labels=lab1)+
#        ggtitle(bquote(~italic(.(symbol0))~" variability in T cell"))+
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
#              axis.title.y=element_text(size=10),
#              plot.title=element_text(hjust=0.5, size=10),
#              legend.position="none")
#
#figfn <- paste("./10_RNA.Variance_output/tmp6/Figure2.5.2_", gene0, ".Tcell.Boxplot.png", sep="")
#png(figfn, width=400, height=600, res=120)
#print(fig2)  
#dev.off()
#
#}
#
#### B, Barplots show DGV, up and down separately, meanwhile with denotions of significance 
#cat("(B).","DGV up and  down separately", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#
#tmp2 <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./10_RNA.Variance_output/tmp6/2_va.", oneMCl, ".meta", sep="")
#   res <- read.table(file=fn,header=T)
#   res0 <- res%>%filter(qval<0.1, abs(beta)>0.5, !is.na(qval)) 
#})
### up and down DGV
#sigs <- tmp2%>%
#        mutate(direction=ifelse(beta>0, "1", "2"))%>%
#        group_by(contrast, MCls, direction)%>%
#        summarise(ngene=n())
#
#
#fmod <- function(subdf) {
#   res <- binom.test(subdf$ngene, 0.5, alternative="two.sided")
#   res$p.value
#}           
#anno_df <- sigs%>%
#              group_by(contrast, MCls)%>%
#              nest()%>%
#              mutate(pval=map_dbl(data,fmod),
#                     y=map_dbl(data,~max((.x)$ngene)))%>%
#              unnest()%>%
#              distinct(contrast, MCls,.keep_all=T)
#xpos <- c("LPS"=0.8, "LPS-DEX"=1.8, "PHA"=2.8, "PHA-DEX"=3.8)
#xmin <- xpos[anno_df$contrast]
#anno_df$xmin <- xmin
#anno_df <- anno_df%>%
#           mutate(xmax=xmin+0.4, y1=y+50)
#
##anno_df <- anno_df%>%mutate(xmin=contrast,xmax=contrast, y1=y+50)
#label <- rep("*", nrow(anno_df))
#label[anno_df$pval<0.01] <- "**"
#label[anno_df$pval<0.001] <- "***"
#anno_df$label <- label
#
#anno_df <- anno_df%>%filter(pval<0.05)
#anno_df$group <- 1:nrow(anno_df)
#
#
##uplab <- c("up", "down")
##names(uplab) <- c("1", "2")
##cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#377eb8", "Tcell"="#e41a1c")
#cols <- c("1"="red","2"="blue")
#mylab <- c("1"="up","2"="down")
#               
#fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=direction))+geom_bar(stat="identity",position="dodge")+  #"stack"
#        scale_fill_manual(values=cols,labels=mylab)+
#        geom_signif(data=anno_df, 
#                    aes(xmin=xmin, xmax=xmax, annotations=label, y_position=y1, group=group),
#                    vjust=0.1, tip_length=0.05, manual=T)+ 
#        facet_wrap(~factor(MCls),ncol=2)+
#        xlab("")+ylab("No. DGV")+ylim(0,900)+
#        theme_bw()+
#        theme(strip.background=element_blank(),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9),
#              strip.text.x=element_text(size=12))
####
#figfn <- "./10_RNA.Variance_output/tmp6/Figure2.3.2_DGV.barplot2.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
#}  ###

#### (2), Barplots show No. DGP, up and down separately, with significant denotion
#cat("(2).","DGV up and down", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#
#tmp2 <- map_dfr(MCls, function(oneMCl){   
#   fn <- paste("./10_RNA.Variance_output/tmp6/3_phi.", oneMCl, ".meta", sep="")
#   res <- read.table(file=fn,header=T)
#   res0 <- res%>%filter(qval<0.1, abs(beta)>0.5, !is.na(qval))   
#})
#
### up and down DGV
#sigs <- tmp2%>%
#        mutate(direction=ifelse(beta>0, "1", "2"))%>%
#        group_by(contrast, MCls, direction)%>%
#        summarise(ngene=n())
#
#
#fmod <- function(subdf) {
#   res <- binom.test(subdf$ngene, 0.5, alternative="two.sided")
#   res$p.value
#}           
#anno_df <- sigs%>%
#              group_by(contrast, MCls)%>%
#              nest()%>%
#              mutate(pval=map_dbl(data,fmod),
#                     y=map_dbl(data,~max((.x)$ngene)))%>%
#              unnest()%>%
#              distinct(contrast, MCls,.keep_all=T)
#xpos <- c("LPS"=0.8, "LPS-DEX"=1.8, "PHA"=2.8, "PHA-DEX"=3.8)
#xmin <- xpos[anno_df$contrast]
#anno_df$xmin <- xmin
#anno_df <- anno_df%>%
#           mutate(xmax=xmin+0.4, y1=y+50)
#
##anno_df <- anno_df%>%mutate(xmin=contrast,xmax=contrast, y1=y+50)
#label <- rep("*", nrow(anno_df))
#label[anno_df$pval<0.01] <- "**"
#label[anno_df$pval<0.001] <- "***"
#anno_df$label <- label
#
#anno_df <- anno_df%>%filter(pval<0.05)
#anno_df$group <- 1:nrow(anno_df)
#
#
##uplab <- c("up", "down")
##names(uplab) <- c("1", "2")
##cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#377eb8", "Tcell"="#e41a1c")
#cols <- c("1"="red","2"="blue")
#mylab <- c("1"="up","2"="down")
#               
#fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=direction))+geom_bar(stat="identity",position="dodge")+  #"stack"
#        scale_fill_manual(values=cols,labels=mylab)+
#        geom_signif(data=anno_df, 
#                    aes(xmin=xmin, xmax=xmax, annotations=label, y_position=y1, group=group),
#                    vjust=0.1, tip_length=0.05, manual=T)+ 
#        facet_wrap(~factor(MCls),ncol=2)+
#        xlab("")+ylab("No. DGP")+ylim(0,600)+
#        theme_bw()+
#        theme(strip.background=element_blank(),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9),
#              strip.text.x=element_text(size=12))
####
#figfn <- "./10_RNA.Variance_output/tmp6/Figure3.4_DGP.barplot.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
#}

#### (2), Barplots show No. DGP, up and down separately, with significant denotion
#cat("(2).","DGV up and down", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#
#tmp2 <- map_dfr(MCls, function(oneMCl){
#   fn <- paste("./10_RNA.Variance_output/tmp6/3_phiNew.", oneMCl, ".meta", sep="")
#   res <- read.table(file=fn,header=T)
#   res0 <- res%>%filter(qval<0.1, abs(beta)>0.5, !is.na(qval)) 
#})
#
### up and down DGV
#sigs <- tmp2%>%
#        mutate(direction=ifelse(beta>0, "1", "2"))%>%
#        group_by(contrast, MCls, direction)%>%
#        summarise(ngene=n())
#
#
#fmod <- function(subdf) {
#   res <- binom.test(subdf$ngene, 0.5, alternative="two.sided")
#   res$p.value
#}           
#anno_df <- sigs%>%
#              group_by(contrast, MCls)%>%
#              nest()%>%
#              mutate(pval=map_dbl(data,fmod),
#                     y=map_dbl(data,~max((.x)$ngene)))%>%
#              unnest()%>%
#              distinct(contrast, MCls,.keep_all=T)
#xpos <- c("LPS"=0.8, "LPS-DEX"=1.8, "PHA"=2.8, "PHA-DEX"=3.8)
#xmin <- xpos[anno_df$contrast]
#anno_df$xmin <- xmin
#anno_df <- anno_df%>%
#           mutate(xmax=xmin+0.4, y1=y+50)
#
##anno_df <- anno_df%>%mutate(xmin=contrast,xmax=contrast, y1=y+50)
#label <- rep("*", nrow(anno_df))
#label[anno_df$pval<0.01] <- "**"
#label[anno_df$pval<0.001] <- "***"
#anno_df$label <- label
#
#anno_df <- anno_df%>%filter(pval<0.05)
#anno_df$group <- 1:nrow(anno_df)
#
#
##uplab <- c("up", "down")
##names(uplab) <- c("1", "2")
##cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#377eb8", "Tcell"="#e41a1c")
#cols <- c("1"="red","2"="blue")
#mylab <- c("1"="up","2"="down")
#               
#fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=direction))+geom_bar(stat="identity",position="dodge")+  #"stack"
#        scale_fill_manual(values=cols, labels=mylab)+
#        geom_signif(data=anno_df, 
#                    aes(xmin=xmin, xmax=xmax, annotations=label, y_position=y1, group=group),
#                    vjust=0.1, tip_length=0.05, manual=T)+ 
#        facet_wrap(~factor(MCls),ncol=2)+
#        xlab("")+ylab("No. DGP")+ylim(0,600)+
#        theme_bw()+
#        theme(strip.background=element_blank(),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9),
#              strip.text.x=element_text(size=12))
####
#figfn <- "./10_RNA.Variance_output/tmp6/Figure3x.3.2_DGP.barplot2.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
#}  ###
#
#### (3), Barplots show DGE, up and down separately, with significant denotions
#cat("(3).","DGV up and down", "\n")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#
#
#fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
#res <- read.table(fn, header=T)%>%drop_na(qval)%>%filter(qval<0.1, abs(beta)>0.5)
#
### up and down DGV
#sigs <- res%>%
#        mutate(direction=ifelse(beta>0, "1", "2"))%>%
#        group_by(contrast, MCls, direction)%>%
#        summarise(ngene=n())
#
#          
#anno_df <- sigs%>%
#              group_by(contrast, MCls)%>%
#              nest()%>%
#              mutate(pval=map_dbl(data,fmod),
#                     y=map_dbl(data,~max((.x)$ngene)))%>%
#              unnest()%>%
#              distinct(contrast, MCls,.keep_all=T)
#xpos <- c("LPS"=0.8, "LPS-DEX"=1.8, "PHA"=2.8, "PHA-DEX"=3.8)
#xmin <- xpos[anno_df$contrast]
#anno_df$xmin <- xmin
#anno_df <- anno_df%>%
#           mutate(xmax=xmin+0.4, y1=y+50)
#
##anno_df <- anno_df%>%mutate(xmin=contrast,xmax=contrast, y1=y+50)
#label <- rep("*", nrow(anno_df))
#label[anno_df$pval<0.01] <- "**"
#label[anno_df$pval<0.001] <- "***"
#anno_df$label <- label
#
#anno_df <- anno_df%>%filter(pval<0.05)
#anno_df$group <- 1:nrow(anno_df)
#
#
##uplab <- c("up", "down")
##names(uplab) <- c("1", "2")
##cols <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
##          "NKcell"="#377eb8", "Tcell"="#e41a1c")
#cols <- c("1"="red","2"="blue")
#mylab <- c("1"="up","2"="down")
#               
#fig0 <- ggplot(sigs,aes(x=contrast, y=ngene, fill=direction))+geom_bar(stat="identity",position="dodge")+  #"stack"
#        scale_fill_manual(values=cols,labels=mylab)+
#        geom_signif(data=anno_df, 
#                    aes(xmin=xmin, xmax=xmax, annotations=label, y_position=y1, group=group),
#                    vjust=0.1, tip_length=0.05, manual=T)+ 
#        facet_wrap(~factor(MCls),ncol=2)+
#        xlab("")+ylab("No. DGE")+ylim(0,900)+
#        theme_bw()+
#        theme(strip.background=element_blank(),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9),
#              strip.text.x=element_text(size=12))
####
#figfn <- "./10_RNA.Variance_output/tmp6/Figure4.4_DGE.barplot2.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig0)
#dev.off()
#}

#########################
### (5). gene example ###
#########################
#if (FALSE){
#
#lab1 <- c("CTRL"="CTRL", 
#          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
#          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#           
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                   sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#gene0 <- "ENSG00000161980"
#symbol0 <- "POLR3K"
#
#cvt2$y <- Vx[grepl(gene0,rownames(Vx)),]
#d2 <- cvt2%>%drop_na(y)
#
#fig2 <- ggplot(d2%>%filter(MCls=="Monocyte"), aes(x=treats, y=log2(y), colour=treats))+
#        geom_boxplot()+ylab(bquote(log[2]~"(variability)"))+
#        scale_colour_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("", labels=lab1)+
#        ggtitle(bquote(~italic(.(symbol0))~" variability in Monocyte cell"))+
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
#              axis.title.y=element_text(size=10),
#              plot.title=element_text(hjust=0.5, size=10),
#              legend.position="none")
#
#figfn <- paste("./10_RNA.Variance_output/tmp6/Figure2.5_", gene0, ".Mono.Boxplot.png", sep="")
#png(figfn, width=400, height=600, res=120)
#print(fig2)  
#dev.off()
#
#}

#####################################################
### (3) plot show number of genes per combination ###
#####################################################
#
#if(FALSE){
#cat("(3).", "Number of effective genes", "\n")
#rm(list=ls())
#
#### (II), Distribution of effective number of genes after removing large outliers
#
#load("./10_RNA.Variance_output/tmp7/1_RNA.Vx.RData")
#load("./10_RNA.Variance_output/tmp7/1_RNA.Sx.mu.RData")
#load("./10_RNA.Variance_output/tmp7/1_RNA.Sx.phi.RData")
#
#Vx[Sx.mu>0.345] <- NA
#Vx[Sx.phi>270.86] <- NA
#xx <- !is.na(Vx)
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], 
#                  treats=gsub("-EtOH", "", cvt0[,2]), sampleID=cvt0[,3], 
#                  Batch2=cvt0[,4], ngene=colSums(xx))  
#                  
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#lab1 <- c("CTRL"="CTRL",
#          "LPS"="LPS", "LPS-DEX"="LPS+DEX", 
#          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#                  
#fig1 <- ggplot(cvt, aes(x=treats, y=ngene, fill=treats))+
#        geom_violin()+ylab("")+
#        facet_wrap(~factor(MCls), nrow=2)+
#        ggtitle("#Genes for each combination")+
#        scale_fill_manual(values=col1)+
#        scale_x_discrete("",labels=lab1)+
#        theme_bw()+
#        theme(legend.position="none",
#              plot.title=element_text(hjust=0.5),
#              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
#                             
#figfn <- "./10_RNA.Variance_output/tmp7/Figure1.3.1.violin.png"
#png(figfn, width=600, height=600, res=120)
#print(fig1)
#dev.off()
#                  
#                  
#fig2 <- ggplot(cvt, aes(x=ngene))+
#     geom_histogram(alpha=0.5,color="grey30")+
#     xlab("#Genes for each combination")+ylab("Counts")+
#     facet_wrap(~factor(MCls), nrow=2, scales="free")+
#     theme_bw()
#     
#figfn <- "./10_RNA.Variance_output/tmp7/Figure1.3.2.hist2.png"
#png(figfn, width=600, height=600, res=120)
#print(fig2)
#dev.off()
### 
#} ### End



#################################################################
#### (4), Box plot, Figures for variance, mean and dispersion ###
#################################################################
#
#if(FALSE){
#
#cat("(4).", "box plot shows distribution of variance, mean and dispersion", "\n")
#
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Bx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Phx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.PhxNew.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.mu.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.phi.RData")
#
#d1 <- melt(Vx)
#d2 <- melt(Bx)
#d3a <- melt(Phx)
#d3b <- melt(PhxNew2)
#d4 <- melt(Sx.mu)
#d5 <- melt(Sx.phi)
#
#dd <- data.frame(X1=d1$X1, X2=d1$X2,
#                 va=d1$value, mu=d2$value, phi=d3a$value, phiNew=d3b$value,
#                 se.mu=d4$value, se.phi=d5$value)
#ddx <- dd%>%
#      drop_na(va, mu, phi, phiNew, se.mu, se.phi)%>%
#      filter(se.mu<0.43, se.phi<546.143)
#
#cvt0 <- str_split(ddx$X2, "_", simplify=T)
#
#dd2 <- ddx%>%
#       mutate(MCls=cvt0[,1],treats=gsub("-EtOH", "", cvt0[,2]))%>%
#       dplyr::rename(ensgene=X1)
#      
#dd2 <- dd2%>%
#       group_by(ensgene, MCls, treats)%>%
#       summarise(va=mean(va), mu=mean(mu), phi=mean(phi), phiNew=mean(phiNew))
#            
#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#            
#col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
#          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
#           
###variance
###va+1e-04
##phi+0.1
### mean value
#fig1 <- ggplot(dd2)+
#        geom_boxplot(aes(x=MCls, y=log2(mu), color=MCls))+
#        scale_color_manual(values=col2)+
#        ylab(bquote(log[2]~"(Mean)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#fig2 <- ggplot(dd2)+
#        geom_boxplot(aes(x=treats, y=log2(mu), color=treats))+
#        scale_color_manual(values=col1)+
#        ylab(bquote(log[2]~"(Mean)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#figfn <- "./10_RNA.Variance_output/tmp6/Figure1.4.3.mu.boxplot.png"
#png(figfn, width=650, height=300, res=110)
#print(plot_grid(fig1, fig2, ncol=2, align="hv",
#                labels="AUTO", label_fontface="plain"))
#dev.off()        
####
#
#### variance facet by MCls and colored by treats
##dd3 <- dd2%>%mutate(x=log2(va+1e-04))
##dd4 <- dd3%>%group_by(treats)%>%nest()
##fig0 <- ggplot(dd2)+
##        geom_boxplot(aes(x=treats, y=log2(va+1e-04), color=treats))+
##        scale_color_manual(values=cols1)+
##        ylab(bquote(log[2]~"(variance)"))+
##        facet_wrap(~MCls,nrow=2,scales="free")+ 
##        theme_bw()+
##        theme(axis.title.x=element_blank(),
##              legend.position="none")
#        
##figfn <- "./10_RNA.Variance_output/tmp5/Figure1.4.3.mu.boxplot1.png"
##png(figfn, width=650, height=600, res=120)
##print(fig0)
##dev.off() 
#
#
#} ###(3), End 



######################################
#### 3, scatter plot specific gene ###
######################################
#### fig3, scatter plots
#rm(list=ls())
#load("./10_RNA.Variance2_output/10.RNA.Vx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Bx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Phx.RData")
#
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#gene0 <- "ENSG00000169442"
#x <- Bx[gene0,]
#y <- Vx[gene0,]
#cvt$x <- x
#cvt$y <- y
#subi <- (!is.na(x))&(!is.na(y))
##dd <- data.frame(x=x[subi],y=y[subi])
#
#cols1 <- c("CTRL"="#828282", 
#           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4") 
#dd <- cvt %>%filter((!is.na(x))&(!is.na(y)))
####
#fig1 <- ggplot(dd, aes(x=x, y=y))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab("mean value")+
#        ylab("variance")+
#        geom_smooth(method="lm",formula=y~poly(x,2))+
#        #facet_wrap(~MCls,nrow=2)+
#        theme_bw()+theme(legend.position="none")
#        
###        
#dd0 <- dd%>%filter(x>0.001)
#fig2 <- ggplot(dd0,aes(x=log10(x), y=log10(y+0.001)))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        #xlab("mean value")+
#        #ylab("variance")+
#        #geom_smooth(method="lm",formula=y~poly(x,2))+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"(va"~+~0.001~")"))+
#        #facet_wrap(~MCls,nrow=2)+
#        theme_bw()+
#        theme(legend.position=c(0.8,0.3),legend.background=element_blank())
#                
#figfn <- paste("./10_RNA.Variance2_output/Figure0.3.va.", gene0, ".scatter.pdf", sep="")
#pdf(figfn,width=10,height=6)
#plot_grid(fig1, fig2, ncol=2)
#dev.off()   
#
#
#### fig4, scatter plots for dispersion
#gene0 <- "ENSG00000169442"
#x <- Bx[gene0,]
#y <- Phx[gene0,]
#cvt$x <- x
#cvt$y <- y
##subi <- (!is.na(x))&(!is.na(y))
##dd <- data.frame(x=x[subi],y=y[subi])
#dd <- cvt %>%filter((!is.na(x))&(!is.na(y)))
#fig1 <- ggplot(dd,aes(x=x, y=y))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab("mean value")+ylab("dispersion")+
#        theme_bw()+theme(legend.position="none")
#
###        
#dd0 <- dd%>%filter(x>0.001)
#fig2 <- ggplot(dd0,aes(x=log10(x), y=log10(y+0.001)))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"("~phi~+~0.001~")"))+
#        theme_bw()+
#        theme(legend.position=c(0.8,0.3),legend.background=element_blank())
#        
#        
#        
#figfn <- paste("./10_RNA.Variance2_output/Figure0.3.phi.", gene0, ".scatter.pdf", sep="")
#pdf(figfn,width=10,height=6)
#plot_grid(fig1, fig2, ncol=2)
#dev.off()  
#} ###End
# 





###########################################################################################
##### 7, scatter plot for beta between LPS vs LPS-DEX, PHA vs PHA-DEX in va, phi and mu ###
###########################################################################################
#
#### 
#if(FALSE){
#
#cat("6.3", "\n")
#rm(list=ls())
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#
#      
#### label function       
#feq <- function(x){
#  r <- format(as.numeric(x$estimate),digits=3)
#  p <- x$p.value
#  if(p<0.001) symb <- "***"
#  if(p>=0.001 & p<0.01) symb <- "**"
#  if (p>=0.01 & p<0.05) symb <- "*"
#  if(p>0.05) symb <- "NS"
#  
#  eq <- bquote(italic(R)==.(r)~","~.(symb))
#  eq 
#}       
#
#### (1), variance
#cat("(1)", "compare beta of variance between no DEX and DEX", "\n")
###va data frame
#dfr <- map_dfr(MCls, function(x){
#   fn <- paste("./10_RNA.Variance_output/tmp5/2_Batch_va.", x, ".meta", sep="")
#   res <- read.table(file=fn,header=T)%>%mutate(MCls=x)  
#})
#dfr <- dfr%>%mutate(rn=paste(MCls, gene, sep="_"))
#
#             
#### (2), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH
#dfa <- dfr%>%filter(contrast=="LPS")    
#dfb <- dfr%>%filter(contrast=="LPS-DEX")%>%
#       dplyr::select(rn, beta, pval, qval)
#       
#df1 <- dfa%>%inner_join(dfb,by="rn")
#
#anno_df1 <- df1%>%group_by(MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#           #unnest()
#     
#fig1 <- ggplot(df1, aes(x=beta.x,y=beta.y))+
#        geom_point(size=0.3,color="grey30")+ 
#        geom_text(data=anno_df1, x=4.8, y=7.5, aes(label=eq), colour="blue", size=2.5, parse=T)+ 
#        facet_wrap(~MCls, nrow=2)+
#        xlab(bquote(Beta~"from LPS-EtOH vs CTRL"))+ylab(bquote(Beta~"from LPS-DEX vs LPS-EtOH"))+
#        theme_bw()+
#        theme(strip.background=element_blank())
#                           
#figfn <- "./10_RNA.Variance_output/tmp5/Figure6.1a.va.png"
#png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
#print(fig1)
#dev.off()
#
#
#### (3), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
#
#dfa <- dfr%>%filter(contrast=="PHA")    
#dfb <- dfr%>%filter(contrast=="PHA-DEX")%>%
#       dplyr::select(rn, beta, pval, qval)
#       
#df2 <- dfa%>%inner_join(dfb,by="rn")
#
#anno_df2 <- df2%>%group_by(MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#           #unnest()
#     
#fig2 <- ggplot(df2, aes(x=beta.x,y=beta.y))+
#        geom_point(size=0.3,color="grey30")+ 
#        geom_text(data=anno_df2, x=5, y=7.5, aes(label=eq), colour="blue", size=2.5, parse=T)+ 
#        facet_wrap(~MCls, nrow=2)+
#        xlab(bquote(Beta~"from PHA-EtOH vs CTRL"))+ylab(bquote(Beta~"from PHA-DEX vs LPS-EtOH"))+
#        theme_bw()+
#        theme(strip.background=element_blank())
#                           
#figfn <- "./10_RNA.Variance_output/tmp5/Figure6.1b.va.png"
#png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
#print(fig2)
#dev.off()
#
#
#
#####
#### (3), mean expression
#cat("(3)", "compare beta of mean expression between no DEX and DEX", "\n")
###mean(mu), data frame
#df3 <- map_dfr(MCls, function(x){
#   fn <- paste("./10_RNA.Variance_output/tmp4/4_Batch_mu.", x, ".meta", sep="")
#   res <- read.table(file=fn,header=T)%>%mutate(MCls=x)  
#})
#df3 <- df3%>%
#       dplyr::select(beta, gene, contrast, qval, MCls)%>%
#       mutate(rn=paste(MCls, contrast, gene, sep="_"))
#
### LPS
#dfa <- df3%>%
#       filter(contrast=="LPS")%>%
#       mutate(rn1=paste(MCls, gene, sep="_"))
#dfb <- df3%>%
#       filter(contrast=="LPS-DEX")%>%
#       mutate(rn1=paste(MCls, gene, sep="_"))%>%
#       dplyr::select(beta, contrast, rn1)
#       
#dfcomb1 <- dfa%>%inner_join(dfb,by="rn1")%>%mutate(compare="LPS vs LPS-DEX")
#
### PHA
#dfa <- df3%>%
#       filter(contrast=="PHA")%>%
#       mutate(rn1=paste(MCls, gene, sep="_"))
#dfb <- df3%>%
#       filter(contrast=="PHA-DEX")%>%
#       mutate(rn1=paste(MCls, gene, sep="_"))%>%
#       dplyr::select(beta, contrast, rn1)
#dfcomb2 <- dfa%>%inner_join(dfb,by="rn1")%>%mutate(compare="PHA vs PHA-DEX")
#
#dfcomb <- rbind(dfcomb1, dfcomb2)
#
#anno_df3 <- dfcomb%>%
#           group_by(compare, MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#           #unnest()
#
#fig3 <- ggplot(dfcomb, aes(x=beta.x,y=beta.y))+
#        geom_point(size=0.3)+
#        geom_text(data=anno_df3, x=0.7, y=0.85, aes(label=eq), colour="blue", size=2.5, parse=T)+
#        facet_grid(compare~MCls)+
#        xlab(bquote(beta~"no DEX"))+ylab(bquote(beta~"with DEX"))+
#        theme_bw()+
#        theme(#strip.background=element_blank(),
#              #strip.text.x=element_text(size=12),
#              legend.title=element_blank(),
#              legend.text=element_text(color="black",size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text.x=element_text(color="black",size=9),
#              axis.text.y=element_text(color="black",size=9))
#
#figfn <- "./10_RNA.Variance_output/tmp4/Figure6.3.mu.png"
#png(filename=figfn, width=900, height=500, pointsize=12, res=130)  
#print(fig3)
#dev.off()
#
#} ### 6.1, 6.2 and 6.3, End
#        
#       
       
              
       
##############################################################
### Part 2, Show Figures for variance, mean and dispersion ###
##############################################################


###############################################################
#### 1, Box plot, Figures for variance, mean and dispersion ###
###############################################################

#if(FALSE){

#rm(list=ls())
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Bx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Phx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.PhxNew.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.mu.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.phi.RData")

#d1 <- melt(Vx)
#d2 <- melt(Bx)
#d3a <- melt(Phx)
#d3b <- melt(PhxNew2)
#d4 <- melt(Sx.mu)
#d5 <- melt(Sx.phi)

#dd <- data.frame(X1=d1$X1, X2=d1$X2,
#                 va=d1$value, mu=d2$value, phi=d3a$value, phiNew=d3b$value,
#                 se.mu=d4$value, se.phi=d5$value)
#ddx <- dd%>%
#      drop_na(va, mu, phi, phiNew, se.mu, se.phi)%>%
#      filter(se.mu<0.43, se.phi<546.143)

# <- cvt0 <- str_split(ddx$X2, "_", simplify=T)

#dd2 <- ddx%>%
#       mutate(MCls=cvt0[,1],treats=gsub("-EtOH", "", cvt0[,2]))%>%
#       dplyr::rename(ensgene=X1)
      
#dd2 <- dd2%>%
#       group_by(ensgene, MCls, treats)%>%
#       summarise(va=mean(va), mu=mean(mu), phi=mean(phi),phiNew=mean(phiNew))
#            
#cols1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4") 
#cols2 <- c("Tcell"="#e41a1c", "NKcell"="#377eb8", 
#          "Bcell"="#4daf4a", "Monocyte"="#984ea3")
           
#fig1 <- ggplot(dd2,aes(x=factor(treats),y=log10(va+1e-04),fill=factor(treats)))+
#        geom_boxplot()+
#        scale_fill_manual("", values=cols1)+
#        facet_wrap(~MCls, nrow=2)+
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              axis.text.x=element_text(angle=90),
#              legend.position="none")
#figfn <- "./10_RNA.Variance_output/tmp5/Figure1.4.1.va.boxplot2.png"
#png(figfn, width=800, height=800, res=120)
#print(fig1)
#dev.off()
##mean value            
#fig2 <- ggplot(dd2)+
#        geom_boxplot(aes(x=factor(treats),y=log10(mu),fill=factor(treats)))+
#        scale_fill_manual("", values=cols1)+
#        facet_wrap(~MCls, nrow=2)+
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              axis.text.x=element_text(angle=90),
#              legend.position="none")
#figfn <- "./10_RNA.Variance_output/tmp5/Figure1.4.2.mu.boxplot2.png"
#png(figfn, width=800, height=800, res=120)
#print(fig2)
#dev.off()

##dispersion
#fig3 <- ggplot(dd2)+
#        geom_boxplot(aes(x=factor(treats),y=log10(phi+0.1),fill=factor(treats)))+
#        scale_fill_manual("", values=cols1)+
#        facet_wrap(~MCls, nrow=2)+
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              axis.text.x=element_text(angle=90),
#              legend.position="none")
#figfn <- "./10_RNA.Variance_output/tmp5/Figure1.4.3.phi.boxplot2.png"
#png(figfn, width=800, height=800, res=120)
#print(fig3)
#dev.off()

##variance
##va+1e-04
#phi+0.1
#fig1 <- ggplot(dd2)+
#        geom_boxplot(aes(x=MCls, y=log2(mu), color=MCls))+
#        scale_color_manual(values=cols2)+
#        ylab(bquote(log[2]~"(Mean)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#fig2 <- ggplot(dd2)+
#        geom_boxplot(aes(x=treats, y=log2(mu), color=treats))+
#        scale_color_manual(values=cols1)+
#        ylab(bquote(log[2]~"(Mean)"))+ 
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              legend.position="none")
#        
#figfn <- "./10_RNA.Variance_output/tmp6/Figure1.4.3.mu.boxplot.png"
#png(figfn, width=650, height=300, res=110)
#print(plot_grid(fig1, fig2, ncol=2, align="hv",
#                labels="AUTO", label_fontface="plain"))
#dev.off()        
####
#
#### variance facet by MCls and colored by treats
##dd3 <- dd2%>%mutate(x=log2(va+1e-04))
##dd4 <- dd3%>%group_by(treats)%>%nest()
##fig0 <- ggplot(dd2)+
##        geom_boxplot(aes(x=treats, y=log2(va+1e-04), color=treats))+
##        scale_color_manual(values=cols1)+
##        ylab(bquote(log[2]~"(variance)"))+
##        facet_wrap(~MCls,nrow=2,scales="free")+ 
##        theme_bw()+
##        theme(axis.title.x=element_blank(),
##              legend.position="none")
#        
##figfn <- "./10_RNA.Variance_output/tmp5/Figure1.4.3.mu.boxplot1.png"
##png(figfn, width=650, height=600, res=120)
##print(fig0)
##dev.off() 
#
#
#} ###(3), End 


############################################################################################
### 2, scatter plots, showing relations between mean variance, mean dispersion ###
############################################################################################  
#if (FALSE){
#
#rm(list=ls())
##MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
##Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Vx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Bx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Phx.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.PhxNew.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.mu.RData")
#load("./10_RNA.Variance_output/tmp6/1_RNA.Sx.phi.RData")
#
#d1 <- melt(Vx)
#d2 <- melt(Bx)
#d3a <- melt(Phx)
#d3b <- melt(PhxNew2)
#d4 <- melt(Sx.mu)
#d5 <- melt(Sx.phi)
#
#dd <- data.frame(X1=d1$X1, X2=d1$X2,
#                 va=d1$value, mu=d2$value, phi=d3a$value, phiNew=d3b$value,
#                 se.mu=d4$value, se.phi=d5$value)
#ddx <- dd%>%
#      drop_na(va, mu, phi, phiNew, se.mu, se.phi)%>%
#      filter(se.mu<0.43, se.phi<546.143)
#
#cvt0 <- str_split(ddx$X2, "_", simplify=T)
#
#dd2 <- ddx%>%
#       mutate(MCls=cvt0[,1],treats=gsub("-EtOH", "", cvt0[,2]))%>%
#       dplyr::rename(ensgene=X1)
#      
##dd2 <- dd2%>%
##       group_by(ensgene, MCls, treats)%>%
##       summarise(va=mean(va), mu=mean(mu), phi=mean(phi),phiNew=mean(phiNew))

 


#dd2 <- dd%>%drop_na(va, mu, phi)%>%
#       group_by(ensgene, MCls, treats)%>%
#       summarize(va.mean=mean(va), mu.mean=mean(mu), phi.mean=mean(phi))

#ddnew <- ddnew%>%drop_na(va,mu,phi)%>%filter(va<1.1e+05,mu>0.001&mu<753) 
#cols1 <- c("CTRL"="#828282", 
#           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4") 


#############################
### option-1.5, log10     ###
#############################                    
###figure 1 
#dd3 <- dd2%>%mutate(x=log10(mu), y=log10(va+1e-04))  
#
#anno_df1 <- dd3%>%
#            group_by(treats, MCls)%>%
#            nest()%>%
#            mutate(corr=map_dbl(data, ~cor((.x)$x, (.x)$y, method="pearson")))
#
#fig1 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log["10"]~"("~mu~")"))+
#        ylab(bquote(log["10"]~"(va"~+~1e-04~")"))+          
#        #facet_wrap(~MCls, nrow=2, scales="free")+
#        facet_grid(treats~MCls)+
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()
                 
#fig1 <- ggplot(dd3, aes(x=x, y=y))+
#        geom_point(size=0.05,alpha=0.01)+
        #scale_color_manual("",values=cols1, guide=guide_legend(override.aes=list(size=3)))+
#        xlab(bquote(log["10"]~"("~mu~")"))+
#        ylab(bquote(log["10"]~"(va"~+~1e-08~")"))+          
#        facet_wrap(~MCls,nrow=2)+    
#        geom_smooth(method="loess",formula=y~poly(x,2),span=0.3)+

#figfn <- "./10_RNA.Variance_output/tmp6/Figure1.5.1.va_mu.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig1)
#dev.off() 
  

#### figure 2, scatter plots for dispersion
#dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phiNew))  ##phi+0.01
#fmod <- function(df){
#   lm0 <- lm(y~x, data=df)
#   summary(lm0)$r.squared
#}
#anno_df2 <- dd3%>%
#            group_by(treats, MCls)%>%
#            nest()%>%
#            mutate(corr=map_dbl(data, ~cor((.x)$x, (.x)$y, method="pearson")),
#                   r2=map_dbl(data, fmod))
#fig2 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log["10"]~"("~mu~")"))+
#        ylab(bquote(log["10"]~"("~phi~")"))+
#        ylim(-2,1.5)+
#        facet_grid(treats~MCls)+          
#        #facet_wrap(~MCls, nrow=2, scales="free")+
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()
#        
##fig2 <- ggplot(dd2%>%filter(mu.mean<2.6e+06),aes(x=log10(mu.mean), y=log10(phi.mean)))+
##        geom_point(size=0.05)+
##        #scale_color_manual("",values=cols1)+
##        xlab(bquote(log[10]~"("~mu~")"))+
##        ylab(bquote(log[10]~"("~phi~+~0.001~")"))+
##        geom_smooth(method="lm",formula=y~x)+
##        facet_wrap(~MCls,nrow=2)+
##        theme_bw()+theme(legend.position="none")
#figfn <- "./10_RNA.Variance_output/tmp6/Figure1.5.2.phi_mu.scatter2.png"
#png(figfn, width=900, height=800, res=150)
#print(fig2)
#dev.off()
#
#
#### label function       
#feq <- function(x){
#  r <- format(as.numeric(x$estimate),digits=3)
#  p <- x$p.value
#  if(p<0.001) symb <- "***"
#  if(p>=0.001 & p<0.01) symb <- "**"
#  if (p>=0.01 & p<0.05) symb <- "*"
#  if(p>0.05) symb <- "NS"
#  
#  eq <- bquote(italic(R)==.(r)~","~.(symb))
#  as.character(as.expression(eq)) 
#} 
#
#############################
#### option 1.6, log      ###
#############################
#### (1)
#if(FALSE){
#cat("1.6.1", "va vs mu", "\n")
#dd3 <- dd2%>%mutate(x=log10(mu),y=log10(va+1e-06))
#corr <- cor.test(dd3$x,dd3$y)
#fig1 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        geom_text(x=-1.2, y=0, label=feq(corr), parse=T, size=3.5)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"(va"~+~1e-06~")"))+          
#        #facet_wrap(~MCls,nrow=2,scales="free")+
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()
#figfn <- "./10_RNA.Variance_output/tmp4/Figure1.6.1.va.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig1)
#dev.off()
#
#### (2)
#cat("1.6.2", "mu vs phi", "\n")
#dd3 <- dd2%>%mutate(x=log10(mu),y=log10(phi+1e-04))
#corr <- cor.test(dd3$x, dd3$y)
#fig2 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        geom_text(x=0.2, y=0, label=feq(corr), parse=T, size=3.5)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log~"("~mu~")"))+
#        ylab(bquote(log[10]~"("~phi~+~1e-04~")"))+         
#       #facet_wrap(~MCls,nrow=2,scales="free")+
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()
#figfn <- "./10_RNA.Variance_output/tmp4/Figure1.6.2.phi.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig2)
#dev.off()
#
#### (3)
#cat("1.6.3", "va vs phi")
#dd3 <- dd2%>%mutate(x=log10(va+1e-06),y=log10(phi+1e-04))
#corr <- cor.test(dd3$x, dd3$y)
#fig3 <- ggplot(dd3, aes(x=x,y=y))+
#        stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+
#        geom_text(x=0, y=0, label=feq(corr), parse=T, size=3.5)+
#        scale_fill_viridis_c()+
#        xlab(bquote(log~"(va"~+~1e-06~")"))+
#        ylab(bquote(log[10]~"("~phi~+~1e-04~")"))+          
#       #facet_wrap(~MCls,nrow=2,scales="free")+
#        #geom_smooth(method="lm",formula=y~x)+
#        theme_bw()
#figfn <- "./10_RNA.Variance_output/tmp4/Figure1.6.3.va_phi.scatter.png"
#png(figfn, width=900, height=800, res=150)
#print(fig3)
#dev.off()
#} ###  End, 1.6
#
#} ###End   
                           #

##########################################
#### 3.7, show plots for specific gene ###
##########################################
#if(FALSE){
#rm(list=ls())
#load("./10_RNA.Variance2_output/10.RNA.Vx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Bx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Phx.RData")
#
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
#
##Naii <- is.na(Vx)
##Num <- rowSums(Naii)
##Vtmp <- Vx[Num==0,]
##range(Vtmp)
##cols <- pal_npg("nrc",alpha=0.6)(5)
#cols1 <- c("CTRL"="#828282", 
#           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")   
#
### fig2, hist plot for variance
#gene0 <- "ENSG00000169442"
#cvt$v <- Vx[gene0,]
#dd1 <- cvt %>% filter(!is.na(v))
#fig1 <- ggplot(dd1,aes(x=log(v)))+
#        geom_histogram(color="#e9ecef", alpha=0.6, position="identity")+
#        xlab("log(variance)")+facet_grid(MCls~treats)+theme_bw()
####
#figfn <- paste("./10_RNA.Variance_output/Figure2.1_va.", gene0, ".hist.log.pdf", sep="")
#pdf(figfn)
#fig1
#dev.off()
#
#### fig2. box plot for variance
#gene0 <- "ENSG00000169442"
#cvt$v <- Vx[gene0,]
#dd1 <- cvt %>% filter(!is.na(v))
#fig2 <- ggplot(dd1, aes(x=factor(treats),y=v, fill=factor(treats)))+
#        geom_boxplot()+ylab("Variance")+
#        scale_fill_manual("", values=cols1)+
#        facet_wrap(~MCls,nrow=2)+
#        theme_bw()+
#        theme(axis.title.x=element_blank(),
#              axis.text.x=element_text(angle=90),
#              legend.position="none")
#              #plot.title = element_text(hjust = 0.5))
#              
#figfn <- paste("./10_RNA.Variance_output/Figure2.2_va.", gene0, ".box.pdf", sep="")
#pdf(figfn)
#fig2
#dev.off()
#
#
#
######################################
#### 3, scatter plot specific gene ###
######################################
#### fig3, scatter plots
#rm(list=ls())
#load("./10_RNA.Variance2_output/10.RNA.Vx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Bx.RData")
#load("./10_RNA.Variance2_output/10.RNA.Phx.RData")
#
#bti2 <- colnames(Vx)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treats=cvt0[,2], sampleID=cvt0[,3], Batch2=cvt0[,4])
#
#gene0 <- "ENSG00000169442"
#x <- Bx[gene0,]
#y <- Vx[gene0,]
#cvt$x <- x
#cvt$y <- y
#subi <- (!is.na(x))&(!is.na(y))
##dd <- data.frame(x=x[subi],y=y[subi])
#
#cols1 <- c("CTRL"="#828282", 
#           "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4") 
#dd <- cvt %>%filter((!is.na(x))&(!is.na(y)))
####
#fig1 <- ggplot(dd, aes(x=x, y=y))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab("mean value")+
#        ylab("variance")+
#        geom_smooth(method="lm",formula=y~poly(x,2))+
#        #facet_wrap(~MCls,nrow=2)+
#        theme_bw()+theme(legend.position="none")
#        
###        
#dd0 <- dd%>%filter(x>0.001)
#fig2 <- ggplot(dd0,aes(x=log10(x), y=log10(y+0.001)))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        #xlab("mean value")+
#        #ylab("variance")+
#        #geom_smooth(method="lm",formula=y~poly(x,2))+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"(va"~+~0.001~")"))+
#        #facet_wrap(~MCls,nrow=2)+
#        theme_bw()+
#        theme(legend.position=c(0.8,0.3),legend.background=element_blank())
#                
#figfn <- paste("./10_RNA.Variance2_output/Figure0.3.va.", gene0, ".scatter.pdf", sep="")
#pdf(figfn,width=10,height=6)
#plot_grid(fig1, fig2, ncol=2)
#dev.off()   
#
#
#### fig4, scatter plots for dispersion
#gene0 <- "ENSG00000169442"
#x <- Bx[gene0,]
#y <- Phx[gene0,]
#cvt$x <- x
#cvt$y <- y
##subi <- (!is.na(x))&(!is.na(y))
##dd <- data.frame(x=x[subi],y=y[subi])
#dd <- cvt %>%filter((!is.na(x))&(!is.na(y)))
#fig1 <- ggplot(dd,aes(x=x, y=y))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab("mean value")+ylab("dispersion")+
#        theme_bw()+theme(legend.position="none")
#
###        
#dd0 <- dd%>%filter(x>0.001)
#fig2 <- ggplot(dd0,aes(x=log10(x), y=log10(y+0.001)))+
#        geom_point(aes(colour=factor(treats)))+scale_color_manual("",values=cols1)+
#        xlab(bquote(log[10]~"("~mu~")"))+
#        ylab(bquote(log[10]~"("~phi~+~0.001~")"))+
#        theme_bw()+
#        theme(legend.position=c(0.8,0.3),legend.background=element_blank())
#        
#        
#        
#figfn <- paste("./10_RNA.Variance2_output/Figure0.3.phi.", gene0, ".scatter.pdf", sep="")
#pdf(figfn,width=10,height=6)
#plot_grid(fig1, fig2, ncol=2)
#dev.off()  
#} ###End
# 


#########################################
#### 4.6, scatter plot of LPS and PHA ###
#########################################
#
#if(FALSE){
#cat("4.6", "\n")
#rm(list=ls())
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#
#      
#### label function       
#feq <- function(x){
#  r <- format(as.numeric(x$estimate),digits=3)
#  p <- x$p.value
#  if(p<0.001) symb <- "***"
#  if(p>=0.001 & p<0.01) symb <- "**"
#  if (p>=0.01 & p<0.05) symb <- "*"
#  if(p>0.05) symb <- "NS"
#  
#  eq <- bquote(italic(R)==.(r)~","~.(symb))
#  eq 
#}       
#
#### (1), residual dispersion 
#cat("(1)", "compare beta of variance between no DEX and DEX", "\n")
###va data frame
#dfr <- map_dfr(MCls, function(x){
#   fn <- paste("./10_RNA.Variance_output/tmp6/3_phiNew.", x, ".meta", sep="")
#   res <- read.table(file=fn,header=T)%>%mutate(MCls=x)  
#})
#dfr <- dfr%>%mutate(rn=paste(MCls, gene, sep="_"))
#
#             
#### (2), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH
#dfa <- dfr%>%filter(contrast=="LPS")    
#dfb <- dfr%>%filter(contrast=="LPS-DEX")%>%
#       dplyr::select(rn, beta, pval, qval)
#       
#df1 <- dfa%>%inner_join(dfb,by="rn")
#
#anno_df1 <- df1%>%group_by(MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#           #unnest()
#     
#fig1 <- ggplot(df1, aes(x=beta.x,y=beta.y))+
#        geom_point(size=0.3,color="grey30")+ 
#        geom_text(data=anno_df1, x=4.8, y=7.5, aes(label=eq), colour="blue", size=2.5, parse=T)+ 
#        facet_wrap(~MCls, nrow=2)+
#        xlab(bquote(Beta~"from LPS-EtOH vs CTRL"))+ylab(bquote(Beta~"from LPS-DEX vs LPS-EtOH"))+
#        theme_bw()+
#        theme(strip.background=element_blank())
#                           
#figfn <- "./10_RNA.Variance_output/tmp6/Figure3x.5_LPS_phi.png"
#png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
#print(fig1)
#dev.off()
#
#
#### (3), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
#
#dfa <- dfr%>%filter(contrast=="PHA")    
#dfb <- dfr%>%filter(contrast=="PHA-DEX")%>%
#       dplyr::select(rn, beta, pval, qval)
#       
#df2 <- dfa%>%inner_join(dfb,by="rn")
#
#anno_df2 <- df2%>%group_by(MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#           #unnest()
#     
#fig2 <- ggplot(df2, aes(x=beta.x,y=beta.y))+
#        geom_point(size=0.3,color="grey30")+ 
#        geom_text(data=anno_df2, x=5, y=7.5, aes(label=eq), colour="blue", size=2.5, parse=T)+ 
#        facet_wrap(~MCls, nrow=2)+
#        xlab(bquote(Beta~"from PHA-EtOH vs CTRL"))+ylab(bquote(Beta~"from PHA-DEX vs LPS-EtOH"))+
#        theme_bw()+
#        theme(strip.background=element_blank())
#                           
#figfn <- "./10_RNA.Variance_output/tmp6/Figure3x.5_PHA.phi.png"
#png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
#print(fig2)
#dev.off()
#} ###

###

###
### overlap 
#if(FALSE){
#cat("6.2", "scatter plot to compare beta between expression, va, phi and mu", "\n")
#
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#Newcon2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
#           "PHA"="PHA", "PHA-DEX"="PHA+DEX")
#
#
##### (1). pseudo-bulk differential 
##df1 <- read_rds("./6_DEG.CelltypeNew_output/2_meta.rds")%>%
##       drop_na(beta, qval)%>%
##       mutate(zscore=beta/stderr, rn=paste(MCls, contrast, gene, sep="_"))%>%
##       arrange(desc(abs(zscore)))
##
##
##### (2). variance
##fn <- "./10_RNA.Variance_output/tmp7/2_va2.meta"
##df2 <- read.table(file=fn,header=T)%>%mutate(zscore=beta/stderr)%>%drop_na(beta,qval)  
##df2 <- df2%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))
##df2 <- df2%>%arrange(desc(abs(zscore)))
##
##
##### (3). dispersion
##fn <- "./10_RNA.Variance_output/tmp7/3_phi.meta"
##df3 <- read.table(file=fn,header=T)%>%mutate(zscore=beta/stderr)%>%drop_na(beta,qval)  
##df3 <- df3%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))
##df3 <- df3%>%arrange(desc(abs(zscore)))
#
#
#### (4). residual dispersion
#fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
#df4 <- read.table(file=fn,header=T)%>%
#       mutate(zscore=beta/stderr,rn=paste(MCls, contrast, gene, sep="_"))%>%drop_na(beta,qval) 
#       
#
#### (5). mean expression
#fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
#df5 <- read.table(file=fn, header=T)%>%
#       mutate(zscore=beta/stderr, rn=paste(MCls, contrast, gene, sep="_"))%>%
#       drop_na(beta,qval)  
#     
#
#
#
#### my fun 1, generate data frame used for plots 
#myDFxy <- function(dfx, dfy){
#
#   dfx <- dfx%>%dplyr::select(zscore, beta, qval, rn, contrast, MCls)
#   dfy <- dfy%>%dplyr::select(zscore, beta, qval, rn)       
#   dfxy <- dfx%>%
#           inner_join(dfy, by="rn")
#        
#   x <- dfxy$qval.x
#   y <- dfxy$qval.y
#   Bx <- abs(dfxy$beta.x)
#   By <- abs(dfxy$beta.y)
#   gr <- rep(1, nrow(dfxy))
#   gr[x<0.1&Bx>0.5] <- 2
#   gr[y<0.1&By>0.5] <- 3
#   gr[(x<0.1&Bx>0.5)&(y<0.1&By>0.5)] <- 4
#   dfxy$gr <- gr
####
#   dfxy
#}
#
#### label expression function       
#feq <- function(x){
#  #r <- format(as.numeric(x$estimate),digits=1)
#  r <- round(as.numeric(x$estimate), digits=3)
#  p <- x$p.value
#  if(p<0.001) symb <- "***"
#  if(p>=0.001 & p<0.01) symb <- "**"
#  if (p>=0.01 & p<0.05) symb <- "*"
#  if(p>0.05) symb <- "NS"
#  
#  eq <- bquote(italic(R)==.(r)~","~.(symb))
#  eq
#  #r 
#}
#
#mycol <- c("1"="grey50", "2"="red", "3"="blue", "4"="black")
#
##### (1). bulk expression vs variance
##cat("(1).", "expression vs variance", "\n")
##dfxy1 <- myDFxy(df1, df2)
##
##
##mylabel <- c("1"="NS", "2"="DEG", "3"="DVG(variance)", "4"="Both")
##
###green for DEG, blue for DVG, and Black for Both.
##
##anno_df1 <- dfxy1%>%
##           group_by(contrast, MCls)%>%
##           nest()%>%
##           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##                  eq=map(corr,feq))%>%
##           dplyr::select(-data,-corr)#%>%
##        
##fig1 <- ggplot(dfxy1, aes(x=zscore.x, y=zscore.y))+
##        geom_point(aes(colour=factor(gr)), size=0.3)+ #ylim(-12,12)+
##        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
##        geom_text(data=anno_df1, x=-1.3, y=30, aes(label=eq), size=2.5, parse=T)+ 
##        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
##        xlab("z score of gene expression")+
##        ylab("z score of gene variance")+
##        theme_bw()+
##        theme(#strip.background=element_blank(),
##              #strip.text.x=element_text(size=12),
##              legend.title=element_blank(),
##              legend.text=element_text(size=9),
##              #legend.key.size=unit(0.4,units="cm"),
##              axis.text=element_text(size=9))
##
#####              
##figfn <- "./10_RNA.Variance_output/tmp7/Figure5.1_DEGvsDGP.png"
##png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
##print(fig1)
##dev.off() 
##       
##### (2). bulk expression vs dispersion
##cat("(2).", "expression vs dispersion", "\n")
##
##dfxy2 <- myDFxy(df1, df3)
##
##mycol <- c("1"="grey50", "2"="green", "3"="blue", "4"="red")
##mylabel <- c("1"="NS", "2"="DEG", "3"="DVG(dispersion)", "4"="Both")
##
##anno_df2 <- dfxy2%>%
##           group_by(contrast, MCls)%>%
##           nest()%>%
##           mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##                  eq=map(corr,feq))%>%
##           dplyr::select(-data,-corr)#%>%
##        
##fig2 <- ggplot(dfxy2, aes(x=beta.x, y=beta.y))+
##        geom_point(aes(colour=factor(gr)), size=0.3)+
##        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
##        geom_text(data=anno_df2, x=-0.5, y=10, aes(label=eq), size=2.5, parse=T)+ 
##        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
##        #xlab("zscore of gene expression")+ylab("zscore of gene dispersion")+
##        xlab(bquote(beta~"of gene expression"))+
##        ylab(bquote(beta~"of gene dispersion"))+
##        theme_bw()+
##        theme(#strip.background=element_blank(),
##              #strip.text.x=element_text(size=12),
##              legend.title=element_blank(),
##              legend.text=element_text(size=9),
##              #legend.key.size=unit(0.4,units="cm"),
##              axis.text=element_text(size=9))
##
#####              
##figfn <- "./10_RNA.Variance_output/tmp7/Figure7.2_DEGvsDGP.png"
##png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
##print(fig2)
##dev.off()
##
##
##### (3). bulk expression vs residual dispersion 
#####
##cat("(3).", "expression vs residual disperison", "\n")
##dfxy3 <- myDFxy(df5, df4)
##  
##mycol <- c("1"="grey50", "2"="green", "3"="blue", "4"="red")
##mylabel <- c("1"="NS", "2"="DMG(only)", "3"="DVG(only)", "4"="Both")
##
##anno_df3 <- dfxy3%>%
##           group_by(contrast, MCls)%>%
##           nest()%>%
##           mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
##                  eq=map(corr,feq))%>%
##           dplyr::select(-data,-corr)#%>%
##        
##fig3 <- ggplot(dfxy3, aes(x=zscore.x, y=zscore.y))+
##        geom_point(aes(colour=factor(gr)), size=0.3)+
##        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
##        geom_text(data=anno_df3, x=-0.5, y=-10, aes(label=eq), size=2.5, parse=T)+  #-0.5, 28
##        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
##        xlab("zscore of gene mean")+ylab("zscore of gene dispersion")+
##        #xlab(bquote(beta~"of gene expression"))+
##        #ylab(bquote(beta~"of gene dispersion"))+
##        theme_bw()+
##        theme(#strip.background=element_blank(),
##              #strip.text.x=element_text(size=12),
##              legend.title=element_blank(),
##              legend.text=element_text(size=9),
##              #legend.key.size=unit(0.4,units="cm"),
##              axis.text=element_text(size=9))
##
#####              
##figfn <- "./10_RNA.Variance_output/tmp7/tmp1_Monocyte/Figure8.1_DGPvsDMG.png"
###png(filename=figfn, width=900, height=800, pointsize=12, res=130) 
##png(filename=figfn, width=600, height=800, pointsize=12, res=130) 
##print(fig3)
##dev.off()
##
##
##### (4). bulk expression vs mean value
##cat("(4).", "expression vs mean value", "\n")
##dfxy4 <- myDFxy(df1, df5)
##
##mycol <- c("1"="grey50", "2"="green", "3"="blue", "4"="red")
##mylabel <- c("1"="NS", "2"="DEG(only)", "3"="DMG(only)", "4"="Both")
##
##anno_df4 <- dfxy4%>%
##           group_by(contrast, MCls)%>%
##           nest()%>%
##           mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
##                  eq=map(corr,feq))%>%
##           dplyr::select(-data,-corr)#%>%
##        
##fig4 <- ggplot(dfxy4, aes(x=zscore.x, y=zscore.y))+
##        geom_point(aes(colour=factor(gr)), size=0.3)+
##        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
##        geom_text(data=anno_df4, x=0, y=55, aes(label=eq), size=2.5, parse=T)+
##        xlab("zscore of gene expression")+ylab("zscore of gene mean")+
##        #geom_text(data=anno_df4, x=-1.3, y=3, aes(label=eq), size=2.5, parse=T)+ 
##        #xlab(bquote(beta~"of gene expression"))+ylab(bquote(beta~"of gene mean"))+
##        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
##        theme_bw()+
##        theme(#strip.background=element_blank(),
##              #strip.text.x=element_text(size=12),
##              legend.title=element_blank(),
##              legend.text=element_text(size=9),
##              #legend.key.size=unit(0.4,units="cm"),
##              axis.text=element_text(size=9))
##
#####              
##figfn <- "./10_RNA.Variance_output/tmp7/Figure9.1_DMGvsDEG.png"
##png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
##print(fig4)
##dev.off()
##
##x <- dfxy4%>%group_by(MCls, contrast)%>%nest()%>%mutate(ngene=map_dbl(data,nrow))
##
##
#
####
#dfxy5 <- myDFxy(df1, df6)
#
#mycol <- c("1"="grey50", "2"="green", "3"="blue", "4"="red")
#mylabel <- c("1"="NS", "2"="DEG(only)", "3"="DMG2(only)", "4"="Both")
#
#anno_df5 <- dfxy5%>%
#           group_by(contrast, MCls)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#        
#fig5 <- ggplot(dfxy5, aes(x=zscore.x, y=zscore.y))+
#        geom_point(aes(colour=factor(gr)), size=0.3)+
#        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
#        geom_text(data=anno_df5, x=0, y=22, aes(label=eq), size=2.5, parse=T)+
#        xlab("zscore of gene expression(DESeq)")+ylab("zscore of gene mean(Method 2)")+
#        #geom_text(data=anno_df4, x=-1.3, y=3, aes(label=eq), size=2.5, parse=T)+ 
#        #xlab(bquote(beta~"of gene expression"))+ylab(bquote(beta~"of gene mean"))+
#        facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
#        theme_bw()+
#        theme(#strip.background=element_blank(),
#              #strip.text.x=element_text(size=12),
#              legend.title=element_blank(),
#              legend.text=element_text(size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text=element_text(size=9))
#
####              
#figfn <- "./10_RNA.Variance_output/tmp7/Figure9.3_DMG2vsDEG.png"
#png(filename=figfn, width=900, height=800, pointsize=12, res=130)  
#print(fig5)
#dev.off() 
#
#x <- dfxy5%>%group_by(MCls, contrast)%>%nest()%>%mutate(ngene=map_dbl(data,nrow))
#
#
####
#### (4). residual dispersion
#fn <- "./10_RNA.Variance_output/tmp7/tmp1_Monocyte/3_phiNew.meta"
#df4 <- read.table(file=fn, header=T)%>%
#       mutate(zscore=beta/stderr, rn=paste(contrast, gene, sep="_"))%>%drop_na(beta,qval) 
#  
#fn <- "./10_RNA.Variance_output/tmp7/3_phiNew.meta"
#df5 <- read.table(file=fn, header=T)%>%filter(MCls=="Monocyte")%>%
#       mutate(zscore=beta/stderr,rn=paste(contrast, gene, sep="_"))%>%drop_na(beta,qval)
#
#dfxy3 <- myDFxy(df5, df4)
#  
#mycol <- c("1"="grey50", "2"="green", "3"="blue", "4"="red")
#mylabel <- c("1"="NS", "2"="DVG(monocyte)", "3"="DVG(cluster3)", "4"="Both")
#
#anno_df3 <- dfxy3%>%
#           group_by(contrast)%>%
#           nest()%>%
#           mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")),
#                  eq=map(corr,feq))%>%
#           dplyr::select(-data,-corr)#%>%
#        
#fig3 <- ggplot(dfxy3, aes(x=zscore.x, y=zscore.y))+
#        geom_point(aes(colour=factor(gr)), size=0.3)+
#        scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
#        geom_text(data=anno_df3, x=-0.5, y=-10, aes(label=eq), size=2.5, parse=T)+  #-0.5, 28
#        facet_grid(contrast~., labeller=labeller(contrast=Newcon2))+
#        xlab("zscore of dispersion from Monocyte")+ylab("zscore of dispersion from cluster 3")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.text=element_text(size=9),
#              #legend.key.size=unit(0.4,units="cm"),
#              axis.text=element_text(size=9))
#
####              
#figfn <- "./10_RNA.Variance_output/tmp7/tmp1_Monocyte/Figure8.2_DGPmvsDGP3.png"
##png(filename=figfn, width=900, height=800, pointsize=12, res=130) 
#png(filename=figfn, width=600, height=800, pointsize=12, res=130) 
#print(fig3)
#dev.off()
#
#} ### 6.2, End


#####################################
### Part 3, enrichment analysis   ###
#####################################