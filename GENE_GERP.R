############################################
### Analyze genes and GERP info together ###
############################################

### Timothy M. Beissinger
### 11-17-2014

### Load gene data
load("Diversity_vs_Gene_Density.RData")
div0MaizeGenes <- div0Maize
div0TeoGenes <- div0Teo
rm(list=ls()[which(ls() != "div0MaizeGenes" & ls() != "div0TeoGenes")])

### Load GERP data
load("Diversity_vs_GERP.RData")
div0MaizeGerp <- div0Maize
div0TeoGerp <- div0Teo
rm(list=ls()[which(ls() != "div0MaizeGenes" & ls() != "div0TeoGenes" & ls() != "div0MaizeGerp" & ls() != "div0TeoGerp")])

### Make new data frame with distance to putative functional element
#div0Maize <- div0MaizeGenes[,c(1,2,3,6,4,5,13)]
div0Maize <- div0MaizeGenes[1:9]
#div0Teo <- div0TeoGenes[,c(1,2,3,6,4,5,13)]
div0Teo <- div0TeoGenes[1:9]
div0Maize$distanceToGerp <- div0MaizeGerp$distanceToGerp
div0Teo$distanceToGerp <- div0TeoGerp$distanceToGerp

div0Maize$distanceToElement <- apply(cbind(div0Maize$distanceToGene,div0Maize$distanceToGerp),1,min)
div0Teo$distanceToElement <- apply(cbind(div0Teo$distanceToGene,div0Teo$distanceToGerp),1,min)


### Standardize by neutral diversity
pi_maizeN <- mean(div0Maize$tP[which(div0Maize$distanceToElement>=0.01)])
pi_teoN <- mean(div0Teo$tP[which(div0Teo$distanceToElement>=0.01)])

div0Maize$tP.standardized.by.N <- div0Maize$tP/pi_maizeN
div0Teo$tP.standardized.by.N <- div0Teo$tP/pi_teoN

### Neutral-Standardized Spline
TeoSpline.neutral.standardized <- smooth.spline(div0Teo$distanceToElement,div0Teo$tP.standardized.by.N)
MaizeSpline.neutral.standardized <- smooth.spline(div0Maize$distanceToElement,div0Maize$tP.standardized.by.N)


### Plot with inset
#plot(div0Teo$distanceToGene,div0Teo$tP.standardized.by.N,xlim=c(0,0.2),col="blue",cex=0.25)
#points(div0Maize$distanceToGene,div0Maize$tP.standardized.by.N,xlim=c(0,0.2),col="red",cex=0.25)
#png("distanceToElement_neutral-standardized.png",height=8,width=8,units="in",res=300)
plot(NULL,xlim=c(0,0.01),ylim=c(0.3,1.2),xlab="Distance to Gene or Element (cM)",ylab="Diversity / Neutral Diversity (1kb windows)")
lines(TeoSpline.neutral.standardized,col="blue",lwd=3)
lines(MaizeSpline.neutral.standardized,col="red",lwd=3)
legend("topright","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(NULL,xlim=c(0,0.002),ylim=c(0.5,1.3),xlab="",ylab="",main="",cex.axis=0.75)
lines(TeoSpline.neutral.standardized,col="blue",lwd=2)
lines(MaizeSpline.neutral.standardized,col="red",lwd=2)
#text(0.004,-0.55,"Zoom",cex=1.5)

#dev.off()




### Look at ratio of the two
ratio <- predict(TeoSpline.neutral.standardized,seq(0,0.1,length.out=100000))$y / predict(MaizeSpline.neutral.standardized,seq(0,0.1,length.out=100000))$y
plot(seq(0,0.1,length.out=100000),ratio,type="l",xlim=c(0,0.005))


#############################################################################
### Assess significance!           ##########################################
#############################################################################
TeoMat <- matrix(NA,nrow=1000000,ncol=100)
MaizeMat <- matrix(NA,nrow=1000000,ncol=100)
for(i in 1:100){
    print(i)
    maizeBoot <- sample(nrow(div0Maize),replace=T)
    teoBoot <- sample(nrow(div0Teo),replace=T)

    div0Maize.boot <- div0Maize[maizeBoot,]
    div0Teo.boot <- div0Teo[teoBoot,]

    ### Standardize by neutral diversity
    pi_maizeN.boot <- mean(div0Maize.boot$tP[which(div0Maize.boot$distanceToElement>=0.01)])
    pi_teoN.boot <- mean(div0Teo.boot$tP[which(div0Teo.boot$distanceToElement>=0.01)])

    div0Maize.boot$tP.standardized.by.N <- div0Maize.boot$tP/pi_maizeN.boot
    div0Teo.boot$tP.standardized.by.N <- div0Teo.boot$tP/pi_teoN.boot

    ### Neutral-Standardized Spline
    TeoSpline.neutral.standardized.boot <- smooth.spline(div0Teo.boot$distanceToElement,div0Teo.boot$tP.standardized.by.N)
    MaizeSpline.neutral.standardized.boot <- smooth.spline(div0Maize.boot$distanceToElement,div0Maize.boot$tP.standardized.by.N)

    TeoMat[,i] <- predict(TeoSpline.neutral.standardized.boot,seq(0,0.1,length.out=1000000))$y
    MaizeMat[,i] <- predict(MaizeSpline.neutral.standardized.boot,seq(0,0.1,length.out=1000000))$y
}

### Determine 95% intervals
TeoLow <- apply(TeoMat,1,quantile,0.025)
TeoHigh <- apply(TeoMat,1,quantile,0.975)
MaizeLow <- apply(MaizeMat,1,quantile,0.025)
MaizeHigh <- apply(MaizeMat,1,quantile,0.975)

options(scipen=100)

pdf("distanceToElement_WithSignificance.pdf",height=10,width=10)
#png("distanceToElement_WithSignificance_Folded2.png")
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.003),ylim=c(0.5,1.1),ylab="Diversity / Neutral Diversity (1kb Windows)",xlab="Distance to Gene or Element (cM)")
#lines(x,TeoHigh,col="blue")
#lines(x,MaizeLow,col="red")
#lines(x,MaizeHigh,col="red")
polygon(c(x,rev(x)), c(TeoLow,rev(TeoHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeLow,rev(MaizeHigh)),col=rgb(1, 0, 0,0.5),border=NA)
observedTeo <- predict(TeoSpline.neutral.standardized,seq(0,0.1,length.out=1000000))$y
observedMaize <- predict(MaizeSpline.neutral.standardized,seq(0,0.1,length.out=1000000))$y
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
legend("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.6,1.1),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
axis(2,c(0.6,0.7,0.8,0.9,1,1.1),at=c(0.6,0.7,0.8,0.9,1,1.1))
polygon(c(x,rev(x)), c(TeoLow,rev(TeoHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeLow,rev(MaizeHigh)),col=rgb(1, 0, 0,0.5),border=NA)
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
dev.off()

save.image("GENE_GERP.RData")







############################################################################################
### Analyze Genes only                                                                   ###
############################################################################################


### Standardize by neutral diversity
tP_maize_gene_N <- mean(div0Maize$tP[which(div0Maize$distanceToGene>=0.01)],na.rm=T)
tP_teo_gene_N <- mean(div0Teo$tP[which(div0Teo$distanceToGene>=0.01)],na.rm=T)

div0Maize$tP.standardized.by.gene.N <- div0Maize$tP/tP_maize_gene_N
div0Teo$tP.standardized.by.gene.N <- div0Teo$tP/tP_teo_gene_N

### Neutral-Standardized Spline
TeoSpline.neutral.gene.standardized <- smooth.spline(div0Teo$distanceToGene,div0Teo$tP.standardized.by.gene.N)
MaizeSpline.neutral.gene.standardized <- smooth.spline(div0Maize$distanceToGene,div0Maize$tP.standardized.by.gene.N)

### Plot
plot(NULL,xlim=c(0,0.1),ylim=c(0.4,1.2),xlab="Distance to Gene or Element (cM)",ylab="Singleton Diversity / Neutral Singleton Diversity (1kb windows)")
lines(TeoSpline.neutral.gene.standardized,col="blue",lwd=3)
lines(MaizeSpline.neutral.gene.standardized,col="red",lwd=3)
legend("topright","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(NULL,xlim=c(0,0.002),ylim=c(0.4,1),xlab="",ylab="",main="",cex.axis=0.75)
lines(TeoSpline.neutral.gene.standardized,col="blue",lwd=2)
lines(MaizeSpline.neutral.gene.standardized,col="red",lwd=2)
#text(0.004,-0.55,"Zoom",cex=1.5)

### Play with rollapply instead of cubic smoothing spline.
maizeData <- cbind(div0Maize$distanceToGene,div0Maize$tP.standardized.by.gene.N)
maizeData <- maizeData[order(maizeData[,1]),]
library(zoo)
rolled<-rollapply(maizeData[,2],width=1000,mean,na.pad=T)
rolledDis<-rollapply(maizeData[,1],width=1000,mean,na.pad=T)

#############################################################################
### Assess significance!           ##########################################
#############################################################################
TeoGeneMat <- matrix(NA,nrow=1000000,ncol=100)
MaizeGeneMat <- matrix(NA,nrow=1000000,ncol=100)
for(i in 1:100){
  print(i)
  maizeBoot <- sample(nrow(div0Maize),replace=T)
  teoBoot <- sample(nrow(div0Teo),replace=T)

  div0Maize.boot <- div0Maize[maizeBoot,]
  div0Teo.boot <- div0Teo[teoBoot,]

  ### Standardize by neutral diversity
  tP_maize_gene_N.boot <- mean(div0Maize.boot$tP[which(div0Maize.boot$distanceToGene>=0.01)])
  tP_teo_gene_N.boot <- mean(div0Teo.boot$tP[which(div0Teo.boot$distanceToGene>=0.01)])

  div0Maize.boot$tP.standardized.by.gene.N <- div0Maize.boot$tP/tP_maize_gene_N.boot
  div0Teo.boot$tP.standardized.by.gene.N <- div0Teo.boot$tP/tP_teo_gene_N.boot

  ### Neutral-Standardized Spline
  TeoSpline.neutral.gene.standardized.boot <- smooth.spline(div0Teo.boot$distanceToGene,div0Teo.boot$tP.standardized.by.gene.N)
  MaizeSpline.neutral.gene.standardized.boot <- smooth.spline(div0Maize.boot$distanceToGene,div0Maize.boot$tP.standardized.by.gene.N)

  TeoGeneMat[,i] <- predict(TeoSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=1000000))$y
  MaizeGeneMat[,i] <- predict(MaizeSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=1000000))$y
}

### Determine 95% intervals
TeoGeneLow <- apply(TeoGeneMat,1,quantile,0.025)
TeoGeneHigh <- apply(TeoGeneMat,1,quantile,0.975)
MaizeGeneLow <- apply(MaizeGeneMat,1,quantile,0.025)
MaizeGeneHigh <- apply(MaizeGeneMat,1,quantile,0.975)

options(scipen=100)




pdf("distanceToGene_WithSignificance_Folded2.pdf",height=10,width=10,pointsize=18)
#png("distanceToGene_WithSignificance_Folded2.png",height=10,width=10,units="in",res=300)
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.01),ylim=c(0.6,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="Pairwise diversity surrounding genes")
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
observedTeo <- predict(TeoSpline.neutral.gene.standardized,x)$y
observedMaize <- predict(MaizeSpline.neutral.gene.standardized,x)$y
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
legend("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.6,1),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
axis(2,c(0.6,0.7,0.8,0.9,1),at=c(0.6,0.7,0.8,0.9,1))
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
#lines(observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
dev.off()



png("distanceToGene_WithSignificance_MaizeOnly_Folded2.png",height=10,width=10,pointsize=18,units="in",res=300)
#png("distanceToGene_WithSignificance_Folded2.png",height=10,width=10,units="in",res=300)
x <- seq(0,0.1,length.out=1000000)
plot(x,MaizeGeneLow,col="blue",type="n",xlim=c(0,0.01),ylim=c(0.6,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="Pairwise maize diversity surrounding genes")
#polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
#observedTeo <- predict(TeoSpline.neutral.gene.standardized,x)$y
observedMaize <- predict(MaizeSpline.neutral.gene.standardized,x)$y
#lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
#legend("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
#par(new = TRUE)
#par(fig = c(0.47, 0.97, 0.1, 0.6))
#plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.6,1),ylab="",xlab="",xaxt="n",yaxt="n")
#axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
#axis(2,c(0.6,0.7,0.8,0.9,1),at=c(0.6,0.7,0.8,0.9,1))
#polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
#polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
#lines(observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
dev.off()


save(TeoGeneLow,TeoGeneHigh,MaizeGeneLow,MaizeGeneHigh,TeoSpline.neutral.gene.standardized,MaizeSpline.neutral.gene.standardized,file="pi_gene_objects.Robj")

save.image("GENE_GERP.RData")

####################### Manuscript Figures #####################


options(scipen=100)
#pdf("distanceToGene_WithSignificance_Folded2_manuscript.pdf",height=10,width=10,pointsize=18)
png("distanceToGene_WithSignificance_Folded2_manuscript.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.01),ylim=c(0.4,1.05),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)")#,main="Pairwise diversity surrounding genes")
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
observedTeo <- predict(TeoSpline.neutral.gene.standardized,x)$y
observedMaize <- predict(MaizeSpline.neutral.gene.standardized,x)$y
lines(x,observedTeo,col=rgb(0, 0, 1,1),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,1),lwd=2)
text(x=0.0005,y=1.02,labels="A",cex=4)
source("legend.v2.R")
legend.v2("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=2,fill=c(rgb(1,0,0,0.5),rgb(0, 0, 1,0.5)), bty="n",border=c(NA,NA,NA))
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.6,1),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
axis(2,c(0.6,0.7,0.8,0.9,1),at=c(0.6,0.7,0.8,0.9,1))
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
#lines(x,observedTeo,col=rgb(0, 0, 1,1),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,1),lwd=2)
dev.off()
