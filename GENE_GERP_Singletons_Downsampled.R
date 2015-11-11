############################################
### Analyze genes and GERP info together ###
### FOCUS ON SINGLETON SITES             ###
############################################

### Timothy M. Beissinger
### 12-11-2014

### Load singleton info
load("../OUTS/tF_maize_downsampled.Robj")
tF.maize <- allChromosomes
load("../OUTS/tF_teo.Robj")
tF.teo <- allChromosomes

### Compute tF per site
names(tF.maize)[5] <- "tF.win"
names(tF.teo)[5] <- "tF.win"
tF.maize$tF <- tF.maize$tF/tF.maize$nSites
tF.teo$tF <- tF.teo$tF/tF.teo$nSites

### Trim tF matrices
tF.maize.backup <- tF.maize
tF.teo.backup <- tF.teo
tF.maize <- tF.maize[which(tF.maize$nSites>=100),]
tF.teo <- tF.teo[which(tF.teo$nSites>=100),]

#########################################################################
### Determine cM position for each gene by interpolating from panzea map
map <- read.table("../SNPs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs.hmp.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="")
map <- map[,1:5]
ensemblUp <- map[,c(3,4,4)]
#write.table(file="../SNPs/ensemblUp.txt",ensemblUp,quote=F,col.names=F,row.names=F) # upload this file to ensembl to convert to maize v3
ensemblDown <- read.table("../SNPs/ensemblDown.gff",header=F,stringsAsFactors=F,sep="\t")
rem <- which(abs(as.numeric(ensemblUp[,2])-as.numeric(ensemblDown[,4])) > 2000000) #identify positions with massive shifts
map$posV3 <- ensemblDown[,4]
map <- map[-rem,] #remove positions with massive shifts



### Interpolate genetic position for every tF.maize position with info
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
maizeCm <- rep(NA,nrow(tF.maize))
rows <- nrow(tF.maize)

for(i in 1:nrow(tF.maize)){
cat( 100*i/rows, "% done", "\r")
lowerIndex <- which(map$chrom == tF.maize$chr[i] & map$posV3 <= tF.maize$center[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==tF.maize$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == tF.maize$chr[i] & map$posV3 > tF.maize$center[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[tF.maize$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==tF.maize$chr[i])][length(which(map$chrom==tF.maize$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {tF.maize$center[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

maizeCm[i] <- newGen
}

tF.maize$cm <- maizeCm


save.image("GENE_GERP_Singletons_Downsampled.RData")

### Interpolate genetic position for every tF.teo position with info
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
teoCm <- rep(NA,nrow(tF.teo))
rows <- nrow(tF.teo)

for(i in 1:nrow(tF.teo)){
cat( 100*i/rows, "% done", "\r")
lowerIndex <- which(map$chrom == tF.teo$chr[i] & map$posV3 <= tF.teo$center[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==tF.teo$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == tF.teo$chr[i] & map$posV3 > tF.teo$center[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[tF.teo$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==tF.teo$chr[i])][length(which(map$chrom==tF.teo$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {tF.teo$center[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

teoCm[i] <- newGen
}

tF.teo$cm <- teoCm

save.image("GENE_GERP_Singletons_Downsampled.RData")


### Load gene data, previously calculated in Diversity_vs_Gene_Density.R script
load("genesMod.Robj")

### Determine gene center
genesMod$center.cM <- {genesMod$end.cM - genesMod$start.cM}/2 + genesMod$start.cM
names(genesMod)[7]

geneList <- list()
for(i in 1:10){
  geneList[[i]] <- genesMod$center.cM[which(genesMod$chromosome_name == i)]
}

### Determine distance to gene center maize
geneDist<-function(tF,geneList){ min(abs(tF[8] - geneList[[tF[1]]]),na.rm=T)}

tF.maize$distanceToGene <-  apply(tF.maize,1,geneDist,geneList=geneList)

save.image("GENE_GERP_Singletons_Downsampled.RData")

### Determine distance to gene center teosinte
tF.teo$distanceToGene <-  apply(tF.teo,1,geneDist,geneList=geneList)

save.image("GENE_GERP_Singletons_Downsampled.RData")

### Load gerp data, previously calculated in Diversity_vs_GERP.R script
load("gerpElements.Robj")


### Determine gerp element center
gerpElements$center.cM <- {gerpElements$end.cM - gerpElements$start.cM}/2 + gerpElements$start.cM

gerpList <- list()
for(i in 1:10){
  gerpList[[i]] <- gerpElements$center.cM[which(gerpElements$chr == i)]
}

### Determine distance to gerp center maize
gerpDist<-function(tF,gerpList){ min(abs(tF[8] - gerpList[[tF[1]]]),na.rm=T)}

tF.maize$distanceToGerp <-  apply(tF.maize,1,gerpDist,gerpList=gerpList)
save.image("GENE_GERP_Singletons_Downsampled.RData")

### Determine distance to gerp center teosinte
tF.teo$distanceToGerp <-  apply(tF.teo,1,gerpDist,gerpList=gerpList)

save.image("GENE_GERP_Singletons_Downsampled.RData")


### Make new data frame with distance to putative functional element
tF.maize$distanceToElement <- apply(cbind(tF.maize$distanceToGene,tF.maize$distanceToGerp),1,min)
tF.teo$distanceToElement <- apply(cbind(tF.teo$distanceToGene,tF.teo$distanceToGerp),1,min)

### Standardize by neutral diversity
tF_maizeN <- mean(tF.maize$tF[which(tF.maize$distanceToElement>=0.01)],na.rm=T)
tF_teoN <- mean(tF.teo$tF[which(tF.teo$distanceToElement>=0.01)],na.rm=T)

tF.maize$tF.standardized.by.N <- tF.maize$tF/tF_maizeN
tF.teo$tF.standardized.by.N <- tF.teo$tF/tF_teoN

### Neutral-Standardized Spline
TeoSpline.neutral.standardized <- smooth.spline(tF.teo$distanceToElement,tF.teo$tF.standardized.by.N)
MaizeSpline.neutral.standardized <- smooth.spline(tF.maize$distanceToElement,tF.maize$tF.standardized.by.N)

### Plot
plot(NULL,xlim=c(0,0.1),ylim=c(0.4,1.2),xlab="Distance to Gene or Element (cM)",ylab="Singleton Diversity / Neutral Singleton Diversity (1kb windows)")
lines(TeoSpline.neutral.standardized,col="blue",lwd=3)
lines(MaizeSpline.neutral.standardized,col="red",lwd=3)
legend("topright","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(NULL,xlim=c(0,0.002),ylim=c(0.4,1),xlab="",ylab="",main="",cex.axis=0.75)
lines(TeoSpline.neutral.standardized,col="blue",lwd=2)
lines(MaizeSpline.neutral.standardized,col="red",lwd=2)
#text(0.004,-0.55,"Zoom",cex=1.5)



#############################################################################
### Assess significance!           ##########################################
#############################################################################
TeoMat <- matrix(NA,nrow=1000000,ncol=100)
MaizeMat <- matrix(NA,nrow=1000000,ncol=100)
for(i in 1:100){
  print(i)
  maizeBoot <- sample(nrow(tF.maize),replace=T)
  teoBoot <- sample(nrow(tF.teo),replace=T)

  tF.maize.boot <- tF.maize[maizeBoot,]
  tF.teo.boot <- tF.teo[teoBoot,]

  ### Standardize by neutral diversity
  tF_maizeN.boot <- mean(tF.maize.boot$tF[which(tF.maize.boot$distanceToElement>=0.01)])
  tF_teoN.boot <- mean(tF.teo.boot$tF[which(tF.teo.boot$distanceToElement>=0.01)])

  tF.maize.boot$tF.standardized.by.N <- tF.maize.boot$tF/tF_maizeN.boot
  tF.teo.boot$tF.standardized.by.N <- tF.teo.boot$tF/tF_teoN.boot

  ### Neutral-Standardized Spline
  TeoSpline.neutral.standardized.boot <- smooth.spline(tF.teo.boot$distanceToElement,tF.teo.boot$tF.standardized.by.N)
  MaizeSpline.neutral.standardized.boot <- smooth.spline(tF.maize.boot$distanceToElement,tF.maize.boot$tF.standardized.by.N)

  TeoMat[,i] <- predict(TeoSpline.neutral.standardized.boot,seq(0,0.1,length.out=1000000))$y
  MaizeMat[,i] <- predict(MaizeSpline.neutral.standardized.boot,seq(0,0.1,length.out=1000000))$y
}

### Determine 95% intervals
TeoLow <- apply(TeoMat,1,quantile,0.025)
TeoHigh <- apply(TeoMat,1,quantile,0.975)
MaizeLow <- apply(MaizeMat,1,quantile,0.025)
MaizeHigh <- apply(MaizeMat,1,quantile,0.975)

options(scipen=100)




pdf("distanceToElement_WithSignificance_Singletons.pdf",height=10,width=10)
#png("distanceToElement_WithSignificance_Singletons.png")
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.05),ylim=c(0.4,1.1),ylab="Singleton Diversity / Neutral Singleton Diversity",xlab="Distance to Gene or Element (cM)")
polygon(c(x,rev(x)), c(TeoLow,rev(TeoHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeLow,rev(MaizeHigh)),col=rgb(1, 0, 0,0.5),border=NA)
observedTeo <- predict(TeoSpline.neutral.standardized,seq(0,0.1,length.out=1000000))$y
observedMaize <- predict(MaizeSpline.neutral.standardized,seq(0,0.1,length.out=1000000))$y
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
legend("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.4,0.9),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
axis(2,c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1),at=c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1))
polygon(c(x,rev(x)), c(TeoLow,rev(TeoHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeLow,rev(MaizeHigh)),col=rgb(1, 0, 0,0.5),border=NA)
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
dev.off()

save.image("GENE_GERP_Singletons_Downsampled.RData")





####################################################################################
####################################################################################
### Repeat for distanct to gene -- more interesting                              ###
####################################################################################
####################################################################################

### Standardize by neutral diversity
tF_maize_gene_N <- mean(tF.maize$tF[which(tF.maize$distanceToGene>=0.02)],na.rm=T)
tF_teo_gene_N <- mean(tF.teo$tF[which(tF.teo$distanceToGene>=0.02)],na.rm=T)

tF.maize$tF.standardized.by.gene.N <- tF.maize$tF/tF_maize_gene_N
tF.teo$tF.standardized.by.gene.N <- tF.teo$tF/tF_teo_gene_N

### Neutral-Standardized Spline
TeoSpline.neutral.gene.standardized <- smooth.spline(tF.teo$distanceToGene,tF.teo$tF.standardized.by.gene.N)
MaizeSpline.neutral.gene.standardized <- smooth.spline(tF.maize$distanceToGene,tF.maize$tF.standardized.by.gene.N)

### SAVE FOR B-Calcs
#save(TeoSpline.neutral.gene.standardized,MaizeSpline.neutral.gene.standardized,file="../../Bcalcs/splines_Downsampled.Robj")

### Plot
plot(NULL,xlim=c(0,0.1),ylim=c(0.4,1.2),xlab="Distance to Gene or Element (cM)",ylab="Singleton Diversity / Neutral Singleton Diversity (1kb windows)")
lines(TeoSpline.neutral.gene.standardized,col="blue",lwd=3)
lines(MaizeSpline.neutral.gene.standardized,col="red",lwd=3)
legend("topright","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(NULL,xlim=c(0,0.0005),ylim=c(0.4,1),xlab="",ylab="",main="",cex.axis=0.75)
lines(TeoSpline.neutral.gene.standardized,col="blue",lwd=2)
lines(MaizeSpline.neutral.gene.standardized,col="red",lwd=2)
#text(0.004,-0.55,"Zoom",cex=1.5)



#############################################################################
### Assess significance!           ##########################################
#############################################################################
TeoGeneMat <- matrix(NA,nrow=1000000,ncol=100)
MaizeGeneMat <- matrix(NA,nrow=1000000,ncol=100)
for(i in 1:100){
  print(i)
  maizeBoot <- sample(nrow(tF.maize),replace=T)
  teoBoot <- sample(nrow(tF.teo),replace=T)

  tF.maize.boot <- tF.maize[maizeBoot,]
  tF.teo.boot <- tF.teo[teoBoot,]

  ### Standardize by neutral diversity
  tF_maize_gene_N.boot <- mean(tF.maize.boot$tF[which(tF.maize.boot$distanceToGene>=0.02)])
  tF_teo_gene_N.boot <- mean(tF.teo.boot$tF[which(tF.teo.boot$distanceToGene>=0.02)])

  tF.maize.boot$tF.standardized.by.gene.N <- tF.maize.boot$tF/tF_maize_gene_N.boot
  tF.teo.boot$tF.standardized.by.gene.N <- tF.teo.boot$tF/tF_teo_gene_N.boot

  ### Neutral-Standardized Spline
  TeoSpline.neutral.gene.standardized.boot <- smooth.spline(tF.teo.boot$distanceToGene,tF.teo.boot$tF.standardized.by.gene.N)
  MaizeSpline.neutral.gene.standardized.boot <- smooth.spline(tF.maize.boot$distanceToGene,tF.maize.boot$tF.standardized.by.gene.N)

  TeoGeneMat[,i] <- predict(TeoSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=1000000))$y
  MaizeGeneMat[,i] <- predict(MaizeSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=1000000))$y
}

### Determine 95% intervals
TeoGeneLow <- apply(TeoGeneMat,1,quantile,0.025)
TeoGeneHigh <- apply(TeoGeneMat,1,quantile,0.975)
MaizeGeneLow <- apply(MaizeGeneMat,1,quantile,0.025)
MaizeGeneHigh <- apply(MaizeGeneMat,1,quantile,0.975)

options(scipen=100)




#pdf("distanceToGene_WithSignificance_Singletons.pdf",height=10,width=10)
png("distanceToGene_WithSignificance_Singletons_Downsampled.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.05),ylim=c(0.4,1.1),ylab="Singleton Diversity / Neutral Singleton Diversity",xlab="Distance to Gene (cM)",main="Singleton diversity surrounding genes")
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
observedTeo <- predict(TeoSpline.neutral.gene.standardized,seq(0,0.1,length.out=1000000))$y
observedMaize <- predict(MaizeSpline.neutral.gene.standardized,seq(0,0.1,length.out=1000000))$y
lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
legend("bottomleft","(x,y)",col=c("red","blue"),c("Maize","Teosinte"),lwd=3)
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.0005),ylim=c(0.4,0.9),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0001,0.0002,0.0003,0.0004,0.0005),at=c(0,0.0001,0.0002,0.0003,0.0004,0.0005))
axis(2,c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1),at=c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1))
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.5),border=NA)
#lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
dev.off()

save.image("GENE_GERP_Singletons_Downsampled.RData")


################ Manuscript figs ###################
load("SingletonSplines.Robj")

options(scipen=100)
#pdf("distanceToGene_WithSignificance_Singletons_manuscript.pdf",height=10,width=10)
png("distanceToGene_WithSignificance_Singletons_Downsampled_threeLines_manuscript.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=1000000)
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.03),ylim=c(0.4,1.05),ylab="Singleton Diversity / Neutral Singleton Diversity",xlab="Distance to Gene (cM)")#,main="Singleton diversity surrounding genes")
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.25),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow_Single,rev(MaizeGeneHigh_Single)),col=rgb(1, 0, 0,0.5),border=NA) #not downsampled
observedTeo <- predict(TeoSpline.neutral.gene.standardized,seq(0,0.1,length.out=1000000))$y
observedMaize <- predict(MaizeSpline.neutral.gene.standardized,seq(0,0.1,length.out=1000000))$y
lines(x,observedTeo,col=rgb(0, 0, 1,1),lwd=2)
lines(x,observedMaize,col=rgb(1, 0, 0,0.25),lwd=2)
lines(x,observedMaize_Single,col=rgb(1,0,0,1),lwd=2)
text(x=0.002,y=1,labels="B",cex=4)
source("legend.v2.R")
legend.v2("bottomleft","(x,y)",col=c("red",rgb(1, 0, 0,0.25),"blue"),c("Maize","Downsampled Maize","Teosinte"),lwd=2,fill=c(rgb(1,0,0,0.5), rgb(1, 0, 0,0.25),rgb(0, 0, 1,0.5)),bty="n",border=c(NA,NA,NA))
par(new = TRUE)
par(fig = c(0.47, 0.97, 0.1, 0.6))
plot(x,TeoGeneLow,col="blue",type="n",xlim=c(0,0.002),ylim=c(0.4,0.9),ylab="",xlab="",xaxt="n",yaxt="n")
axis(1,c(0,0.0005,0.001,0.0015,0.002),at=c(0,0.0005,0.001,0.0015,0.002))
axis(2,c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1),at=c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1))
polygon(c(x,rev(x)), c(TeoGeneLow,rev(TeoGeneHigh)),col=rgb(0, 0, 1,0.5),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow,rev(MaizeGeneHigh)),col=rgb(1, 0, 0,0.25),border=NA)
polygon(c(x,rev(x)), c(MaizeGeneLow_Single,rev(MaizeGeneHigh_Single)),col=rgb(1, 0, 0,0.5),border=NA) #not downsampled
#lines(x,observedTeo,col=rgb(0, 0, 1,0.75),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,0.75),lwd=2)
#lines(x,observedTeo,col=rgb(0, 0, 1,1),lwd=2)
#lines(x,observedMaize,col=rgb(1, 0, 0,0.25),lwd=2)
#lines(x,observedMaize_Single,col=rgb(1,0,0,1),lwd=2)
dev.off()





save.image("GENE_GERP_Singletons_Downsampled.RData")
