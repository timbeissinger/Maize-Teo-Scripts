#####################################################################
### Use this script to plot Singletons vs. Pi in Maize and Teo    ###
#####################################################################


### Timothy M. Beissinger
### 5-15-2015
### updated 8/17/2015

####################################################################
####################################################################
####################################################################
### Load Pi data and re-fit to neutral >= 0.02 cM
####################################################################
####################################################################
####################################################################

load("GENE_GERP.RData")



### Standardize by neutral diversity
tP_maize_gene_N <- mean(div0Maize$tP[which(div0Maize$distanceToGene>=0.02)],na.rm=T) # NOTICE THE 0.02!!! Different from regular script!!!
tP_teo_gene_N <- mean(div0Teo$tP[which(div0Teo$distanceToGene>=0.02)],na.rm=T)# NOTICE THE 0.02!!! Different from regular script!!!

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




#############################################################################
### Assess significance!           ##########################################
#############################################################################
TeoGeneMat <- matrix(NA,nrow=100000,ncol=100)
MaizeGeneMat <- matrix(NA,nrow=100000,ncol=100)
for(i in 1:100){
  print(i)
  maizeBoot <- sample(nrow(div0Maize),replace=T)
  teoBoot <- sample(nrow(div0Teo),replace=T)

  div0Maize.boot <- div0Maize[maizeBoot,]
  div0Teo.boot <- div0Teo[teoBoot,]

  ### Standardize by neutral diversity
  tP_maize_gene_N.boot <- mean(div0Maize.boot$tP[which(div0Maize.boot$distanceToGene>=0.02)])
  tP_teo_gene_N.boot <- mean(div0Teo.boot$tP[which(div0Teo.boot$distanceToGene>=0.02)])

  div0Maize.boot$tP.standardized.by.gene.N <- div0Maize.boot$tP/tP_maize_gene_N.boot
  div0Teo.boot$tP.standardized.by.gene.N <- div0Teo.boot$tP/tP_teo_gene_N.boot

  ### Neutral-Standardized Spline
  TeoSpline.neutral.gene.standardized.boot <- smooth.spline(div0Teo.boot$distanceToGene,div0Teo.boot$tP.standardized.by.gene.N)
  MaizeSpline.neutral.gene.standardized.boot <- smooth.spline(div0Maize.boot$distanceToGene,div0Maize.boot$tP.standardized.by.gene.N)

  TeoGeneMat[,i] <- predict(TeoSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=100000))$y
  MaizeGeneMat[,i] <- predict(MaizeSpline.neutral.gene.standardized.boot,seq(0,0.1,length.out=100000))$y
}

###Checkpoint###
save.image("SingletonsVsPi.RData")

### Determine 95% intervals
TeoGeneLow <- apply(TeoGeneMat,1,quantile,0.025)
TeoGeneHigh <- apply(TeoGeneMat,1,quantile,0.975)
MaizeGeneLow <- apply(MaizeGeneMat,1,quantile,0.025)
MaizeGeneHigh <- apply(MaizeGeneMat,1,quantile,0.975)

options(scipen=100)



### Append splines so it is clear that these are pi
TeoGeneLow_Pi <- TeoGeneLow
TeoGeneHigh_Pi <- TeoGeneHigh
MaizeGeneLow_Pi <- MaizeGeneLow
MaizeGeneHigh_Pi <- MaizeGeneHigh
TeoSpline.neutral.gene.standardized_Pi <- TeoSpline.neutral.gene.standardized
MaizeSpline.neutral.gene.standardized_Pi <- MaizeSpline.neutral.gene.standardized

### Load Singleton curves
load("GENE_GERP_Singletons.RData")

### Append splines so it is clear that these are pi
TeoGeneLow_Single <- TeoGeneLow
TeoGeneHigh_Single <- TeoGeneHigh
MaizeGeneLow_Single <- MaizeGeneLow
MaizeGeneHigh_Single <- MaizeGeneHigh
TeoSpline.neutral.gene.standardized_Single <- TeoSpline.neutral.gene.standardized
MaizeSpline.neutral.gene.standardized_Single <- MaizeSpline.neutral.gene.standardized



################ Maize Figure #################
png("distanceToGene_WithSignificance_Folded2_maizeSingleVsPi.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=100000)
plot(x,MaizeGeneLow_Pi,col="blue",type="n",xlim=c(0,0.02),ylim=c(0.4,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="Maize Diversity Away from Genes")
polygon(c(x,rev(x)), c(MaizeGeneLow_Pi,rev(MaizeGeneHigh_Pi)),col=rgb(1, 0, 0,0.8),border=NA)
x <- seq(0,0.1,length.out=1000000)
polygon(c(x,rev(x)), c(MaizeGeneLow_Single,rev(MaizeGeneHigh_Single)),col=rgb(1, 0, 0,0.5),border=NA)
observedMaize_Pi <- predict(MaizeSpline.neutral.gene.standardized_Pi,x)$y
observedMaize_Single <- predict(MaizeSpline.neutral.gene.standardized_Single,x)$y
lines(x,observedMaize_Pi,col=rgb(1, 0, 0,1),lwd=2)
lines(x,observedMaize_Single,col=rgb(1, 0, 0,0.75),lwd=2)
legend("bottomright","(x,y)",col=c(rgb(1, 0, 0,0.9),rgb(1, 0, 0,0.5)),c("Pi","Singletons"),lwd=3)
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

################ Teosinte Figure #################
png("distanceToGene_WithSignificance_Folded2_teoSingleVsPi.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=100000)
plot(x,TeoGeneLow_Pi,col="blue",type="n",xlim=c(0,0.02),ylim=c(0.4,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="Teosinte Diversity Away from Genes")
polygon(c(x,rev(x)), c(TeoGeneLow_Pi,rev(TeoGeneHigh_Pi)),col=rgb(0, 0, 1,0.8),border=NA)
x <- seq(0,0.1,length.out=1000000)
polygon(c(x,rev(x)), c(TeoGeneLow_Single,rev(TeoGeneHigh_Single)),col=rgb(0, 0, 1,0.5),border=NA)
observedTeo_Pi <- predict(TeoSpline.neutral.gene.standardized_Pi,x)$y
observedTeo_Single <- predict(TeoSpline.neutral.gene.standardized_Single,x)$y
lines(x,observedTeo_Pi,col=rgb(0, 0, 1,1),lwd=2)
lines(x,observedTeo_Single,col=rgb(0, 0, 1,0.75),lwd=2)
legend("bottomright","(x,y)",col=c(rgb(0, 0, 1,0.9),rgb(0, 0, 1,0.5)),c("Pi","Singletons"),lwd=3)
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






##########################################
##### Everything Figure            #######
##########################################

png("distanceToGene_WithSignificance_Folded2_maizeAndTeoSingleVsPi.png",height=10,width=10,units="in",res=300,pointsize=18)
x <- seq(0,0.1,length.out=100000)
plot(x,MaizeGeneLow_Pi,col="blue",type="n",xlim=c(0,0.01),ylim=c(0.4,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="")
polygon(c(x,rev(x)), c(MaizeGeneLow_Pi,rev(MaizeGeneHigh_Pi)),col=rgb(1, 0, 0,0.8),border=NA)
x <- seq(0,0.1,length.out=1000000)
polygon(c(x,rev(x)), c(MaizeGeneLow_Single,rev(MaizeGeneHigh_Single)),col=rgb(1, 0, 0,0.4),border=NA)
observedMaize_Pi <- predict(MaizeSpline.neutral.gene.standardized_Pi,x)$y
observedMaize_Single <- predict(MaizeSpline.neutral.gene.standardized_Single,x)$y
lines(x,observedMaize_Pi,col=rgb(1, 0, 0,1),lwd=2)
lines(x,observedMaize_Single,col=rgb(1, 0, 0,0.75),lwd=2)


### Teo
x <- seq(0,0.1,length.out=100000)
#plot(x,TeoGeneLow_Pi,col="blue",type="n",xlim=c(0,0.02),ylim=c(0.4,1.1),ylab="Diversity / Neutral Diversity",xlab="Distance to Gene (cM)",main="Teosinte Diversity Away from Genes")
polygon(c(x,rev(x)), c(TeoGeneLow_Pi,rev(TeoGeneHigh_Pi)),col=rgb(0, 0, 1,0.8),border=NA)
x <- seq(0,0.1,length.out=1000000)
polygon(c(x,rev(x)), c(TeoGeneLow_Single,rev(TeoGeneHigh_Single)),col=rgb(0, 0, 1,0.4),border=NA)
observedTeo_Pi <- predict(TeoSpline.neutral.gene.standardized_Pi,x)$y
observedTeo_Single <- predict(TeoSpline.neutral.gene.standardized_Single,x)$y
lines(x,observedTeo_Pi,col=rgb(0, 0, 1,1),lwd=2)
lines(x,observedTeo_Single,col=rgb(0, 0, 1,0.75),lwd=2)

source("legend.v2.R")
legend.v2("bottomright","(x,y)",col=c("red","blue","red","blue"),fill=c(rgb(1, 0, 0,0.8),rgb(0, 0, 1,0.8),rgb(1, 0, 0,0.4),rgb(0, 0, 1,0.4)),c("Pi maize","Pi teosinte","Singletons maize","Singletons teosinto"),lwd=3,bty="n",border=NA)

dev.off()








save.image("SingletonsVsPi.RData")
