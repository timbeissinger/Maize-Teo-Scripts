################################################################################
### Use this script to make a plot of diversity surrounding Syn, non-syn     ###
### substitutions between tripsicum and maize.  Look specifically at the     ###
### 10% most and 10% least diverse genes. This script was adapted from       ###
### plotDiversity_TvM.R in June, 2015.                                       ###
################################################################################

### 6/22/2015

### Load effects estimates
effects <- read.table("../SNPs/TvMeffects.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="",na.strings="-",skip=8)
levels(as.factor(effects$Consequence))

effects <- effects[,c(1,2,3,4,5,6,7)]

### Make a set of syn, non variants
syn <- effects[which(effects$Consequence == "synonymous_variant"),]
mis <- effects[which(effects$Consequence == "missense_variant"),]
int <- effects[which(effects$Consequence == "intergenic_variant"),]

### Remove ambiguous subs (positions in both syn, mis) #no ambiguous subs in int
amb <- intersect(syn$Location,mis$Location)
syn <- syn[-which(syn$Location %in% amb),]
mis <- mis[-which(mis$Location %in% amb),]


### Remove duplicate positions syn (multiple transcripts)
syn0 <- syn[NULL,]
levels <- levels(as.factor(syn$Location))
nlevels <- length(levels(as.factor(syn$Location)))

for(i in 1:nlevels){
print(i)
uniqueRow <- which(syn$Location==levels[i])[1]
syn0[nrow(syn0)+1,] <- syn[uniqueRow,]
}

### Remove duplicate positions mis (multiple transcripts)
mis0 <- mis[NULL,]
levels <- levels(as.factor(mis$Location))
nlevels <- length(levels(as.factor(mis$Location)))

for(i in 1:nlevels){
print(i)
uniqueRow <- which(mis$Location==levels[i])[1]
mis0[nrow(mis0)+1,] <- mis[uniqueRow,]
}


### There are no duplicate positions in int
int0 <- int

### Put syn0m, mis0, and int0 in order
options(scipen=10)
syn0$chr <- as.numeric(unlist(strsplit(syn0$Location,split=":"))[seq(1,2*nrow(syn0),2)])
syn0$pos <- as.numeric(unlist(strsplit(syn0$Location,split=":"))[seq(2,2*nrow(syn0),2)])

mis0$chr <- as.numeric(unlist(strsplit(mis0$Location,split=":"))[seq(1,2*nrow(mis0),2)])
mis0$pos <- as.numeric(unlist(strsplit(mis0$Location,split=":"))[seq(2,2*nrow(mis0),2)])

int0$chr <- as.numeric(unlist(strsplit(int0$Location,split=":"))[seq(1,2*nrow(int0),2)])
int0$pos <- as.numeric(unlist(strsplit(int0$Location,split=":"))[seq(2,2*nrow(int0),2)])

syn0 <- syn0[order(syn0$chr,syn0$pos),]
mis0 <- mis0[order(mis0$chr,mis0$pos),]
int0 <- int0[order(int0$chr,int0$pos),]


### Load genetic map
map <- read.table("../SNPs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs/NAM_phasedImputed_1cM_AllZeaGBSv2.3_allChrs.hmp.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="")
map <- map[,1:5]
ensemblUp <- map[,c(3,4,4)]
#write.table(file="../SNPs/ensemblUp.txt",ensemblUp,quote=F,col.names=F,row.names=F) # upload this file to ensembl to convert to maize v3
ensemblDown <- read.table("../SNPs/ensemblDown.gff",header=F,stringsAsFactors=F,sep="\t")

rem <- which(abs(as.numeric(ensemblUp[,2])-as.numeric(ensemblDown[,4])) > 2000000) #identify positions with massive shifts

map$posV3 <- ensemblDown[,4]

map <- map[-rem,] #remove positions with massive shifts

### CHECKPOINT ###
save.image("plotDiversity_Weak_Strong.RData")

### Interpolate genetic position for every syn0 SNP
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
syn0$cm <- NA

for(i in 1:nrow(syn0)){
print(i)
lowerIndex <- which(map$chrom == syn0$chr[i] & map$posV3 <= syn0$pos[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==syn0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == syn0$chr[i] & map$posV3 >= syn0$pos[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[syn0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==syn0$chr[i])][length(which(map$chrom==syn0$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {syn0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

syn0$cm[i] <- newGen
}


### Interpolate genetic position for every mis0 SNP
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
mis0$cm <- NA

for(i in 1:nrow(mis0)){
print(i)
lowerIndex <- which(map$chrom == mis0$chr[i] & map$posV3 <= mis0$pos[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==mis0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == mis0$chr[i] & map$posV3 >= mis0$pos[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[mis0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==mis0$chr[i])][length(which(map$chrom==mis0$chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {mis0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

mis0$cm[i] <- newGen
}


# ### Interpolate genetic position for every int0 SNP
# chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
# int0$cm <- NA
# 
# for(i in 1:nrow(int0)){
# print(i)
# lowerIndex <- which(map$chrom == int0$chr[i] & map$posV3 <= int0$pos[i]) #find index of map anchors smaller than observed position
# belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
# belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==int0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position
# 
# higherIndex <- which(map$chrom == int0$chr[i] & map$posV3 >= int0$pos[i]) #find index of map anchors larger than observed position
# abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[int0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
# aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==int0$chr[i])][length(which(map$chrom==int0$chr[i]))]+1,na.rm=T) #take corresponding genetic position
# 
# scale <- {int0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
# newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position
# 
# int0$cm[i] <- newGen
# }


### CHECKPOINT ###
save.image("plotDiversity_Weak_Strong.RData")


### Load diversity data (from angsd) [December: updated to whole genome FOLDED--more data]
diversity <- read.table("../../WholeGenomeFolded2/OUTS/BKN_WholeGenome_windows.thetas.gz.pestPG",comment.char="",skip=1,stringsAsFactors=F,header=T)
div0 <- diversity[,c(2,3,5,14)]

### Interpolate genetic position for every angsd position with diversity info.
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
divCm <- rep(NA,nrow(div0))
rows <- nrow(div0)

for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
lowerIndex <- which(map$chrom == div0$Chr[i] & map$posV3 <= div0$WinCenter[i]) #find index of map anchors smaller than observed position
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==div0$Chr[i])][1]-1,na.rm=T) #take corresponding genetic position

higherIndex <- which(map$chrom == div0$Chr[i] & map$posV3 > div0$WinCenter[i]) #find index of map anchors larger than observed position
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[div0$Chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==div0$Chr[i])][length(which(map$chrom==div0$Chr[i]))]+1,na.rm=T) #take corresponding genetic position

scale <- {div0$WinCenter[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

divCm[i] <- newGen
}

div0$cm <- divCm

### Checkpoint ###
save.image("plotDiversity_Weak_Strong.RData")

#### Remove windows with no data from div0
div0with0 <- div0 #backup
div0 <- div0with0[which(div0with0$nSites>=100),]

### Correct tP by nSites
div0$tP.win <- div0$tP
div0$tP <- div0$tP.win/div0$nSites

### For div0 windows compute distance to nearest mis0 substitution
misDis <- rep(NA,nrow(div0))
misGene <- rep(NA,nrow(div0))
misLoc <- rep(NA,nrow(div0))
rows <- nrow(div0)
for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
misTemp <- mis0[which(mis0$chr==div0$Chr[i]),]
dist <- abs(div0$cm[i]-misTemp$cm) # distance to nearest sub
sub <- which(dist==min(dist))[1]
misDis[i] <- div0$cm[i]-misTemp$cm[sub]
misGene[i] <- misTemp$Gene[sub]
misLoc[i] <- paste(misTemp$chr[sub],misTemp$pos[sub],sep=":")
}

### For div0 windows compute distance to nearest syn0 substitution
synDis <- rep(NA,nrow(div0))
synGene <- rep(NA,nrow(div0))
synLoc <- rep(NA,nrow(div0))
rows <- nrow(div0)
for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
synTemp <- syn0[which(syn0$chr==div0$Chr[i]),]
dist <- abs(div0$cm[i]-synTemp$cm) # distance to nearest sub
sub <- which(dist==min(dist))[1]
synDis[i] <- div0$cm[i]-synTemp$cm[sub]
synGene[i] <- synTemp$Gene[sub]
synLoc[i] <- paste(synTemp$chr[sub],synTemp$pos[sub],sep=":")
}


### Checkpoint ###
save.image("plotDiversity_Weak_Strong.RData")



# # ### For div0 windows compute distance to nearest int0 substitution
# # intDis <- rep(NA,nrow(div0))
# # rows <- nrow(div0)
# # for(i in 1:nrow(div0)){
# # cat( 100*i/rows, "% done", "\r")
# # intTemp <- int0[which(int0$chr==div0$Chr[i]),]
# # dist <- abs(div0$cm[i]-intTemp$cm) # distance to nearest sub
# # sub <- which(dist==min(dist))[1]
# # intDis[i] <- div0$cm[i]-intTemp$cm[sub]
# # }
# # 
# # ### Checkpoint ###
# # #save.image("plotDiversity_Weak_Strong.RData")
# 
# 


###################################################################
### Plot diversity standardized by neutral diversity 2/26/2015 ####
###################################################################

### Standardize by diversity far from genes
pi_maizeN <- 0.007341639 # calculated in GENE_GERP.R
#pi_teoN <- 0.01208945 # calculated in GENE_GERP.R
div0$tP.N <- div0$tP/pi_maizeN
#div0Teo$tP.N <- div0Teo$tP/pi_teoN


### Add synDis, synGene, synLoc, misDis, misGene, misLoc, to div0 objects
div0$synDis <- synDis
div0$synGene <- synGene
div0$synLoc <- synLoc
div0$misDis <-misDis
div0$misGene <-misGene
div0$misLoc <-misLoc



#####################################################################
### High quality Figs ###############################################
#####################################################################
#png("plotDiversity_TvM_Folded2_unNeutralized_June.png",width=8,height=7,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
legend("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(3,3),pch=NA)
text(x=-0.0018,y=1.2,labels="A",cex=4)
#dev.off()

#####################################################################
########## End high quality Figs (6-4-2015)##########################
#####################################################################



### Checkpoint ###
save.image("plotDiversity_Weak_Strong.RData")



################################################################################################
### Perform analysis using only sites near genes that are among the most or least conserved ####
################################################################################################

### Load gene data from Simon -- info for mean GERP score across genes.
load("aveGerp.RData")

### isolate most conserved 10% of genes
conservedGenes <- names(aveGerp)[which(aveGerp >= quantile(aveGerp,0.9))]
unconservedGenes <- names(aveGerp)[which(aveGerp <= quantile(aveGerp,0.1))]

### Trim div0
div0MisCons <- div0[which(div0$misGene %in% conservedGenes),]
div0SynCons <- div0[which(div0$synGene %in% conservedGenes),]

div0MisUncons <- div0[which(div0$misGene %in% unconservedGenes),]
div0SynUncons <- div0[which(div0$synGene %in% unconservedGenes),]



### Maize Conserved Plot
png("plotDiversity_Maize_conserved_June",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),main="Maize diversity around subs \n in GERP conserved regions",cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,"foo",cex.axis=1.2)
mtext("Pairwise diversity",side=4,line=3,cex=1.2)
synLowNeutCons <- loess(div0SynCons$tP.N~div0SynCons$synDis,span=0.01)
lines(synLowNeutCons$x[order(synLowNeutCons$x)],synLowNeutCons$fitted[order(synLowNeutCons$x)],col="darkgray",lwd=4)
misLowNeutCons <- loess(div0MisCons$tP.N~div0MisCons$misDis,span=0.01)
lines(misLowNeutCons$x[order(misLowNeutCons$x)],misLowNeutCons$fitted[order(misLowNeutCons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()


### Maize Unconserved Plot
png("plotDiversity_Maize_unconserved_June",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.02,.02),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),main="Maize diversity around subs \n in GERP Unconserved regions",cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,"foo",cex.axis=1.2)
mtext("Pairwise diversity",side=4,line=3,cex=1.2)
synLowNeutUncons <- loess(div0SynUncons$tP.N~div0SynUncons$synDis,span=0.01)
lines(synLowNeutUncons$x[order(synLowNeutUncons$x)],synLowNeutUncons$fitted[order(synLowNeutUncons$x)],col="darkgray",lwd=4)
misLowNeutUncons <- loess(div0MisUncons$tP.N~div0MisUncons$misDis,span=0.01)
lines(misLowNeutUncons$x[order(misLowNeutUncons$x)],misLowNeutUncons$fitted[order(misLowNeutUncons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()

#######
### Checkpoint ###
save.image("plotDiversity_Weak_Strong.RData")






#####################################################################
### Confidence intervals ###########################################
#####################################################################
x<-seq(-0.2,0.2,length.out=100000)

SynCons <- matrix(NA,nrow=length(x),ncol=100)
MisCons <- matrix(NA,nrow=length(x),ncol=100)
SynUncons <- matrix(NA,nrow=length(x),ncol=100)
MisUncons <- matrix(NA,nrow=length(x),ncol=100)

for(i in 1:100){
  print(i)
  div0MisCons_boot <- div0MisCons[sample(nrow(div0MisCons),replace=T),]
  div0SynCons_boot <- div0SynCons[sample(nrow(div0SynCons),replace=T),]
  div0MisUncons_boot <- div0MisUncons[sample(nrow(div0MisUncons),replace=T),]
  div0SynUncons_boot <- div0SynUncons[sample(nrow(div0SynUncons),replace=T),]
  
  synLowNeutCons_boot <- loess(div0SynCons_boot$tP.N~div0SynCons_boot$synDis,span=0.01)
  misLowNeutCons_boot <- loess(div0MisCons_boot$tP.N~div0MisCons_boot$misDis,span=0.01)
  SynCons[,i] <- predict(synLowNeutCons_boot,x)
  MisCons[,i] <- predict(misLowNeutCons_boot,x)
  
  synLowNeutUncons_boot <- loess(div0SynUncons_boot$tP.N~div0SynUncons_boot$synDis,span=0.01)
  misLowNeutUncons_boot <- loess(div0MisUncons_boot$tP.N~div0MisUncons_boot$misDis,span=0.01)
  SynUncons[,i] <- predict(synLowNeutUncons_boot,x)
  MisUncons[,i] <- predict(misLowNeutUncons_boot,x)
}

SynCons_low <- apply(SynCons,1,quantile,0.025)
SynCons_high <- apply(SynCons,1,quantile,0.975)

MisCons_low <- apply(MisCons,1,quantile,0.025)
MisCons_high <- apply(MisCons,1,quantile,0.975)

SynUncons_low <- apply(SynUncons,1,quantile,0.025)
SynUncons_high <- apply(SynUncons,1,quantile,0.975)

MisUncons_low <- apply(MisUncons,1,quantile,0.025)
MisUncons_high <- apply(MisUncons,1,quantile,0.975)

##################################################################################
### Manuscript Supplemental Figs #################################################
##################################################################################
### Maize conserved fig ###
x<-seq(-0.2,0.2,length.out=100000)
png("plotDiversity_TvM_Conserved_Significance_June.png",width=8,height=7,units="in",res=300,pointsize=12)

par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.1),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowNeutCons <- loess(div0SynCons$tP.N~div0SynCons$synDis,span=0.015)
polygon(c(x,rev(x)),c(SynCons_high,rev(SynCons_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeutCons$x[order(synLowNeutCons$x)],synLowNeutCons$fitted[order(synLowNeutCons$x)],col="darkgray",lwd=2)
#misLowNeutCons <- loess(div0MisCons$tP.N~div0MisCons$misDis,span=0.015)
polygon(c(x,rev(x)),c(MisCons_high,rev(MisCons_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeutCons$x[order(misLowNeutCons$x)],misLowNeutCons$fitted[order(misLowNeutCons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=2,fill=c(rgb(.863,.663,.663,0.5),rgb(1,0,0,0.5)),bty="n",border=c(NA,NA))
text(x=-0.0018,y=1.05,labels="A",cex=4)
par(new=TRUE)
par(fig = c(0.17, 0.42, 0.22, 0.47))
par(mar=c(0,0,0,0),mgp=c(0,0.6,0))
plot(NULL,xlim=c(-.15,.15),xlab="",ylab="",ylim=c(0.4,1.3),cex.lab=1,cex.axis=0.9)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=0.9)
polygon(c(x,rev(x)),c(SynCons_high,rev(SynCons_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeutCons$x[order(synLowNeutCons$x)],synLowNeutCons$fitted[order(synLowNeutCons$x)],col="darkgray",lwd=1)
polygon(c(x,rev(x)),c(MisCons_high,rev(MisCons_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeutCons$x[order(misLowNeutCons$x)],misLowNeutCons$fitted[order(misLowNeutCons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=1)
dev.off()

### Maize unconserved fig ###
x<-seq(-0.2,0.2,length.out=100000)
png("plotDiversity_TvM_Unconserved_Significance_June.png",width=8,height=7,units="in",res=300,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.01,.01),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0,1.5),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011), at=c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011)/pi_maizeN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowNeutUncons <- loess(div0SynUncons$tP.N~div0SynUncons$synDis,span=0.4)
polygon(c(x,rev(x)),c(SynUncons_high,rev(SynUncons_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeutUncons$x[order(synLowNeutUncons$x)],synLowNeutUncons$fitted[order(synLowNeutUncons$x)],col="darkgray",lwd=2)
#misLowNeutUncons <- loess(div0MisUncons$tP.N~div0MisUncons$misDis,span=0.4)
polygon(c(x,rev(x)),c(MisUncons_high,rev(MisUncons_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeutUncons$x[order(misLowNeutUncons$x)],misLowNeutUncons$fitted[order(misLowNeutUncons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=2,fill=c(rgb(.863,.663,.663,0.5),rgb(1,0,0,0.5)),bty="n",border=c(NA,NA))
text(x=-0.01,y=1.45,labels="B",cex=4)
par(new=TRUE)
par(fig = c(0.17, 0.42, 0.22, 0.47))
par(mar=c(0,0,0,0),mgp=c(0,0.6,0))
plot(NULL,xlim=c(-.15,.15),xlab="",ylab="",ylim=c(0.6,1.6),cex.lab=1,cex.axis=0.9)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=0.9)
polygon(c(x,rev(x)),c(SynUncons_high,rev(SynUncons_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeutUncons$x[order(synLowNeutUncons$x)],synLowNeutUncons$fitted[order(synLowNeutUncons$x)],col="darkgray",lwd=1)
polygon(c(x,rev(x)),c(MisUncons_high,rev(MisUncons_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeutUncons$x[order(misLowNeutUncons$x)],misLowNeutUncons$fitted[order(misLowNeutUncons$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=1)
dev.off()



##################################################################################
### End Supplemental Manuscript Figs #############################################
##################################################################################

### Checkpoint ###
save.image("plotDiversity_Weak_Strong.RData")

