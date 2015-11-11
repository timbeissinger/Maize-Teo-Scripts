################################################################################
### Use this script to make a plot of SINGLETON diversity surrounding Syn,   ###
### non-syn substitutions between tripsicum and maize.                       ###
################################################################################

### Timothy M. Beissinger
### 1/22/2015

### Set working directory
setwd("~/Documents/DomesticationBottleneck/Dom_Bot/Selection/scripts/")

### Load Trip vs. Maize effects estimates
effects <- read.table("../SNPs/TvMeffects.txt",header=T,stringsAsFactors=F,sep="\t",comment.char="",na.strings="-",skip=8)
levels(as.factor(effects$Consequence))


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
save.image("plotDiversity_TvM_Singletons.RData")

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

### Interpolate genetic position for every int0 SNP
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
int0$cm <- NA

for(i in 1:nrow(int0)){
  print(i)
  lowerIndex <- which(map$chrom == int0$chr[i] & map$posV3 <= int0$pos[i]) #find index of map anchors smaller than observed position
  belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
  belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==int0$chr[i])][1]-1,na.rm=T) #take corresponding genetic position
  
  higherIndex <- which(map$chrom == int0$chr[i] & map$posV3 >= int0$pos[i]) #find index of map anchors larger than observed position
  abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[int0$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
  aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==int0$chr[i])][length(which(map$chrom==int0$chr[i]))]+1,na.rm=T) #take corresponding genetic position
  
  scale <- {int0$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
  newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position
  
  int0$cm[i] <- newGen
}


### CHECKPOINT ###
save.image("plotDiversity_TvM_Singletons.RData")


### Load singleton info
load("../OUTS/tF_maize.Robj")
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

### CHECKPOINT ###
save.image("plotDiversity_TvM_Singletons.RData")


#### Interpolate genetic position for every tF.teo position with info
#chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
#teoCm <- rep(NA,nrow(tF.teo))
#rows <- nrow(tF.teo)
#
#for(i in 1:nrow(tF.teo)){
#  cat( 100*i/rows, "% done", "\r")
#  lowerIndex <- which(map$chrom == tF.teo$chr[i] & map$posV3 <= tF.teo$center[i]) #find index of map anchors smaller than observed position
#  belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
#  belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==tF.teo$chr[i])][1]-1,na.rm=T) #take corresponding genetic position
#  
#  higherIndex <- which(map$chrom == tF.teo$chr[i] & map$posV3 > tF.teo$center[i]) #find index of map anchors larger than observed position
#  abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[tF.teo$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
#  aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==tF.teo$chr[i])][length(which(map$chrom==tF.teo$chr[i]))]+1,na.rm=T) #take corresponding genetic position
#  
#  scale <- {tF.teo$center[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
#  newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position
#  
#  teoCm[i] <- newGen
#}
#
#tF.teo$cm <- teoCm
#
#### CHECKPOINT ###
#save.image("plotDiversity_TvM_Singletons.RData")

### For div0 windows compute distance to nearest mis0 substitution
misDis <- rep(NA,nrow(tF.maize))
rows <- nrow(tF.maize)
for(i in 1:nrow(tF.maize)){
  cat( 100*i/rows, "% done", "\r")
  misTemp <- mis0[which(mis0$chr==tF.maize$chr[i]),]
  dist <- abs(tF.maize$cm[i]-misTemp$cm) # distance to nearest sub
  sub <- which(dist==min(dist))[1]
  misDis[i] <- tF.maize$cm[i]-misTemp$cm[sub]
}


### For div0 windows compute distance to nearest syn0 substitution
synDis <- rep(NA,nrow(tF.maize))
rows <- nrow(tF.maize)
for(i in 1:nrow(tF.maize)){
  cat( 100*i/rows, "% done", "\r")
  synTemp <- syn0[which(syn0$chr==tF.maize$chr[i]),]
  dist <- abs(tF.maize$cm[i]-synTemp$cm) # distance to nearest sub
  sub <- which(dist==min(dist))[1]
  synDis[i] <- tF.maize$cm[i]-synTemp$cm[sub]
}


### For div0 windows compute distance to nearest int0 substitution
intDis <- rep(NA,nrow(tF.maize))
rows <- nrow(tF.maize)
for(i in 1:nrow(tF.maize)){
  cat( 100*i/rows, "% done", "\r")
  intTemp <- int0[which(int0$chr==tF.maize$chr[i]),]
  dist <- abs(tF.maize$cm[i]-intTemp$cm) # distance to nearest sub
  sub <- which(dist==min(dist))[1]
  intDis[i] <- tF.maize$cm[i]-intTemp$cm[sub]
}



### Loess plot
png("plotDiversity_TvM_Singletons.png",width=8,height=6,units="in",res=300)
plot(NULL,xlim=c(-.005,.005),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(0.00,0.025))
#synLow <- loess(tF.maize$tF~synDis,span=0.01)
lines(synLow$x[order(synLow$x)],synLow$fitted[order(synLow$x)],col="darkgray",lwd=3)
#misLow <- loess(tF.maize$tF~misDis,span=0.01)
lines(misLow$x[order(misLow$x)],misLow$fitted[order(misLow$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
#intLow <- loess(tF.maize$tF~intDis,span=0.01)
#lines(intLow$x[order(intLow$x)],intLow$fitted[order(intLow$x)],col="green",lwd=3)
legend("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(3,3,3),pch=NA)
dev.off()

### CHECKPOINT ###
save.image("plotDiversity_TvM_Singletons.RData")


