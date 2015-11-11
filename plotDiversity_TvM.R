################################################################################
### Use this script to make a plot of diversity surrounding Syn, non-syn     ###
### substitutions between tripsicum and maize.                               ###
################################################################################

### 9/16/2014
### Updated 12-9-2014

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
save.image("plotDiversity_TvM.RData")

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
save.image("plotDiversity_TvM.RData")


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

#### Remove windows with no data from div0
div0with0 <- div0 #backup
div0 <- div0with0[which(div0with0$nSites>=100),]

### Correct tP by nSites
div0$tP.win <- div0$tP
div0$tP <- div0$tP.win/div0$nSites

### CHECKPOINT ###
save.image("plotDiversity_TvM.RData")
#
# ### Determine p.theta up- and down- stream of each sub in mis0
# divMisRight <- matrix(NA,nrow=nrow(mis0),ncol=200)
# colnames(divMisRight) <- seq(0.001,0.2,0.001)
# divMisLeft <-  matrix(NA,nrow=nrow(mis0),ncol=200)
# colnames(divMisLeft) <- seq(0.001,0.2,0.001)
#
# for(i in 1:nrow(mis0)){
# print(i)
# pos <- mis0$cm[i]
# chr <- mis0$chr[i]
# divTemp <- div0[which(div0$Chr==chr & abs(div0$cm-pos)<= 1 ),]
#
# for(j in 1:200){
# maxRight <- j/1000
# minRight <- j/1000-0.001
# divMisRight[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}<=maxRight & {divTemp$cm-pos}>=minRight)])
#
# maxLeft <- -j/1000
# minLeft <- -j/1000+0.001
# divMisLeft[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}>=maxLeft & {divTemp$cm-pos}<=minLeft)])
# }
# }
#
#
#
# ### Determine p.theta up- and down- stream of each sub in syn0
# divSynRight <- matrix(NA,nrow=nrow(syn0),ncol=200)
# colnames(divSynRight) <- seq(0.001,0.2,0.001)
# divSynLeft <-  matrix(NA,nrow=nrow(syn0),ncol=200)
# colnames(divSynLeft) <- seq(0.001,0.2,0.001)
#
# for(i in 1:nrow(syn0)){
# print(i)
# pos <- syn0$cm[i]
# chr <- syn0$chr[i]
# divTemp <- div0[which(div0$Chr==chr & abs(div0$cm-pos)<= 1 ),]
#
# for(j in 1:200){
# maxRight <- j/1000
# minRight <- j/1000-0.001
# divSynRight[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}<=maxRight & {divTemp$cm-pos}>=minRight)])
#
# maxLeft <- -j/1000
# minLeft <- -j/1000+0.001
# divSynLeft[i,j] <- mean(divTemp$tP[which({divTemp$cm-pos}>=maxLeft & {divTemp$cm-pos}<=minLeft)])
# }
# }
#
#
# ### CHECKPOINT ###
# save.image("plotDiversity_TvM.RData")
#
# ### Run on Slurm and load back from slurm ###
# load("plotDiversity_TvM_Slurm.RData")
#
# ### Make mis plot ###
#
# divMis <- cbind(divMisLeft[,200:1],divMisRight)
# x <- rep(c(seq(-0.2,-0.001,0.001),seq(0.001,0.2,0.001)),nrow(mis0))
# y <- as.vector(t(divMis))
# plot(x,y)
# lines(spline(x,y),col="red",lwd=3)
#
#
# ### Zoom plot, -0.05 cM to 0.05 cM
# divMis <- cbind(divMisLeft[,10:1],divMisRight[,1:10])
# x <- rep(c(seq(-0.01,-0.001,0.001),seq(0.001,0.01,0.001)),nrow(mis0))
# y <- as.vector(t(divMis))
# plot(x,y)
# lines(spline(x,y),col="red",lwd=3)
#
# ### Make syn plot
# divSyn <- cbind(divSynLeft[,200:1],divSynRight)
# x <- rep(c(seq(-0.2,-0.001,0.001),seq(0.001,0.2,0.001)),nrow(syn0))
# y <- as.vector(t(divSyn))
# plot(x,y)
# lines(spline(x,y),col="red",lwd=3)
#
#
# ### SCRATCH ###
# divMis <- cbind(divMisLeft[,100:1],divMisRight)
# x <- rep(c(seq(-1,-0.01,0.01),seq(0.01,1,0.01)),nrow(mis0))
# y <- as.vector(t(divMis))
# plot(x,y)
# lines(spline(x,y),col="red",lwd=3)
#
# divSyn <- cbind(divSynLeft[,10:1],divSynRight[,1:10])
# x <- rep(c(seq(-.1,-0.01,0.01),seq(0.01,.1,0.01)),nrow(syn0))
# y <- as.vector(t(divSyn))
# plot(x,y)
# lines(spline(x,y),col="red",lwd=3)


### For div0 windows compute distance to nearest mis0 substitution
misDis <- rep(NA,nrow(div0))
rows <- nrow(div0)
for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
misTemp <- mis0[which(mis0$chr==div0$Chr[i]),]
dist <- abs(div0$cm[i]-misTemp$cm) # distance to nearest sub
sub <- which(dist==min(dist))[1]
misDis[i] <- div0$cm[i]-misTemp$cm[sub]
}

### Checkpoint ###
save.image("plotDiversity_TvM.RData")

png("../MissenseDiversity.png")
plot(misDis,div0$tP,xlim=c(-1,1),xlab="Distance to nearest missense substitution",ylab="Diversity")
#misSpline <- spline(misDis,div0$tP)
misSpline <- smooth.spline(misDis,div0$tP)
lines(misSpline$x,misSpline$y,col="red",lwd=3)
dev.off()

#scatter.smooth(misDis,div0$tP,span=0.001,col="red",xlim=c(0,1))

### For div0 windows compute distance to nearest syn0 substitution
synDis <- rep(NA,nrow(div0))
rows <- nrow(div0)
for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
synTemp <- syn0[which(syn0$chr==div0$Chr[i]),]
dist <- abs(div0$cm[i]-synTemp$cm) # distance to nearest sub
sub <- which(dist==min(dist))[1]
synDis[i] <- div0$cm[i]-synTemp$cm[sub]
}

### Checkpoint ###
save.image("plotDiversity_TvM.RData")

png("../SynonymousDiversity.png")
plot(synDis,div0$tP,xlim=c(-.03,.03),xlab="Distance to nearest synonymous substitution",ylab="Diversity")
#synSpline <- spline(synDis,div0$tP)
synSpline <- smooth.spline(synDis,div0$tP)
lines(synSpline$x,synSpline$y,col="blue",lwd=3)
misSpline <- smooth.spline(misDis,div0$tP)
lines(misSpline$x,misSpline$y,col="red",lwd=3)
dev.off()



### For div0 windows compute distance to nearest int0 substitution
intDis <- rep(NA,nrow(div0))
rows <- nrow(div0)
for(i in 1:nrow(div0)){
cat( 100*i/rows, "% done", "\r")
intTemp <- int0[which(int0$chr==div0$Chr[i]),]
dist <- abs(div0$cm[i]-intTemp$cm) # distance to nearest sub
sub <- which(dist==min(dist))[1]
intDis[i] <- div0$cm[i]-intTemp$cm[sub]
}

### Checkpoint ###
save.image("plotDiversity_TvM.RData")


### Loess
png("plotDiversity_TvM_including_noncoding.png",width=10,height=8,units="in",res=300)
#plot(intDis,div0$tP,xlim=c(-.005,.005),xlab="Distance to nearest intergenic substitution",ylab="Diversity",ylim=c(0,0.015))
#smoothScatter(intDis,div0$tP,xlim=c(-.005,.005),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(0,0.015),col="red")
plot(NULL,xlim=c(-.005,.005),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(0,0.015))
#synLow <- loess(div0$tP~synDis,span=0.03)
lines(synLow$x[order(synLow$x)],synLow$fitted[order(synLow$x)],col="darkgray",lwd=3)
#misLow <- loess(div0$tP~misDis,span=0.03)
lines(misLow$x[order(misLow$x)],misLow$fitted[order(misLow$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
#intLow <- loess(div0$tP~intDis,span=0.03)
lines(intLow$x[order(intLow$x)],intLow$fitted[order(intLow$x)],col="green",lwd=3)
legend("bottomright","(x,y)", c("Synonymous","Missense","Intergenic","*Note: background shading indicates the \n density of the underlying data"),col=c("darkgray","darkred","green",NA),lwd=c(3,3,3,NA),pch=NA)
dev.off()


png("plotDiversity_TvM_including_noncoding_Folded2.png",width=8,height=6,units="in",res=300)
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(0.003,0.011))
synLow <- loess(div0$tP~synDis,span=0.01)
lines(synLow$x[order(synLow$x)],synLow$fitted[order(synLow$x)],col="darkgray",lwd=3)
misLow <- loess(div0$tP~misDis,span=0.01)
lines(misLow$x[order(misLow$x)],misLow$fitted[order(misLow$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
intLow <- loess(div0$tP~intDis,span=0.01)
lines(intLow$x[order(intLow$x)],intLow$fitted[order(intLow$x)],col="green",lwd=3)
legend("bottomright","(x,y)", c("Synonymous","Missense","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
dev.off()



################
### Scale to mean 0, variance 1, and plot (for comparison with teosinte)
div0$tP.normalized <- (div0$tP-mean(div0$tP,na.rm=T))/sd(div0$tP)

png("plotDiversity_TvM_including_noncoding_NORMALIZED.png",width=10,height=8,units="in",res=400)
#smoothScatter(intDis,div0$tP.normalized,xlim=c(-.005,.005),xlab="Distance to nearest intergenic substitution",ylab="Diversity",ylim=c(-1,1),col="red")
plot(NULL,xlim=c(-.005,.005),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(-1,1))
#synLowNorm <- loess(div0$tP.normalized~synDis,span=0.03)
lines(synLowNorm$x[order(synLowNorm$x)],synLowNorm$fitted[order(synLowNorm$x)],col="darkgray",lwd=3)
#misLowNorm <- loess(div0$tP.normalized~misDis,span=0.03)
lines(misLowNorm$x[order(misLowNorm$x)],misLowNorm$fitted[order(misLowNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
#intLowNorm <- loess(div0$tP.normalized~intDis,span=0.03)
lines(intLowNorm$x[order(intLowNorm$x)],intLowNorm$fitted[order(intLowNorm$x)],col="green",lwd=3)
legend("top","(x,y)", c("Synonymous","Missense","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
dev.off()

################
### Compare normalized diversity to the same in teosinte
load("TeoDiversityObjects.Robj")

png("plotDiversity_TvT_vs_TvM_including_noncoding_NORMALIZED_Folded2.png",height=1000,width=1000)
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution",ylab="Diversity",ylim=c(-0.8,0.5))
#synLowNorm <- loess(div0$tP.normalized~synDis,span=0.01)
#synLowTeoNorm <- loess(div0Teo$tP.normalized~synDisTeo,span=0.01)
lines(synLowNorm$x[order(synLowNorm$x)],synLowNorm$fitted[order(synLowNorm$x)],col="darkgray",lwd=4)
lines(synLowTeoNorm$x[order(synLowTeoNorm$x)],synLowTeoNorm$fitted[order(synLowTeoNorm$x)],col="darkgray", lwd=2)
#misLowNorm <- loess(div0$tP.normalized~misDis,span=0.01)
#misLowTeoNorm <- loess(div0Teo$tP.normalized~misDisTeo,span=0.01)
lines(misLowNorm$x[order(misLowNorm$x)],misLowNorm$fitted[order(misLowNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
lines(misLowTeoNorm$x[order(misLowTeoNorm$x)],misLowTeoNorm$fitted[order(misLowTeoNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=2)
#intLowNorm <- loess(div0$tP.normalized~intDis,span=0.01)
#intLowTeoNorm <- loess(div0Teo$tP.normalized~intDisTeo,span=0.01)
lines(intLowNorm$x[order(intLowNorm$x)],intLowNorm$fitted[order(intLowNorm$x)],col="green",lwd=4)
lines(intLowTeoNorm$x[order(intLowTeoNorm$x)],intLowTeoNorm$fitted[order(intLowTeoNorm$x)],col="green",lwd=2)
legend("topright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)
dev.off()

png("plotDiversity_TvT_vs_TvM_including_noncoding_NORMALIZED_panels_Folded2.png",height=1000,width=1000,pointsize=20)
par(mfrow=c(2,2))
# All maize
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Standardized Diversity",ylim=c(-0.5,0.7),main="All subs, maize only")
lines(synLowNorm$x[order(synLowNorm$x)],synLowNorm$fitted[order(synLowNorm$x)],col="darkgray",lwd=3)
lines(misLowNorm$x[order(misLowNorm$x)],misLowNorm$fitted[order(misLowNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
lines(intLowNorm$x[order(intLowNorm$x)],intLowNorm$fitted[order(intLowNorm$x)],col="green",lwd=3)
legend("top","(x,y)", c("Synonymous","Nonsynonymous","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
#maize syn vs teo syn
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Standardized Diversity",ylim=c(-0.8,0.5),main="Synonymous substitutions")
lines(synLowTeoNorm$x[order(synLowTeoNorm$x)],synLowTeoNorm$fitted[order(synLowTeoNorm$x)],col="pink", lwd=3)
lines(synLowNorm$x[order(synLowNorm$x)],synLowNorm$fitted[order(synLowNorm$x)],col="purple",lwd=3)
legend("top","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
#maize missense vs teo missense
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Standardized Diversity",ylim=c(-0.8,0.5),main = "Nonsynonymous Substitutions")
lines(misLowTeoNorm$x[order(misLowTeoNorm$x)],misLowTeoNorm$fitted[order(misLowTeoNorm$x)],col="pink", lwd=3)
lines(misLowNorm$x[order(misLowNorm$x)],misLowNorm$fitted[order(misLowNorm$x)],col="purple" ,lwd=3)
legend("top","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
#maize non-coding vs teo noncoding
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Standardized Diversity",ylim=c(-0.8,0.6),main = "Intergenic Substitutions")
lines(intLowNorm$x[order(intLowNorm$x)],intLowNorm$fitted[order(intLowNorm$x)],col="purple",lwd=3)
lines(intLowTeoNorm$x[order(intLowTeoNorm$x)],intLowTeoNorm$fitted[order(intLowTeoNorm$x)],col="pink",lwd=3)
legend("top","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
dev.off()

###################################################################
### Plot diversity standardized by neutral diversity 2/26/2015 ####
###################################################################

### Standardize by diversity far from genes
pi_maizeN <- 0.007341639 # calculated in GENE_GERP.R
pi_teoN <- 0.01208945 # calculated in GENE_GERP.R
div0$tP.N <- div0$tP/pi_maizeN
div0Teo$tP.N <- div0Teo$tP/pi_teoN


png("plotDiversity_TvM_Folded2_Neutralized.png",width=10,height=8,units="in",res=400)
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution",ylab="Diversity / Neutral Diversity",ylim=c(0.5,1.1))
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=3)
#legend("bottomright","(x,y)", c("Synonymous","Missense","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
legend("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(3,3),pch=NA)
dev.off()

png("plotDiversity_TvM_Folded2_unNeutralized.png",width=8,height=7,units="in",res=400)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.5,1.3))
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN)
mtext("Pairwise Diversity",side=4,line=3)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=3)
#legend("bottomright","(x,y)", c("Synonymous","Missense","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
legend("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(3,3),pch=NA)
dev.off()


png("plotDiversity_TvT_vs_TvM_including_noncoding_Neutralized_Folded2.png",height=8,width=10,units="in",res=400)
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.1))
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=4)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=2,lty=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=2,lty=3)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#intLowTeoNeut <- loess(div0Teo$tP.N~intDisTeo,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=4)
#lines(intLowTeoNeut$x[order(intLowTeoNeut$x)],intLowTeoNeut$fitted[order(intLowTeoNeut$x)],col="green",lwd=2)
#legend("bottomright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)
legend("bottomright","(x,y)", c("Synonymous, maize", "Nonsynonymous, maize","Synonymous, teo","Nonsynonymous, teo"),col=c("darkgray","darkred","darkgray","darkred"),lwd=c(4,4,2,2),lty=c(1,1,3,3),pch=NA)
dev.off()

png("plotDiversity_TvT_vs_TvM_including_noncoding_unNeutralized_Folded2.png",height=8,width=10,units="in",res=400)
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution",ylab="Maize Pairwise Diversity",ylim=c(0.4,1.1),yaxt="n")
axis(2, labels=c(0.003,0.004,0.005,0.006,0.007,0.008), at=c(0.003,0.004,0.005,0.006,0.007,0.008)/pi_maizeN)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013)/pi_teoN)
abline(h=1,lty=2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=4)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=2,lty=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=2,lty=3)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#intLowTeoNeut <- loess(div0Teo$tP.N~intDisTeo,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=4)
#lines(intLowTeoNeut$x[order(intLowTeoNeut$x)],intLowTeoNeut$fitted[order(intLowTeoNeut$x)],col="green",lwd=2)
#legend("bottomright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)
legend("bottomright","(x,y)", c("Synonymous, maize", "Nonsynonymous, maize","Synonymous, teo","Nonsynonymous, teo"),col=c("darkgray","darkred","darkgray","darkred"),lwd=c(4,4,2,2),lty=c(1,1,3,3),pch=NA)
dev.off()

##################################################
### Maize Only ###################################
##################################################

png("plotDiversity_TvM_unNeutralized_Folded2_May.png",height=8,width=10,units="in",res=400)
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution",ylab="Maize Pairwise Diversity",ylim=c(0.4,1.1),yaxt="n")
axis(2, labels=c(0.003,0.004,0.005,0.006,0.007,0.008), at=c(0.003,0.004,0.005,0.006,0.007,0.008)/pi_maizeN)
#axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013)/pi_teoN)
#abline(h=1,lty=2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=4)
#lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=2,lty=3)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
#lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=2,lty=3)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#intLowTeoNeut <- loess(div0Teo$tP.N~intDisTeo,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=4)
#lines(intLowTeoNeut$x[order(intLowTeoNeut$x)],intLowTeoNeut$fitted[order(intLowTeoNeut$x)],col="green",lwd=2)
#legend("bottomright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()


##################################################
### Teosintee Only ###############################
##################################################

png("plotDiversity_TvT_unNeutralized_Folded2_May.png",height=8,width=10,units="in",res=400)
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution",ylab="Teosinte Pairwise Diversity",ylim=c(0.4,1.1),yaxt="n")
#axis(2, labels=c(0.003,0.004,0.005,0.006,0.007,0.008), at=c(0.003,0.004,0.005,0.006,0.007,0.008)/pi_maizeN)
axis(2, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013)/pi_teoN)
#abline(h=1,lty=2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
#lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=4)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=4,lty=1)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
#lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=4,lty=1)
#intLowNeut <- loess(div0$tP.N~intDis,span=0.01)
#intLowTeoNeut <- loess(div0Teo$tP.N~intDisTeo,span=0.01)
#lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=4)
#lines(intLowTeoNeut$x[order(intLowTeoNeut$x)],intLowTeoNeut$fitted[order(intLowTeoNeut$x)],col="green",lwd=2)
#legend("bottomright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()



#####################################################################
### High quality Figs ###############################################
#####################################################################
png("plotDiversity_TvM_Folded2_unNeutralized_June.png",width=8,height=7,units="in",res=400,pointsize=12)
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
dev.off()

png("plotDiversity_TvT_Folded2_unNeutralized_June.png",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), at=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015)/pi_teoN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=4,lty=1)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkblue", alpha.f = 0.8), lwd=4,lty=1)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkblue"),lwd=c(4,4),lty=c(1,1),pch=NA)
text(x=-0.0018,y=1.2,labels="B",cex=4)
dev.off()

#####################################################################
########## End high quality Figs (6-4-2015)##########################
#####################################################################








png("plotDiversity_TvT_vs_TvM_panels_Folded2_Neutralized.png",height=1000,width=1000,pointsize=20)
par(mfrow=c(2,2))
# All maize
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(0.4,1.35),main="All subs, maize only")
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=3)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=3)
lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="green",lwd=3)
legend("bottomright","(x,y)", c("Synonymous","Nonsynonymous","Intergenic"),col=c("darkgray","darkred","green"),lwd=c(3,3,3),pch=NA)
#maize syn vs teo syn
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(.4,1.35),main="Synonymous substitutions")
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="pink", lwd=3)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="purple",lwd=3)
legend("bottomright","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
#maize missense vs teo missense
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(.4,1.35),main = "Nonsynonymous Substitutions")
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col="pink", lwd=3)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col="purple" ,lwd=3)
legend("bottomright","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
#maize non-coding vs teo noncoding
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(.4,1.35),main = "Intergenic Substitutions")
lines(intLowNeut$x[order(intLowNeut$x)],intLowNeut$fitted[order(intLowNeut$x)],col="purple",lwd=3)
lines(intLowTeoNeut$x[order(intLowTeoNeut$x)],intLowTeoNeut$fitted[order(intLowTeoNeut$x)],col="pink",lwd=3)
legend("bottomright","(x,y)", c("maize","teosinte"),col=c("purple","pink"),lwd=c(3,3,3),pch=NA)
dev.off()

png("plotDiversity_TvT_vs_TvM_TwoPanels_Folded2_Neutralized.png",height=500,width=1000,pointsize=20)
par(mfrow=c(1,2))
#maize syn vs teo syn
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(.4,1.35),main="Synonymous substitutions")
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="blue", lwd=2)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="red",lwd=2)
legend("bottomright","(x,y)", c("maize","teosinte"),col=c("red","blue"),lwd=c(3,3,3),pch=NA)
#maize missense vs teo missense
plot(NULL,xlim=c(-.003,.003),xlab="Distance to nearest substitution (cM)",ylab="Diversity/Neutral Diversity",ylim=c(.4,1.35),main = "Nonsynonymous Substitutions")
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col="blue", lwd=2)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col="red" ,lwd=2)
legend("bottomright","(x,y)", c("maize","teosinte"),col=c("red","blue"),lwd=c(3,3,3),pch=NA)
dev.off()

### Checkpoint ###
save.image("plotDiversity_TvM_Neutralized.RData")


### Checkpoint ###
save.image("plotDiversity_TvM.RData")
save(div0,file="div0_Folded2.Robj")



#########################################################################################
### Perform analysis using only sites within Hufford candidate domestication regions ####
#########################################################################################

### 6-4-2015

### Load workspace
load("plotDiversity_TvM_Neutralized.RData")

### Load huff's table
candidates <- read.csv("../huffordDomesticationCandidates.csv",header=T,stringsAsFactors=F,sep=",")

### Clean huff's table
candidates <- candidates[,c(3,4,5,6)]
names(candidates) <- c("chr","start","end","XPCLRscore")
trimmed <- candidates[1,]
for(i in 2:nrow(candidates)){
  if (candidates$start[i] != candidates$start[i-1]) trimmed <- rbind(trimmed,candidates[i,])
  
}

### Add synDis, synDisTeo, misDis, misDisTeo, intDis, intDisTeo, to div0 objects
div0$synDis <- synDis
div0Teo$synDis <- synDisTeo
div0$misDis <-misDis
div0Teo$misDis <- misDisTeo
div0$intDis <- intDis
div0Teo$intDis <- intDisTeo

### Determine sites within huff's regions
div0Trimmed <- div0[0,] #maize
div0TeoTrimmed<-div0Teo[0,] #Teo

for (i in 1:nrow(trimmed)){
  print(i)
  tmpRows <- which(div0$Chr == trimmed$chr[i] & div0$WinCenter >= trimmed$start[i] & div0$WinCenter <= trimmed$end[i])
  div0Trimmed <- rbind(div0Trimmed,div0[tmpRows,])
} #maize

for (i in 1:nrow(trimmed)){
  print(i)
  tmpRows <- which(div0Teo$Chr == trimmed$chr[i] & div0Teo$WinCenter >= trimmed$start[i] & div0Teo$WinCenter <= trimmed$end[i])
  div0TeoTrimmed <- rbind(div0TeoTrimmed,div0Teo[tmpRows,])
} #teo

### Maize Plot
png("plotDiversity_Maize_candidates_June",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.2,.2),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),main="Maize diversity in Hufford domestication regions",cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,"foo",cex.axis=1.2)
mtext("Pairwise diversity",side=4,line=3,cex=1.2)
#synLowNeutTrimmed <- loess(div0Trimmed$tP.N~div0Trimmed$synDis,span=0.2)
lines(synLowNeutTrimmed$x[order(synLowNeutTrimmed$x)],synLowNeutTrimmed$fitted[order(synLowNeutTrimmed$x)],col="darkgray",lwd=4)
#misLowNeutTrimmed <- loess(div0Trimmed$tP.N~div0Trimmed$misDis,span=0.2)
lines(misLowNeutTrimmed$x[order(misLowNeutTrimmed$x)],misLowNeutTrimmed$fitted[order(misLowNeutTrimmed$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=4)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()


### Teo Plot
png("plotDiversity_Teosinte_candidates_June.png",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.2,.2),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),main="Teosinte diversity in Hufford domestication regions",cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), at=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015)/pi_teoN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowTeoNeutTrimmed <- loess(div0TeoTrimmed$tP.N~div0TeoTrimmed$synDis,span=0.2)
lines(synLowTeoNeutTrimmed$x[order(synLowTeoNeutTrimmed$x)],synLowTeoNeutTrimmed$fitted[order(synLowTeoNeutTrimmed$x)],col="darkgray",lwd=4)
#misLowTeoNeutTrimmed <- loess(div0TeoTrimmed$tP.N~div0TeoTrimmed$misDis,span=0.2)
lines(misLowTeoNeutTrimmed$x[order(misLowTeoNeutTrimmed$x)],misLowTeoNeutTrimmed$fitted[order(misLowTeoNeutTrimmed$x)],col=adjustcolor("darkblue", alpha.f = 0.8) ,lwd=4)
legend("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkblue"),lwd=c(4,4),lty=c(1,1),pch=NA)
dev.off()

#######
### Checkpoint ###
save.image("plotDiversity_TvM_Neutralized.RData")






#####################################################################
### Confidence intervals ###########################################
#####################################################################
x<-seq(-0.2,0.2,length.out=100000)

maizeSyn <- matrix(NA,nrow=length(x),ncol=100)
maizeMis <- matrix(NA,nrow=length(x),ncol=100)
teoSyn <- matrix(NA,nrow=length(x),ncol=100)
teoMis <- matrix(NA,nrow=length(x),ncol=100)
maizeTrimSyn <- matrix(NA,nrow=length(x),ncol=100)
maizeTrimMis <- matrix(NA,nrow=length(x),ncol=100)
teoTrimSyn <- matrix(NA,nrow=length(x),ncol=100)
teoTrimMis <- matrix(NA,nrow=length(x),ncol=100)

for(i in 1:100){
  print(i)
  div0_boot <- div0[sample(nrow(div0),replace=T),]
  div0Teo_boot <- div0Teo[sample(nrow(div0Teo),replace=T),]
  div0Trimmed_boot <- div0Trimmed[sample(nrow(div0Trimmed),replace=T),]
  div0TeoTrimmed_boot <- div0TeoTrimmed[sample(nrow(div0TeoTrimmed),replace=T),]
  
  synLowNeut_boot <- loess(div0_boot$tP.N~div0_boot$synDis,span=0.01)
  misLowNeut_boot <- loess(div0_boot$tP.N~div0_boot$misDis,span=0.01)
  maizeSyn[,i] <- predict(synLowNeut_boot,x)
  maizeMis[,i] <- predict(misLowNeut_boot,x)
  
  synLowTeoNeut_boot <- loess(div0Teo_boot$tP.N~div0Teo_boot$synDis,span=0.01)
  misLowTeoNeut_boot <- loess(div0Teo_boot$tP.N~div0Teo_boot$misDis,span=0.01) 
  teoSyn[,i] <- predict(synLowTeoNeut_boot,x)
  teoMis[,i] <- predict(misLowTeoNeut_boot,x)
  
  synLowNeutTrimmed_boot <- loess(div0Trimmed_boot$tP.N~div0Trimmed_boot$synDis,span=0.2)
  misLowNeutTrimmed_boot <- loess(div0Trimmed_boot$tP.N~div0Trimmed_boot$misDis,span=0.2)
  maizeTrimSyn[,i] <- predict(synLowNeutTrimmed_boot,x)
  maizeTrimMis[,i] <- predict(misLowNeutTrimmed_boot,x)
  
  synLowTeoNeutTrimmed_boot <- loess(div0TeoTrimmed_boot$tP.N~div0TeoTrimmed_boot$synDis,span=0.2)
  misLowTeoNeutTrimmed_boot <- loess(div0TeoTrimmed_boot$tP.N~div0TeoTrimmed_boot$misDis,span=0.2)  
  teoTrimSyn[,i] <- predict(synLowTeoNeutTrimmed_boot,x)
  teoTrimMis[,i] <- predict(misLowTeoNeutTrimmed_boot,x)
}

maizeSyn_low <- apply(maizeSyn,1,quantile,0.025)
maizeSyn_high <- apply(maizeSyn,1,quantile,0.975)

maizeMis_low <- apply(maizeMis,1,quantile,0.025)
maizeMis_high <- apply(maizeMis,1,quantile,0.975)

teoSyn_low <- apply(teoSyn,1,quantile,0.025)
teoSyn_high <- apply(teoSyn,1,quantile,0.975)

teoMis_low <- apply(teoMis,1,quantile,0.025)
teoMis_high <- apply(teoMis,1,quantile,0.975)

maizeTrimSyn_low <- apply(maizeTrimSyn,1,quantile,0.025)
maizeTrimSyn_high <- apply(maizeTrimSyn,1,quantile,0.975)

maizeTrimMis_low <- apply(maizeTrimMis,1,quantile,0.025)
maizeTrimMis_high <- apply(maizeTrimMis,1,quantile,0.975)

teoTrimSyn_low <- apply(teoTrimSyn,1,quantile,0.025)
teoTrimSyn_high <- apply(teoTrimSyn,1,quantile,0.975)

teoTrimMis_low <- apply(teoTrimMis,1,quantile,0.025)
teoTrimMis_high <- apply(teoTrimMis,1,quantile,0.975)

#####################################################################
### Manuscript Figs #################################################
#####################################################################
### Maize fig ###
x<-seq(-0.2,0.2,length.out=100000)
png("plotDiversity_TvM_Folded2_Significance_Aug.png",width=8,height=7,units="in",res=300,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.1),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
polygon(c(x,rev(x)),c(maizeSyn_high,rev(maizeSyn_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=2)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
polygon(c(x,rev(x)),c(maizeMis_high,rev(maizeMis_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous","Nonsynonymous"),col=c("darkgray","darkred"),lwd=2,fill=c(rgb(.863,.663,.663,0.5),rgb(1,0,0,0.5)),bty="n",border=c(NA,NA))
text(x=-0.0018,y=1.08,labels="A",cex=4)
par(new=TRUE)
par(fig = c(0.17, 0.42, 0.22, 0.47))
par(mar=c(0,0,0,0),mgp=c(0,0.6,0))
plot(NULL,xlim=c(-.15,.15),xlab="",ylab="",ylim=c(0.4,1.3),cex.lab=1,cex.axis=0.9)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,cex.axis=0.9)
polygon(c(x,rev(x)),c(maizeSyn_high,rev(maizeSyn_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=1)
polygon(c(x,rev(x)),c(maizeMis_high,rev(maizeMis_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=1)
dev.off()

### Teo fig ###
png("plotDiversity_TvT_Folded2_Significance_Aug.png",height=7,width=8,units="in",res=400,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.1),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), at=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015)/pi_teoN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowTeoNeut <- loess(div0Teo$tP.N~synDisTeo,span=0.01)
polygon(c(x,rev(x)),c(teoSyn_high,rev(teoSyn_low)),col=rgb(.663,.663,.863,0.5),border=NA)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray", lwd=2,lty=1)
#misLowTeoNeut <- loess(div0Teo$tP.N~misDisTeo,span=0.01)
polygon(c(x,rev(x)),c(teoMis_high,rev(teoMis_low)),col=rgb(0,0,1,0.5),border=NA)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkblue", alpha.f = 0.8), lwd=2,lty=1)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkblue"),lwd=2,fill=c(rgb(.663,.663,.863,0.5),rgb(0,0,1,0.5)),bty="n",border=c(NA,NA))
text(x=-0.0018,y=1.08,labels="B",cex=4)
par(new=TRUE)
par(fig = c(0.17, 0.39, 0.22, 0.44))
par(mar=c(0,0,0,0),mgp=c(0,0.6,0))
plot(NULL,xlim=c(-.15,.15),xlab="",ylab="",ylim=c(0.4,1.3),cex.lab=1,cex.axis=0.9)
axis(4, labels=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), at=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015)/pi_teoN,cex.axis=0.9)
polygon(c(x,rev(x)),c(teoSyn_high,rev(teoSyn_low)),col=rgb(.663,.663,.863,0.5),border=NA)
lines(synLowTeoNeut$x[order(synLowTeoNeut$x)],synLowTeoNeut$fitted[order(synLowTeoNeut$x)],col="darkgray",lwd=1)
polygon(c(x,rev(x)),c(teoMis_high,rev(teoMis_low)),col=rgb(0,0,1,0.5),border=NA)
lines(misLowTeoNeut$x[order(misLowTeoNeut$x)],misLowTeoNeut$fitted[order(misLowTeoNeut$x)],col=adjustcolor("darkblue", alpha.f = 0.8) ,lwd=1)
dev.off()

### Hufford Regions Maize Fig ##
png("plotDiversity_TvM_Folded2_Candidates_Significance_June.png",width=8,height=7,units="in",res=300,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.15,.15),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009), at=c(0.003,0.004,0.005,0.006,0.007,0.008,0.009)/pi_maizeN,"foo",cex.axis=1.2)
mtext("Pairwise diversity",side=4,line=3,cex=1.2)
polygon(c(x,rev(x)),c(maizeTrimSyn_high,rev(maizeTrimSyn_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeutTrimmed$x[order(synLowNeutTrimmed$x)],synLowNeutTrimmed$fitted[order(synLowNeutTrimmed$x)],col="darkgray",lwd=2)
polygon(c(x,rev(x)),c(maizeTrimMis_high,rev(maizeTrimMis_low)),col=rgb(1,0,0,0.5),border=NA)
lines(misLowNeutTrimmed$x[order(misLowNeutTrimmed$x)],misLowNeutTrimmed$fitted[order(misLowNeutTrimmed$x)],col="darkred" ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkred"),lwd=2,fill=c(rgb(.863,.663,.663,0.5),rgb(1,0,0,0.5)),bty="n",border=c(NA,NA))
dev.off()

### Hufford Regions Teo Fig ##
png("plotDiversity_TvT_Folded2_Candidates_Significance_June.png",width=8,height=7,units="in",res=300,pointsize=12)
par(mar=c(5,4,4,5))
plot(NULL,xlim=c(-.15,.15),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(0.4,1.3),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), at=c(0.003,0.005,0.007,0.009,0.011,0.013,0.015)/pi_teoN,cex.axis=1.2)
mtext("Pairwise diversity",side=4,line=3,cex=1.2)
polygon(c(x,rev(x)),c(teoTrimSyn_high,rev(teoTrimSyn_low)),col=rgb(.663,.663,.863,0.5),border=NA)
lines(synLowTeoNeutTrimmed$x[order(synLowTeoNeutTrimmed$x)],synLowTeoNeutTrimmed$fitted[order(synLowTeoNeutTrimmed$x)],col="darkgray",lwd=2)
polygon(c(x,rev(x)),c(teoTrimMis_high,rev(teoTrimMis_low)),col=rgb(0,0,1,0.5),border=NA)
lines(misLowTeoNeutTrimmed$x[order(misLowTeoNeutTrimmed$x)],misLowTeoNeutTrimmed$fitted[order(misLowTeoNeutTrimmed$x)],col="darkblue" ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous", "Nonsynonymous"),col=c("darkgray","darkblue"),lwd=2,fill=c(rgb(.663,.663,.863,0.5),rgb(0,0,1,0.5)),bty="n",border=c(NA,NA))
dev.off()

#####################################################################
### End Manuscript Figs #############################################
#####################################################################

### Checkpoint ###
save.image("plotDiversity_TvM_Neutralized.RData")


#####################################################################
### Investigate tga1    #############################################
#####################################################################

tga1_chr <- 4
tga1_start <- 44534815#44509046
tga1_end <- 44539478#44512898
tga1_sub <- tga1_start+270

div0_tga1_region <- div0[which(div0$Chr == tga1_chr & div0$WinCenter >= {tga1_start-50000} & div0$WinCenter <= {tga1_end+50000}),]
plot(div0_tga1_region$WinCenter,div0_tga1_region$tP)
lines(c(tga1_start,tga1_end),c(0,0),col="purple",lwd=10)
tga1_spline<-smooth.spline(div0_tga1_region$tP~div0_tga1_region$WinCenter)
lines(tga1_spline)
points({tga1_start+270},0,col="red",pch=16)


### Interpolate genetic position for tga1_start
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
int0$cm <- NA
  lowerIndex <- which(map$chrom == tga1_chr & map$posV3 <= tga1_sub) #find index of map anchors smaller than observed position
  belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
  belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==tga1_chr)][1]-1,na.rm=T) #take corresponding genetic position
  
  higherIndex <- which(map$chrom == tga1_chr & map$posV3 >= tga1_sub) #find index of map anchors larger than observed position
  abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[tga1_chr],na.rm=T) #take smallest position of anchor that is larger than observed position
  aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==tga1_chr)][length(which(map$chrom==tga1_chr))]+1,na.rm=T) #take corresponding genetic position
  
  scale <- {tga1_sub-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
  tga1_sub_gen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

### For div0_tga1_region windows compute distance to tga1 substitution
div0_tga1_region$tga1Dis <- NA
rows <- nrow(div0_tga1_region)
for(i in 1:nrow(div0_tga1_region)){
  cat( 100*i/rows, "% done", "\r")
  dist <- div0_tga1_region$cm[i]-tga1_sub_gen # distance to nearest sub
  div0_tga1_region$tga1Dis[i] <- dist
}

### Make plot
png("plotDiversity_TvM_Folded2_Significance_tga1Supp_June.png",width=8,height=7,units="in",res=300,pointsize=12)
par(mar=c(5,4,4,5))
x<-seq(-0.2,0.2,length.out=100000)
plot(NULL,xlim=c(-.002,.002),xlab="Distance to nearest substitution (cM)",ylab="Diversity / Neutral Diversity",ylim=c(-0.05,1.5),cex.lab=1.2,cex.axis=1.2)
axis(4, labels=c(0,0.002,0.004,0.006,0.008,0.010,0.012), at=c(0,0.002,0.004,0.006,0.008,0.01,0.012)/pi_maizeN,cex.axis=1.2)
mtext("Pairwise Diversity",side=4,line=3,cex=1.2)
#synLowNeut <- loess(div0$tP.N~synDis,span=0.01)
polygon(c(x,rev(x)),c(maizeSyn_high,rev(maizeSyn_low)),col=rgb(.863,.663,.663,0.5),border=NA)
lines(synLowNeut$x[order(synLowNeut$x)],synLowNeut$fitted[order(synLowNeut$x)],col="darkgray",lwd=2)
#misLowNeut <- loess(div0$tP.N~misDis,span=0.01)
#polygon(c(x,rev(x)),c(maizeMis_high,rev(maizeMis_low)),col=rgb(1,0,0,0.5),border=NA)
#lines(misLowNeut$x[order(misLowNeut$x)],misLowNeut$fitted[order(misLowNeut$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=2)
source("legend_plotDiv.R")
legend.v2("bottomright","(x,y)", c("Synonymous","tga1"),col=c("darkgray","darkred"),lwd=2,fill=c(rgb(.863,.663,.663,0.5),NA),bty="n",border=c(NA,NA))
#text(x=-0.0018,y=1.2,labels="A",cex=4)
#plot(div0_tga1_region$tga1Dis,div0_tga1_region$tP.N,xlim=c(-0.002,0.002))
#points(div0_tga1_region$tga1Dis,div0_tga1_region$tP.N)
tga1_loess<-loess(div0_tga1_region$tP.N~div0_tga1_region$tga1Dis,span=0.15)
#x<-seq(-0.002,0.002,length.out=1000)
tga1_predict<-predict(tga1_loess,x)
lines(x,tga1_predict,col="red",lwd=2)
points({tga1_start+270},0,col="red",pch=16)
dev.off()



save.image("plotDiversity_TvM_Neutralized.RData")

#####################################################################################
### As of June 2015, I don't recall that the lines below here do anything useful ####
#####################################################################################

### Binning ###
bins <- seq(-5,5,1e-3)
bins <- cbind(bins,NA,NA,NA,NA)
bins <- data.frame(bins)
names(bins) <- c("binCenter","misNumber","mis","synNumber","syn")
bins$binCenter <- bins$binCenter+0.0005

for(i in 1:nrow(bins)){
  print(i)
  misInd <- which(misDis>={bins$binCenter[i]-0.0005} & misDis<{bins$binCenter[i]+0.0005})
  synInd <- which(synDis>={bins$binCenter[i]-0.0005} & synDis<{bins$binCenter[i]+0.0005})
  intInd <- which(intDis>={bins$binCenter[i]-0.0005} & intDis<{bins$binCenter[i]+0.0005})
  bins$misNumber[i] <- length(misInd)
  bins$synNumber[i] <- length(synInd)
  bins$intNumber[i] <- length(intInd)
  bins$mis[i] <- mean(div0$tP[misInd],na.rm=T)
  bins$syn[i] <- mean(div0$tP[synInd],na.rm=T)
  bins$int[i] <- mean(div0$tP[intInd],na.rm=T)
}


### Compute loess and plot reduction in diversity
pdf("SubDiversity_TvM.pdf")
misLow <- loess(bins$mis[which(bins$misNumber>=10)]~bins$binCenter[which(bins$misNumber>=10)],span=0.1,weights=bins$misNumber[which(bins$misNumber>=10)])
synLow <- loess(bins$syn[which(bins$synNumber>=10)]~bins$binCenter[which(bins$synNumber>=10)],span=0.1,weights=bins$synNumber[which(bins$synNumber>=10)])
intLow <- loess(bins$int[which(bins$intNumber>=10)]~bins$binCenter[which(bins$intNumber>=10)],span=0.1,weights=bins$intNumber[which(bins$intNumber>=10)])
plot(misLow,col="red",cex=0.5,main="Binned diversity",xlab="Distance from substitution",ylab="Diversity",xlim=c(-0.5,0.5))
points(synLow,col="blue",cex=0.5)
points(intLow,col="darkgreen",cex=0.5,pch=19)
lines(misLow$x,misLow$fitted,col="red",lwd=3)
lines(synLow$x,synLow$fitted,col="blue",lwd=3)
lines(intLow$x,intLow$fitted,col="green",lwd=3)
legend("top","(x,y)",c("Synonymous","Missense","Intergenic"),col=c("blue","red","Green"),lwd=c(3,3),pch=NA)
dev.off()

### Plot zoomed reduction in diversity
misLowZoom <- loess(bins$mis[which(bins$misNumber>=10 & abs(bins$binCenter)<=0.1)]~bins$binCenter[which(bins$misNumber>=10 & abs(bins$binCenter)<=0.1)],span=0.2)
synLowZoom <- loess(bins$syn[which(bins$synNumber>=10 & abs(bins$binCenter)<=0.1)]~bins$binCenter[which(bins$synNumber>=10 & abs(bins$binCenter)<=0.1)],span=0.2)
plot(misLowZoom)
points(synLowZoom)
lines(misLowZoom$x,misLowZoom$fitted,col="red",lwd=3)
lines(synLowZoom$x,synLowZoom$fitted,col="blue",lwd=3)

### Checkpoint ###
save.image("plotDiversity_TvM.RData")


### Export binned maize data
binsMaize <- bins
save(binsMaize,file="binsMaize.Robj")




### SCRATCH ###

# NOTES: THERE ARE NOT ENOUGH SYNONYMOUS SUBSITUTIONS FOR THIS MEASURE TO BE USEFUL AS AN INDICATOR OF DIVERGENCE. PERHAPS IT WOULD WORK WITH AN OUTGROUP.

### Determine divergence as number of synonymous subs in each window over number of bp
divergence <- rep(NA,nrow(div0))
for(i in 1:nrow(div0)){
print(i)
index1<-which(syn0$chr==div0$Chr[i] & syn0$pos >= {div0$WinCenter[i]-50000} & syn0$pos < {div0$WinCenter[i]+50000})
index2<- which(div0$Chr==div0$Chr[i] & div0$WinCenter >= {div0$WinCenter[i]-50000} & div0$WinCenter < {div0$WinCenter[i]+50000})
divergence[i] <- length(index1)/sum(div0$nSites[index2])
}

div0$tP.scaled <- div0$tP/divergence

### BINS ###

for(i in 1:nrow(bins)){
  print(i)
  misInd <- which(misDis>={bins$binCenter[i]-0.0005} & misDis<{bins$binCenter[i]+0.0005})
  synInd <- which(synDis>={bins$binCenter[i]-0.0005} & synDis<{bins$binCenter[i]+0.0005})
  bins$misNumber[i] <- length(misInd)
  bins$synNumber[i] <- length(synInd)
  bins$mis[i] <- mean(div0$tP[misInd],na.rm=T)
  bins$syn[i] <- mean(div0$tP[synInd],na.rm=T)
}








points(misDis,div0$tP.normalized,xlim=c(-0.25,0.25),xlab="Distance to nearest substitution",ylab="Diversity")
#synLowNorm <- loess(div0$tP.normalized~synDis,span=0.03)
plot(misDisTeo,div0Teo$tP.normalized,xlim=c(-0.25,0.25),xlab="Distance to nearest substitution",ylab="Diversity",col="red")
#synLowTeoNorm <- loess(div0Teo$tP.normalized~synDisTeo,span=0.03)
lines(synLowNorm$x[order(synLowNorm$x)],synLowNorm$fitted[order(synLowNorm$x)],col="darkgray",lwd=1)
lines(synLowTeoNorm$x[order(synLowTeoNorm$x)],synLowTeoNorm$fitted[order(synLowTeoNorm$x)],col="darkgray", lwd=1)
#misLowNorm <- loess(div0$tP.normalized~misDis,span=0.03)
#misLowTeoNorm <- loess(div0Teo$tP.normalized~misDisTeo,span=0.03)
lines(misLowNorm$x[order(misLowNorm$x)],misLowNorm$fitted[order(misLowNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8) ,lwd=1)
lines(misLowTeoNorm$x[order(misLowTeoNorm$x)],misLowTeoNorm$fitted[order(misLowTeoNorm$x)],col=adjustcolor("darkred", alpha.f = 0.8), lwd=1)
#intLowNorm <- loess(div0$tP.normalized~intDis,span=0.03)
#intLowTeoNorm <- loess(div0Teo$tP.normalized~intDisTeo,span=0.03)
lines(intLowNorm$x[order(intLowNorm$x)],intLowNorm$fitted[order(intLowNorm$x)],col="brown",lwd=1)
lines(intLowTeoNorm$x[order(intLowTeoNorm$x)],intLowTeoNorm$fitted[order(intLowTeoNorm$x)],col="green",lwd=1)
legend("topright","(x,y)", c("Synonymous, maize", "Synonymous, teo","Missense, maize","Missense, teo","Intergenic, maize","Intergenic, teo"),col=c("darkgray","darkgray","darkred","darkred","green","green"),lwd=c(4,2,4,2,4,2),pch=NA)


