##############################################################################
### Use this script to identify substitutions between trip, teo, and maize ###
### from the SNPs called by angsd.                                         ###
##############################################################################

### Timothy M. Beissinger
### 8-15-2014

### Set working directory
getwd()
setwd("../SNPs")

### Load SNP calls
rawSNPs <- read.table("fooSNPs.txt",header=F,stringsAsFactors=F,sep="\t",na.strings="NN")
SNPs <- rawSNPs[,1:41] #remove extra column

### Load and format header
head <- read.table("../INS/snpCallList.txt",stringsAsFactors=F) #load
head <- unlist(strsplit(head[,1],split="mergedBams/"))[seq(2,74,2)] #seperate directory
head <- unlist(strsplit(head,split="_merged.bam")) #remove tail
head <- c("chr","pos","Major","Minor",head) #add early columns before taxa

### Add header
names(SNPs) <- head

### Function to identify substitutions
whichAreSubs <- function(data,outCols=5,popCols=6:41,missing.limit=0.5){
    if(length(which(is.na(data[popCols]))) > (missing.limit*length(popCols)) ) return ("population NA Issue")
    if(length(which(is.na(data[outCols]))) > (missing.limit*length(outCols)) ) return ("outgroup NA Issue")

    out <- unique(unlist(strsplit(as.character(data[outCols]),split="")))
    out <- out[is.na(out)==F]
    pop <- unique(unlist(strsplit(as.character(data[popCols]),split="")))
    pop <- pop[is.na(pop)==F]
    if(length(out)>0 & length(pop)>0 & length(intersect(out,pop))==0) return ("Sub")
}

### Index by millions
length <- nrow(SNPs) # count number of SNPs
floor <- floor(length/1000000) # find lowest million
starts <-c(1, c(1:floor)*1000000 +1)
stops <- c(starts[1:{length(starts)-1}]+999999,length)

# Trip vs. Maize & Teo
TvMT_result <- c()
for(c in 1:length(starts)){
    print(paste(starts[c],":",stops[c],sep=""))
    SNPsub <- SNPs[starts[c]:stops[c],]
    TvMT_result[starts[c]:stops[c]] <- apply(SNPsub,1,whichAreSubs,outCols=5,popCols=6:41,missing.limit=0.2) # trip vs mays
}
#TvMT_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=6:41,missing.limit=0.2) # trip vs mays
TvMT <- SNPs[which(as.character(TvMT_result)=="Sub"),c(1:4,5,6)]
write.table(TvMT,file="TvMT.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Trip vs. Maize
TvM_result <- c()
for(c in 1:length(starts)){
    print(paste(starts[c],":",stops[c],sep=""))
    SNPsub <- SNPs[starts[c]:stops[c],]
    TvM_result[starts[c]:stops[c]] <- apply(SNPsub,1,whichAreSubs,outCols=5,popCols=19:41,missing.limit=0.2) # trip vs maize
}
#TvM_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=19:41,missing.limit=0.2) # trip vs maize
TvM <- SNPs[which(as.character(TvM_result)=="Sub"),c(1:4,5,19)]
write.table(TvM,file="TvM.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Trip vs. Teo
TvT_result <- c()
for(c in 1:length(starts)){
    print(paste(starts[c],":",stops[c],sep=""))
    SNPsub <- SNPs[starts[c]:stops[c],]
    TvT_result[starts[c]:stops[c]] <- apply(SNPsub,1,whichAreSubs,outCols=5,popCols=6:18,missing.limit=0.2) # trip vs teo
}
#TvT_result<- apply(SNPs,1,whichAreSubs,outCols=5,popCols=6:18,missing.limit=0.2) # trip vs teo
TvT <- SNPs[which(as.character(TvT_result)=="Sub"),c(1:4,5,6)]
write.table(TvT,file="TvT.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Teo vs. Maize
TevM_result <- c()
for(c in 1:length(starts)){
    print(paste(starts[c],":",stops[c],sep=""))
    SNPsub <- SNPs[starts[c]:stops[c],]
    TevM_result[starts[c]:stops[c]] <- apply(SNPsub,1,whichAreSubs,outCols=6:18,popCols=19:41,missing.limit=0.2) # teo vs maize
}
#TevM_result <- apply(SNPs,1,whichAreSubs,outCols=6:18,popCols=19:41,missing.limit=0.2) # teo vs maize
TevM <- SNPs[which(as.character(TevM_result)=="Sub"),c(1:4,6,19)]
write.table(TevM,file="TevM.txt",row.names=F,col.names=T,quote=F,sep="\t")

# Print some summary stats
cat(c("Trip-Maize/Teo Subs = ",nrow(TvMT), "\n"))
cat(c("Trip-Maize Subs = ",nrow(TvM), "\n"))
cat(c("Trip-Teo Subs = ",nrow(TvT), "\n"))
cat(c("Teo-Maize Subs = ",nrow(TevM), "\n"))

# Save image
save.image("SubstitutionSpace.RData")


### Identify investigated SNPs in TvMT
TvMT_tested <- SNPs$pos[which(as.character(TvMT_result) == "NULL" | as.character(TvMT_result) =="Sub" )]

### Determine expected distribution of distance between substitutions
uni_subs <- sample(TvMT_tested,nrow(TvMT))
uni_subs <- uni_subs[order(uni_subs)]
uni_gap <- uni_subs-c(1,uni_subs[1:{length(uni_subs)-1}])

### Determine observed distribution of distance between substitutions
gap <- TvMT$pos-c(1,TvMT$pos[1:{length(TvMT$pos)-1}])

### Density plots--observed vs expected

pdf("Observed_vs_expected.pdf")
par(mfrow=c(2,1))
hist(TvMT$pos,breaks=100,xlab="Substitution positions",main="Observed distribution of substitutions at the start of chromosome 1")
hist(uni_subs,breaks=100,xlab="Substitution positions",main="Expected distribution of substitutions at the start of chromosome 1")
par(mfrow=c(1,1))
plot(density(gap),col="red",lwd=3,main="Density of substitution spacing",xlab="Space between subs",xlim=c(0,50000),ylim=c(0,0.003))
lines(density(uni_gap),col="blue",lwd=3)
legend("topright","(x,y)",c("Observed distance between substitutions","Expected distribution between substitutions"),lwd=3,col=c("red","blue"))
dev.off()


### Load large TvMT
TvMT <- read.table("TvMT.txt",header=T,stringsAsFactors=F)
TvMT_1 <- TvMT[which(TvMT$chr==1),]

hist(TvMT_1$pos,breaks=10000)


gap <- TvMT_1$pos-c(1,TvMT_1$pos[1:{length(TvMT_1$pos)-1}])

uni_samp <- sample(max(TvMT_1$pos),nrow(TvMT_1))[1:100]
uni_samp <- uni_samp[order(uni_samp)]
uni_gap <- uni_samp-c(1,uni_samp[1:{length(uni_samp)-1}])

plot(density(uni_gap,bw=500000),col="blue",lwd=3,xlim=c(0,1.5e7))
lines(density(gap,bw=500000),col="red",lwd=3,main="Density of substitution spacing",xlim=c(0,1.5e7))
legend("topright","(x,y)",c("Observed distance between substitutions","Expected distribution between substitutions"),lwd=3,col=c("red","blue"))
