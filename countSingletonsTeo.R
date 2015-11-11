################################################################################
### This script will load angsd mafs output and count singletons             ###
### It should be run on slurm because mafs files are several gigs.           ###
################################################################################

### Tim Beissinger
### 12-11-2014

### Load data
library(data.table)
mafsTeo <- fread("zcat ../WholeGenomeFolded2/OUTS/TIL_WholeGenome_windows.mafs.gz",header=T,stringsAsFactors=F,showProgress=T)


### Remove sites with few individuals
mafsTeo <- mafsTeo[which(mafsTeo$nInd >= 0.8*13),] #Only include sites with 80% of individuals sequenced


### Split on chromosomes
teoList <- list()
for(i in 1:10){
    print(i)
    teoList[[i]] <- mafsTeo[which(mafsTeo$chromo == i),]
}

### Make storage list
teoResults <- list()
for(i in 1:10){
    print(i)
    chr <- rep(i,ceiling(max(teoList[[i]]$position)/1000))
    start <- seq(1,max(teoList[[i]]$position),1000)
    end <-  c(seq(1000,max(teoList[[i]]$position),1000),max(teoList[[i]]$position))
    center <- (start+end)/2
    tF <- rep(NA,length(center))
    nSites <- rep(NA,length(center))
    teoResults[[i]] <- data.frame(cbind(chr,start,end,center,tF,nSites),stringsAsFactors=F)
}

rm(i)

###Analyze

Singletons <- function(x){length(which(x < 3/26 & x >= 1/26))} #teomax = 13

for(chr in 1:10){
    print(chr)
    cuts <- cut(teoList[[chr]]$position,breaks=c(seq(1,max(teoList[[chr]]$position),1000),max(teoList[[chr]]$position)),right=F)
    split <- split(teoList[[chr]]$knownEM,cuts)
    teoResults[[chr]]$nSites <- unlist(sapply(split,length)) 
    teoResults[[chr]]$tF <- unlist(sapply(split,Singletons))
    }

### Combine teo results
allChromosomes <- rbind(teoResults[[1]],teoResults[[2]],teoResults[[3]],teoResults[[4]],teoResults[[5]],teoResults[[6]],teoResults[[7]],teoResults[[8]],teoResults[[9]],teoResults[[10]])

### Write results
write.table(allChromosomes,file="tF_teo.txt")
save(allChromosomes,file="tF_teo.Robj")
