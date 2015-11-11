################################################################################
### This script will load angsd mafs output and count singletons             ###
### It should be run on slurm because mafs files are several gigs.           ###
################################################################################

### Tim Beissinger
### 12-11-2014

### Load data
library(data.table)
mafsMaize <- fread("zcat ../WholeGenomeFolded2/OUTS/BKN_WholeGenome_windows.mafs.gz",header=T,stringsAsFactors=F, showProgress=T)

### Remove sites with few individuals
mafsMaize <- mafsMaize[which(mafsMaize$nInd >= 0.8*23),] #Only include sites with 80% of individuals sequenced

### Split on chromosomes
maizeList <- list()
for(i in 1:10){
    print(i)
    maizeList[[i]] <- mafsMaize[which(mafsMaize$chromo == i),]
}

### Make storage list
maizeResults <- list()
for(i in 1:10){
    print(i)
    chr <- rep(i,ceiling(max(maizeList[[i]]$position)/1000))
    start <- seq(1,max(maizeList[[i]]$position),1000)
    end <-  c(seq(1000,max(maizeList[[i]]$position),1000),max(maizeList[[i]]$position))
    center <- (start+end)/2
    tF <- rep(NA,length(center))
    nSites <- rep(NA,length(center))
    maizeResults[[i]] <- data.frame(cbind(chr,start,end,center,tF,nSites),stringsAsFactors=F)
}

rm(i)

###Analyze

Singletons <- function(x){length(which(x < 3/46 & x >= 1/46))} # maizemax = 23

for(chr in 1:10){
    print(chr)
    cuts <- cut(maizeList[[chr]]$position,breaks=c(seq(1,max(maizeList[[chr]]$position),1000),max(maizeList[[chr]]$position)),right=F)
    split <- split(maizeList[[chr]]$knownEM,cuts)
    maizeResults[[chr]]$nSites <- unlist(sapply(split,length)) 
    maizeResults[[chr]]$tF <- unlist(sapply(split,Singletons))
    }

### Combine maize results
allChromosomes <- rbind(maizeResults[[1]],maizeResults[[2]],maizeResults[[3]],maizeResults[[4]],maizeResults[[5]],maizeResults[[6]],maizeResults[[7]],maizeResults[[8]],maizeResults[[9]],maizeResults[[10]])

### Write results
write.table(allChromosomes,file="tF_maize.txt")
save(allChromosomes,file="tF_maize.Robj")
