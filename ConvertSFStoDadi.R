###############################################################################
### This script converts SFS output from angsd to the format readable by dadi #
###############################################################################

### Timothy M. Beissinger
### 8-19-2014

### set working directory
setwd("~/Documents/DomesticationBottleneck/Dom_Bot/WholeGenomeFolded2/SFS")

## ### Convert BKN_genic.sfs
## sfsBKNgenic <- exp(scan("../OUTS/BKN_genic.sfs"))*96805725 #96,805,725 bp included
## sfsBKNgenic <- sfsBKNgenic[seq(1,length(sfsBKNgenic),2)] # treat as sample from haploid
## sfsBKNgenicDadi <- matrix(nrow=2,ncol=1)
## sfsBKNgenicDadi[1,1] <- paste(length(sfsBKNgenic), "unfolded",sep=" ")
## sfsBKNgenicDadi[2,1] <- paste(sfsBKNgenic,collapse=" ")
## write.table(sfsBKNgenicDadi,file="sfsBKN_genic_unpolarized.dadi",quote=F,row.names=F,col.names=F)

## ### Convert BKN_intergenic.sfs
## sfsBKNintergenic <- exp(scan("../OUTS/BKN_intergenic.sfs"))*412340033  #412,340,033 bp included
## sfsBKNintergenic <- sfsBKNintergenic[seq(1,length(sfsBKNintergenic),2)] # treat as sample from haploid
## sfsBKNintergenicDadi <- matrix(nrow=2,ncol=1)
## sfsBKNintergenicDadi[1,1] <- paste(length(sfsBKNintergenic), "unfolded",sep=" ")
## sfsBKNintergenicDadi[2,1] <- paste(sfsBKNintergenic,collapse=" ")
## write.table(sfsBKNintergenicDadi,file="sfsBKN_intergenic_unpolarized.dadi",quote=F,row.names=F,col.names=F)

## ### Convert BKN_WholeGenome.sfs
## sfsBKNwholeGenome <- exp(scan("../OUTS/BKN_WholeGenome.sfs"))*613848587  #613,848,587  bp included
## sfsBKNwholeGenome <- sfsBKNwholeGenome[seq(1,length(sfsBKNwholeGenome),2)] # treat as sample from haploid
## sfsBKNwholeGenomeDadi <- matrix(nrow=2,ncol=1)
## sfsBKNwholeGenomeDadi[1,1] <- paste(length(sfsBKNwholeGenome), "unfolded",sep=" ")
## sfsBKNwholeGenomeDadi[2,1] <- paste(sfsBKNwholeGenome,collapse=" ")
## write.table(sfsBKNwholeGenomeDadi,file="../SFS/sfsBKN_WholeGenome_unpolarized.dadi",quote=F,row.names=F,col.names=F)


## ### Convert TIL_genic.sfs
## sfsTILgenic <- exp(scan("../OUTS/TIL_genic.sfs"))*102809863 #102,809,863 bp included
## sfsTILgenic <- sfsTILgenic[seq(1,length(sfsTILgenic),2)] # treat as sample from haploid
## sfsTILgenicDadi <- matrix(nrow=2,ncol=1)
## sfsTILgenicDadi[1,1] <- paste(length(sfsTILgenic), "unfolded",sep=" ")
## sfsTILgenicDadi[2,1] <- paste(sfsTILgenic,collapse=" ")
## write.table(sfsTILgenicDadi,file="sfsTIL_genic_unpolarized.dadi",quote=F,row.names=F,col.names=F)

## ### Convert TIL_intergenic.sfs
## sfsTILintergenic <- exp(scan("../OUTS/TIL_intergenic.sfs"))*283746710 #283,746,710 bp included
## sfsTILintergenic <- sfsTILintergenic[seq(1,length(sfsTILintergenic),2)] # treat as sample from haploid
## sfsTILintergenicDadi <- matrix(nrow=2,ncol=1)
## sfsTILintergenicDadi[1,1] <- paste(length(sfsTILintergenic), "unfolded",sep=" ")
## sfsTILintergenicDadi[2,1] <- paste(sfsTILintergenic,collapse=" ")
## write.table(sfsTILintergenicDadi,file="sfsTIL_intergenic_unpolarized.dadi",quote=F,row.names=F,col.names=F)

## ### Convert TIL_WholeGenome.sfs
## sfsTILwholeGenome <- exp(scan("../OUTS/TIL_WholeGenome.sfs"))*487289037 #487,289,037 bp included
## sfsTILwholeGenome <- sfsTILwholeGenome[seq(1,length(sfsTILwholeGenome),2)] # treat as sample from haploid
## sfsTILwholeGenomeDadi <- matrix(nrow=2,ncol=1)
## sfsTILwholeGenomeDadi[1,1] <- paste(length(sfsTILwholeGenome), "unfolded",sep=" ")
## sfsTILwholeGenomeDadi[2,1] <- paste(sfsTILwholeGenome,collapse=" ")
## write.table(sfsTILwholeGenomeDadi,file="../SFS/sfsTIL_WholeGenome_unpolarized.dadi",quote=F,row.names=F,col.names=F)


### Convert 2d sfs. Script adapted from Jacob Crawford, https://github.com/mfumagalli/ngsToolsDev/blob/master/convert.2Dsfs.to.dadi.R

## ############################ Genic first ############################
## # Read in 2D sfs genic
## sfsgenic <- exp(read.table("../OUTS/2dsfs_genic.TIL.BKN.sfs"))*87878483 ##87,878,483 bp included
## sfsgenic <- sfsgenic[seq(1,nrow(sfsgenic),2),seq(1,ncol(sfsgenic),2)]

## # Get sample sizes and make header
## n1=nrow(sfsgenic);
## n2=ncol(sfsgenic);
## ns=paste(n1,n2,collapse=' ');
## ns=paste(ns,"unfolded",collapse=' ');

## # Convert 2D sfs to dadi array format
## dadi=NULL;
## for(i in 1:n1){
## 	dadi=c(dadi,as.numeric(sfsgenic[i,]))
## }

## # Write out dadi format to file with same name as 2D sfs file with .fs appeded
## write.table(ns,file="2d_genic_TIL_BKN_unpolarized.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
## write.table(paste(dadi,collapse=' '),file="2d_genic_TIL_BKN_unpolarized.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);

## ##


## #########################################################################
## ############################ Intergenic next ############################
## #########################################################################
## # Read in 2D sfs intergenic
## sfsintergenic <- exp(read.table("../OUTS/2dsfs_intergenic.TIL.BKN.sfs"))*224987968 #224,987,968 bp included
## sfsintergenic <- sfsintergenic[seq(1,nrow(sfsintergenic),2),seq(1,ncol(sfsintergenic),2)]

## # Get sample sizes and make header
## n1=nrow(sfsintergenic);
## n2=ncol(sfsintergenic);
## ns=paste(n1,n2,collapse=' ');
## ns=paste(ns,"unfolded",collapse=' ');

## # Convert 2D sfs to dadi array format
## dadi=NULL;
## for(i in 1:n1){
## 	dadi=c(dadi,as.numeric(sfsintergenic[i,]))
## }

## # Write out dadi format to file with same name as 2D sfs file with .fs appeded
## write.table(ns,file="../SFS/2d_intergenic_TIL_BKN_unpolarized.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
## write.table(paste(dadi,collapse=' '),file="../SFS/2d_intergenic_TIL_BKN_unpolarized.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);


#########################################################################
############################ Whole genome last ##########################
#########################################################################
# Read in 2D sfs WholeGenome
sfsWholeGenome <- exp(read.table("../OUTS/2dsfs_WholeGenome.TIL.BKN.sfs"))*117887368 #117,887,368 bp included
sfsWholeGenome <- sfsWholeGenome[seq(1,nrow(sfsWholeGenome),2),seq(1,ncol(sfsWholeGenome),2)]

# Get sample sizes and make header
n1=nrow(sfsWholeGenome);
n2=ncol(sfsWholeGenome);
ns=paste(n1,n2,collapse=' ');
ns=paste(ns,"unfolded",collapse=' ');

# Convert 2D sfs to dadi array format
dadi=NULL;
for(i in 1:n1){
	dadi=c(dadi,as.numeric(sfsWholeGenome[i,]))
}

# Write out dadi format to file with same name as 2D sfs file with .fs appeded
write.table(ns,file="../SFS/2d_WholeGenome_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE);
write.table(paste(dadi,collapse=' '),file="../SFS/2d_WholeGenome_TIL_BKN.dadi",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE);
