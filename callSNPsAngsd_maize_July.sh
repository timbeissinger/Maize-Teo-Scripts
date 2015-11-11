#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/WholeGenomeFolded2
#SBATCH -J callSNPsAngsd_maize_july
#SBATCH -o slurm-log/callSNPsAngsd_maize_july%j.out
#SBATCH -p bigmemh
#SBATCH -e slurm-log/callSNPsAngsd_maize_july%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=beissinger@ucdavis.edu
set -e
set -u

# This script will call SNPs in the TILs and BKNs, using Angsd
# 7/6/2015


### Set values
angsdir=/home/beissing/bin/angsd0.610
output=SNPs
ref=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
snpCallList=INS/snpCallList_maize.txt
regionfile="/home/beissing/Dom_Bot_Git/WholeGenome/INS/wholeGenomeRegionFile.txt"
minMapQ=30
minQ=20
nInd=$( wc -l "$snpCallList" | cut -f 1 -d " " )
minperc=0.8
minInd=$( printf "%.0f" $(echo "scale=2;$minperc*$nInd" | bc))
popF=INS/BKN.indF
glikehood=1
minMapQ=30
cpu=6
SNP_pval=1e-6

### Run angsd on bams to call snps

echo "$angsdir/angsd -bam $snpCallList -GL $glikehood -out $output/maize_snps_july -doMaf 1 -indF $popF -doMajorMinor 1 -doGeno 5 -doPost 1 -postCutoff 0.95 -minMapQ $minMapQ -minQ 20 -minInd $minInd -rf $regionfile -P $cpu -SNP_pval $SNP_pval"
echo 
$angsdir/angsd -bam $snpCallList -GL $glikehood -out $output/maize_snps_july -doMaf 1 -indF $popF -doMajorMinor 1 -doGeno 5 -doPost 1 -postCutoff 0.95 -minMapQ $minMapQ -minQ 20 -minInd $minInd -rf $regionfile -P $cpu -SNP_pval $SNP_pval



