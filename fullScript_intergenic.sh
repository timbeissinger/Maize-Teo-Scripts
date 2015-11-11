#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/WholeGenomeFolded2
#SBATCH -J fullScript_intergenic
#SBATCH -o slurm-log/fullScript_intergenic_%j.out
#SBATCH -p bigmemh
#SBATCH -e slurm-log/fullScript_intergenic_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
set -e
set -u

# This script will compute sfs, thetas, tajima's d, and more for the whole INTERGENIC
# genome of maize and teosinte. First, it computes the UNFOLDED sfs for maize and teosinte, pretending that B73 is the ancestor,
# independently. After that, sites that overlap between the two taxa are identified,
# and the UNFOLDED sfs for each taxa is computed again on just those sites. Next, the 2d sfs
# between maize and teosinte is computed. Finally, thetas and Tajima's D over 1000 bp
# windows are computed.


# Set initial values
pop1=TIL
pop2=BKN

angsdir=/home/beissing/bin/angsd0.612
outputdir=/home/beissing/Dom_Bot_Git/WholeGenomeFolded2/OUTS
ref=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
anc=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
pop1List=INS/"$pop1"_list.txt
pop2List=INS/"$pop2"_list.txt
nIndPop1=$( wc -l "$pop1List" | cut -f 1 -d " " )
nIndPop2=$( wc -l "$pop2List" | cut -f 1 -d " " )
nPop1=$( expr 2 \* $nIndPop1 )
nPop2=$( expr 2 \* $nIndPop2 )
minperc=0.8
minIndPop1=$( printf "%.0f" $(echo "scale=2;$minperc*$nIndPop1" | bc))
minIndPop2=$( printf "%.0f" $(echo "scale=2;$minperc*$nIndPop2" | bc))
pop1F=INS/"$pop1".indF
pop2F=INS/"$pop2".indF
glikehood=1
minMapQ=30
cpu=12
regionfile="/home/beissing/Dom_Bot_Git/WholeGenomeFolded2/INS/intergenicRegionFile.txt"
windowsize=1000
step=1000

# Echo input values
echo "pop1 = "$pop1""
echo "pop2 = "$pop2""
echo "angsdir = "$angsdir""
echo "outputdir = "$outputdir""
echo "ref = "$ref""
echo "anc = "$anc""
echo "pop1List = "$pop1List""
echo "pop2List = "$pop2List""
echo "nIndPop1 = "$nIndPop1""
echo "nIndPop2 = "$nIndPop2""
echo "nPop1 = "$nPop1""
echo "nPop2 = "$nPop2""
echo "minperc = "$minperc""
echo "minIndPop1 = "$minIndPop1""
echo "minIndPop2 = "$minIndPop2""
echo "pop1F = "$pop1F""
echo "pop2F = "$pop2F""
echo "glikehood = "$glikehood""
echo "minMapQ = "$minMapQ""
echo "cpu = "$cpu""
echo "regionfile = "$regionfile""
echo "windowsize = "$windowsize""
echo "step = "$step""

# ######################
# ##### Compute pop1 sfs
# command1="-bam "$pop1List" -out "$outputdir"/"$pop1"_intergenic -doMajorMinor 1 -doMaf 1 -indF "$pop1F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop1 -minInd $minIndPop1 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile"
# echo $command1
# echo 
# $angsdir/angsd $command1

# command2=""$outputdir"/"$pop1"_intergenic.saf $nPop1 -P $cpu" 
# echo $command2
# echo
# $angsdir/misc/realSFS $command2 > "$outputdir"/"$pop1"_intergenic.sfs

# #####################
# #### Compute pop2 sfs
# command3="-bam "$pop2List" -out "$outputdir"/"$pop2"_intergenic -doMajorMinor 1 -doMaf 1 -indF "$pop2F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop2 -minInd $minIndPop2 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile"
# echo $command1
# echo 
# $angsdir/angsd $command3

# command4=""$outputdir"/"$pop2"_intergenic.saf $nPop2 -P $cpu" 
# echo $command4
# echo
# $angsdir/misc/realSFS $command4 > "$outputdir"/"$pop2"_intergenic.sfs

# ####################################
# # Next, extract the compressed files
# command5=""$outputdir"/"$pop1"_intergenic.saf.pos.gz "
# command6=""$outputdir"/"$pop2"_intergenic.saf.pos.gz "
# echo gunzip -k $command5
# echo
# echo gunzip -k $command6
# echo
# gunzip -k $command5
# gunzip -k $command6

####################################################################
# Now we find the positions that occur in both populations using the
# uniq POSIX program
command7=" "$outputdir"/"$pop1"_intergenic.saf.pos "$outputdir"/"$pop2"_intergenic.saf.pos|sort|uniq -d " 
echo cat $command7
echo 
cat "$outputdir"/"$pop1"_intergenic.saf.pos "$outputdir"/"$pop2"_intergenic.saf.pos|sort|uniq -d > "$outputdir"/intersect."$pop1"."$pop2"_intergenic.txt

####################################################################
# Now redo angsd sample allele frequency calculation by conditioning
# on the sites that occur in both populations.
command8=" -bam "$pop1List" -out "$outputdir"/"$pop1"_intergenic_conditioned -doMajorMinor 1 -doMaf 1 -indF "$pop1F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop1 -minInd $minIndPop1 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile -sites "$outputdir"/intersect."$pop1"."$pop2"_intergenic.txt"
echo "$angsdir"/angsd $command8
echo 
"$angsdir"/angsd $command8

command9=" -bam "$pop2List" -out "$outputdir"/"$pop2"_intergenic_conditioned -doMajorMinor 1 -doMaf 1 -indF "$pop2F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop2 -minInd $minIndPop2 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile -sites "$outputdir"/intersect."$pop1"."$pop2"_intergenic.txt"
echo "$angsdir"/angsd $command9
echo 
"$angsdir"/angsd $command9


#########################################################
# Now we estimate the joint sfs using the realSFS program
command10="2dsfs "$outputdir"/"$pop1"_intergenic_conditioned.saf "$outputdir"/"$pop2"_intergenic_conditioned.saf $nPop1 $nPop2 -P $cpu"
echo "$angsdir"/misc/realSFS $command10 
echo 
"$angsdir"/misc/realSFS $command10  > "$outputdir"/2dsfs_intergenic."$pop1"."$pop2".sfs


##########################################
### Calculate thetas for each site in pop1.  TREAT AS UNFOLDED!
command11="-bam "$pop1List" -out "$outputdir"/"$pop1"_intergenic_windows -pest "$outputdir"/"$pop1"_intergenic.sfs -indF "$pop1F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop1 -minInd $minIndPop1 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile -doThetas 1 -doMajorMinor 1 -doMaf 1"
echo $command11
echo
$angsdir/angsd $command11

#####################################################
### Estimate Tajima's D and other statistics for pop1.  TREAT AS UNFOLDED!
command12="make_bed "$outputdir"/"$pop1"_intergenic_windows.thetas.gz "
echo $command12
echo 
$angsdir/misc/thetaStat $command12

command13=" do_stat  "$outputdir"/"$pop1"_intergenic_windows.thetas.gz -nChr $nPop1 -win $windowsize -step $step "
echo $command13
echo
$angsdir/misc/thetaStat $command13


##########################################
### Calculate thetas for each site in pop2.  TREAT AS UNFOLDED!
command14="-bam "$pop2List" -out "$outputdir"/"$pop2"_intergenic_windows -pest "$outputdir"/"$pop2"_intergenic.sfs -indF "$pop2F" -doSaf 2 -uniqueOnly 0 -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop2 -minInd $minIndPop2 -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile -doThetas 1 -doMajorMinor 1 -doMaf 1"
echo $command14
echo
$angsdir/angsd $command14

#####################################################
### Estimate Tajima's D and other statistics for pop2.  TREAT AS UNFOLDED!
command15="make_bed "$outputdir"/"$pop2"_intergenic_windows.thetas.gz "
echo $command15
echo 
$angsdir/misc/thetaStat $command15

command16=" do_stat  "$outputdir"/"$pop2"_intergenic_windows.thetas.gz -nChr $nPop2 -win $windowsize -step $step "
echo $command16
echo
$angsdir/misc/thetaStat $command16