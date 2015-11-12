#!/bin/bash

#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e ${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -N msmc
#PBS -q testq
#PBS -l vmem=256Gb,pmem=8Gb,mem=256Gb,nodes=1:ppn=32:ib,walltime=120:00:00
cd $PBS_O_WORKDIR
module use /data004/software/GIF/modules
module load msmc
module load python/3.4.1
module load perl

#cat BKN_chr*.txt | sort -n -k1,1 -n -k2,2 > BKN_allChr.txt
python3 /data004/software/GIF/packages/msmc/20140623/tools/multihetsep_bootstrap.py -n 50 -s 20000000 --chunks_per_chromosome 10 --nr_chromosomes 10 TIL_bootstrap TIL_chr1.txt TIL_chr2.txt TIL_chr3.txt TIL_chr4.txt TIL_chr5.txt TIL_chr6.txt TIL_chr7.txt TIL_chr8.txt TIL_chr9.txt TIL_chr10.txt
python3 /data004/software/GIF/packages/msmc/20140623/tools/multihetsep_bootstrap.py -n 50 -s 20000000 --chunks_per_chromosome 10 --nr_chromosomes 10 BKN_bootstrap BKN_chr1.txt BKN_chr2.txt BKN_chr3.txt BKN_chr4.txt BKN_chr5.txt BKN_chr6.txt BKN_chr7.txt BKN_chr8.txt BKN_chr9.txt BKN_chr10.txt
