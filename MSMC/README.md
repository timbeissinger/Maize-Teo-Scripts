The [msmc](https://github.com/stschiff/msmc) analyses were conducted for maize and teosinte separately. Six samples from each group with high sequencing depth were selected and each sample was treated as a haploid. 

1. Call SNPs on all Hapmpa2 samples with Angsd
```
/data004/software/GIF/packages/ANGSD/0.614/angsd -bam Hapmap2.samples -doMaf 1 -doMajorMinor 1 -uniqueOnly 1 -minMapQ 30 -minQ 20 -minInd 29 -anc /home/lwang/lwang/Zea_mays.AGPv3/Zea_mays.AGPv3.22.dna.genome.fa -GL 1 -doGlf 2 -SNP_pval 1e-6 -indF Hapmap2.samples.F -doGeno 8 -doPost 1 -postCutoff 0.95 -doSaf 2 -P 32 -out Hapmap2 
```

2. Convert Angsd genotype likelihood file into vcf using the script beagleGLF2VCF_Hapmap2.pl

3. selecting the high-depth samples to make fake diploids with the script VCF_phased.pl

4. split the large vcf file into vcfs containing only single chromosome and single individual (the fake diploids) with vcftools. One example as follows:

```
vcftools --vcf ./Hapmap2_landraces/msmcInput/Hapmap2.angsd.phased.vcf --chr 1 --indv TIL09_10 --recode --out ./Hapmap2_landraces/msmcInput/TIL09_10.chr1
```

5. using the msmc built script generate_multihetsep.py to generate the msmc input files. Details of the commands can be found [here](https://github.com/stschiff/msmc-tools).

This step contains three input components:
A. the mask file for each individual and each chromosome, which filters out the too low (less than half of the mean depth) and too high (more than twice of the mean depth) depth nucleotides. Details on generating these files can be found [here](https://github.com/stschiff/msmc-tools).
B. The mappability mask file, which indicates the genomic regions with unique mappability. Details on generating these files can be found [here](https://github.com/HuffordLab/Wang_Private/blob/master/demography/analyses/msmc/mappabilityMask.pdf).
C. the single chromosome single individual vcfs from step 4

6. running msmc 
```
msmc --fixedRecombination -p 20*2+20*4+10*2 -o Hapmap2.TIL.6Hap.p2 TIL_chr1.txt TIL_chr2.txt TIL_chr3.txt TIL_chr4.txt TIL_chr5.txt TIL_chr6.txt TIL_chr7.txt TIL_chr8.txt TIL_chr9.txt TIL_chr10.txt
msmc --fixedRecombination -p 20*2+20*4+10*2 -o Hapmap2.BKN.6Hap.p2 BKN_chr1.txt BKN_chr2.txt BKN_chr3.txt BKN_chr4.txt BKN_chr5.txt BKN_chr6.txt BKN_chr7.txt BKN_chr8.txt BKN_chr9.txt BKN_chr10.txt
```

7. convert the final result file into the real time and Ne with u=3e-8 and generation time as 1 year/generation. (popSize_plot.py)

8. 20 bootstrapping run files were generated with "20150903.msmc_bootstrap.sh" for teosinte and maize, respectively.

9. repeat step 6  and 7 to run msmc on all the 20 bootstrapping files for each group

10. The figure in the final report is generated with R script "msmc.bootstrapping.figures.R".





