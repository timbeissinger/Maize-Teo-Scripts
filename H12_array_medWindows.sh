#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/H12
#SBATCH -J H12_array_medWindows
#SBATCH -o slurm-log/H12_array_medWindows%A_%a.out
#SBATCH -p bigmemh
#SBATCH -e slurm-log/H12_array_medWindows%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=beissinger@ucdavis.edu
#SBATCH --array=1-10
set -e
set -u

# This script will run  H12 scriptfrom Nandita Garud on maize/teo data. Windows will be medium compared to previous runs (200 SNPs)
# 7/15/2015

python SelectionHapStats/scripts/H12_H2H1.py InputData/maizeSNPs_$SLURM_ARRAY_TASK_ID.txt 23 -o H12_OUT/maizeH12_w200j25_$SLURM_ARRAY_TASK_ID.txt -w 200 -j 25
python SelectionHapStats/scripts/H12peakFinder.py H12_OUT/maizeH12_w200j25_$SLURM_ARRAY_TASK_ID.txt -o H12_OUT/maizeH12_w200j25_peaks_$SLURM_ARRAY_TASK_ID.txt
Rscript SelectionHapStats/scripts/H12_viz.R H12_OUT/maizeH12_w200j25_$SLURM_ARRAY_TASK_ID.txt H12_OUT/maizeH12_w200j25_peaks_$SLURM_ARRAY_TASK_ID.txt H12_OUT/maizeH12_w200j25_$SLURM_ARRAY_TASK_ID.pdf 10



python SelectionHapStats/scripts/H12_H2H1.py InputData/teoSNPs_$SLURM_ARRAY_TASK_ID.txt 13 -o H12_OUT/teoH12_w200j25_$SLURM_ARRAY_TASK_ID.txt -w 50 -j 25
python SelectionHapStats/scripts/H12peakFinder.py H12_OUT/teoH12_w200j25_$SLURM_ARRAY_TASK_ID.txt -o H12_OUT/teoH12_w200j25_peaks_$SLURM_ARRAY_TASK_ID.txt
Rscript SelectionHapStats/scripts/H12_viz.R H12_OUT/teoH12_w200j25_$SLURM_ARRAY_TASK_ID.txt H12_OUT/teoH12_w200j25_peaks_$SLURM_ARRAY_TASK_ID.txt H12_OUT/teoH12_w200j25_$SLURM_ARRAY_TASK_ID.pdf 10
