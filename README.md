#Scripts used to analyze linked selection and demography in maize

Below is a description of each script used in Beissinger et al (journal and year TBA). Please feel encouraged to contact the Tim or Jeff (corresponding authors) for further descriptions if needed.

## Part 1: SFS, Genotype likelihoods, and SNP calling

**fullScript_genic.sh** Computes SFS, 2D-SFS, and diversity info for maize and teosinte in genic regions only.

**fullScript_intergenic.sh** Computes SFS, 2D-SFS, and diversity info for maize and teosinte in intergenic regions only.

**fullScript_WholeGenome.sh** Computes SFS, 2D-SFS, and diversity info for maize and teosinte, genome-wide.

**callSNPsAngsd_maize_July.sh** Calls SNPs in maize.

**callSNPsAngsd_teo_july.sh** Calls SNPs in teosinte.

**countSingletonsMaize.R** Counts maize singletons in 1000 bp windows

**countSingletonsTeo.R** Counts teosinte singletons in 1000 bp windows

**countSingletonsMaize_downsample.R** Downsamples maize to a sample size equal to teosinte and then counts maize singletons in 1000 bp windows, 

## Part 2: Code to study demography
**ConvertSFStoDadi.R** Script to convert SFS and 2D-SFS from ANGSD output format to dadi input format

**realisticBottleneck.py** Python function that specifies the domestication model dadi will optimize.

**FINAL_DOMESTICATION_MODEL_SCRIPT_Linux.py** Script to run 1,000 iterations of dadi to estimate maize and teosinte demographic model.

**2d_intergenic_TIL_BKN_unpolarized.dadi** Joint SFS for inputing into dadi for above script. (generated with angsd scripts and format converted with R)

**MSMC** A directory (with its own readme) containing scripts to run MSMC


## Part 2: Code to study selection

**findSubstitutions.R** Identifies substitutions between maize and tripsicum, teosinte and tripsicum, teosinte and maize

**exploreSubsitutions.R** Explores substitutions.

**TvM.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. maize substitutions.

**TvT.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. teosinte substitutions.

**TvMT.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. maize/teosinte (both) substitutions.

**H12_array_medWindows.sh** Script to calculate H12 statistic for our data. Windows are of size 200 SNPs with a step of 25 SNPs.

**plotDiversity_TvM.R** This script includes code to look at pairwise diversity surrounding synonymous and non-synonymous substitutions between maize and tripsicum and between teosinte and tripsicum.

**plotDiversity_TvM_Singletons.R** This script includes code to look at singleton diversity surrounding synonymous and non-synonymous substitutions between maize and tripsicum and between teosinte and tripsicum.

**plotDiversity_Weak_Strong.R** This script includes code to look at diversity surrounding synonymous and non-synonymous substitutions between maize and tripsicum at conserved and unconserved sites.

**GENE_GERP.R** This script includes code to look at pairwise diversity surrounding genes and conserved sites (according to GERP).

**GENE_GERP_Singletons.R** This script includes code to look at singleton diversity surrounding genes and conserved sites.

**GENE_GERP_SINGLETONS_Downsampled.R** This script investigates singleton diversity surrounding genes. It includes code to look at singletons from the downsampled set of maize.

**SingletonsVsPi.R** This script simulatneously investigates singleton and pairwise diversity surrounding genes.
