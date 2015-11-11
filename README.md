#Scripts used to analyze linked selection and demography in maize

Below is a description of each script used in Beissinger et

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



## Part 2: Code to study selection

**findSubstitutions.R** Identifies substitutions between maize and tripsicum, teosinte and tripsicum, teosinte and maize

**exploreSubsitutions.R** Explores substitutions.

**TvM.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. maize substitutions.

**TvT.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. teosinte substitutions.

**TvMT.vep** Use Ensembl Variant Effect Predictor to estimate effects of tripsicum vs. maize/teosinte (both) substitutions.

**H12_array_medWindows.sh** Script to calculate H12 statistic for our data. Windows are of size 200 SNPs with a step of 25 SNPs.