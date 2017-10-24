**Case study: multi-trait mapping of 41 traits in yeast**

The yeast dataset used to show the applicability for multi-trait mapping of a 
large number of traits was first published by Bloom and colleagues in 2013 
[Bloom et al, (2013) "Finding the sources of missing heritability in a yeast
cross" Nature]  and retrieved from
http://genomics-pubs.princeton.edu/YeastCross_BYxRM/.

1. Phenotypes
For analysis with LiMMBo, samples must be fully phenotyped. The missing 
phenotypes in this yeast dataset were imputed when feasible. Before imputation, 
missingness mechanism of data was examined and tested for missingness completely
at random (MCAR) and missingness at random (MAR). The best predictor variables
for each phenotype were chosen based on simulated missingness on the subset of 
fully phenotyped samples. Samples with high missingness rates and phenotypes 
that could not reliably be imputed were excluded from the study. All 
analysis were done in *yeast/phenotypes/phenotypes_yeast.R*, creating Figure
S8-S10 (publication) and Figure 5.1-5.4 (thesis)

1. Genotypes
Genotypes were filtered for samples that passed the phenotype pre-processing
as described above and formated for plink (.ped/.map) via 
*yeast/genotypes/genotypes_yeast.R*.

1. Relationship 
Genetic relationship between the samples in the test cohort was estimated in 
*yeast/genotypes/relationship_yeast.sh* via plink's grm option. Different 
pruning set-ups were tested.

Genotype and phenotype data was formated for the linear mixed model software
limix via *yeast/GWAS/format_data4limix_yeast.sh*.

The variance decomposition of into genetic and noise variance components 
Cg and Cn was done via *runLiMMBo*, the association mapping in 
uni-variate or multi-variate mode providing Cg and Cn via *gwas.py*. 
Parameter set-ups for both function calls can be found in 
*yeast/GWAS/GWAS_yeast.sh*.

The results of the association mapping are analysed and summarised in
plots via *yeast/GWAS/GWAS_analysis.R*, creating Figure 5 (publication) and 
Figure 5.5, 5.6 and B2 (thesis).






