**Case study: multi-trait mapping of 41 traits in yeast**

The yeast dataset to show the applicability for multi-trait mapping of a 
large number of traits was first published by Bloom and colleagues in 2013 
[Bloom et al, (2013) "Finding the sources of missing heritability in a yeast
cross" Nature]  and retrieved from
http://genomics-pubs.princeton.edu/YeastCross_BYxRM/.

For analyis with LiMMBo/mtSet, samples must be fully phenotyped. The missing 
phenotypes in this yeast dataset were imputed when feasable. Samples with
high missingness rates and phenotypes that could not reliably be imputed were
excluded from the study (in *phenotypes_yeast.R*).

Genotypes were filtered for samples that passed the phenotype pre-processing
as described above and formated for plink via *genotypes_yeast.R*. 

Genetic relationship between the samples in the test cohort was estimated in 
*relationship_yeast.sh* via plink's grm option. Different pruning set-ups 
were tested.

Genotype and phenotype data was formated for the linear mixed model software
limix via *format_data4limix_yeast.sh*.

The variance decomposition of into genetic and noise variance components 
Cg and Cn was done via *runLiMMBo*, the association mapping in 
uni-variate or multi-variate mode providing Cg and Cn via *gwas.py*. 
Parameter set-ups for both function calls can be found in *xxx*.

The results of the association mapping are analysed and summarised in
plots via *GWAS_analysis.R*






