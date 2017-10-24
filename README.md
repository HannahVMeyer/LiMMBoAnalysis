# LiMMBo: Linear Mixed Model with Bootstrapping

This repository contains the python module **limmbo** which enables multi-trait genome-wide association analyses of high-dimensional phenotypes as well as all set-up and analyses scripts for **model evaluation**. 

**Background**
Genome-wide association studies have helped to shed light on the genetic architecture of complex traits and diseases. In recent years, phenotyping of the samples has often gone beyond single traits and it has become common to record multi- to high-dimensional phenotypes for individuals. Whilst these rich datasets offer the potential to analyse complex trait structures and pleiotropic effects at a genome-wide level, multi-trait mapping of high-dimensional phenotypes poses computational challenges which have limited the analyses to a moderate number of traits.  

**LiMMBo**
LiMMBo is a novel and computationally efficient approach for multivariate analysis of high-dimensional phenotypes based on linear mixed models with bootstrapping (LiMMBo). 
1. **LiMMBo is accessible as an open source, Python module**. The implementation, requirements and documentation can be found in [*limmbo*](limmbo).

1. **Model performance**: LiMMBo's performance was tested in terms of covariance paramter estimation, calibration, run time and power based on simulated datasets. The data simulation strategy can be found in [*simulateData*](simulateData). Variance decompostion and genetic association analyses were conducted with files in [*setupLiMMBo*](setupLiMMBo) and evaluated with files from [*EvaluateLiMMBo*](EvaluateLiMMBo). The simulation studies demonstrate the statistical validity of LiMMBo and show that it yields results consistent with existing methods when applied to datasets with a moderate number of traits. In these cases, LiMMBo is faster than existing methods and can scale to hundreds of phenotypes. In a systematic power for genetic association studies on phenotypes with different trait architectures and sizes, multivariate genetic mapping via LiMMBo improves power relative to univariate approaches. The actual increase in power depends on phenotype dimension, pleiotropy and the proportion of phenotypic variance explained by genetics.

1. **Case study: Multi-trait mapping of 41 traits in yeast with LiMMBo**: To demonstrate that LiMMBo enables more powerful analyses than previously possible, we analysed a publicly available dataset with a large number of growth traits in yeast. We detect pleiotropic loci (i.e. loci that influence two or more phenotypic traits) for a number of known and novel correlated traits. The data processing and GWAS can be found in [*yeast*](yeast). 
