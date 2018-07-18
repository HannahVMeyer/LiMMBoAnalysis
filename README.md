# LiMMBo: Linear Mixed Model with Bootstrapping

This repository contains the set-up and analyses scripts for the **model evaluation** of the python module [**limmbo**](https://github.com/HannahVMeyer/limmbo) which enables multi-trait genome-wide association analyses of high-dimensional phenotypes.

## Background
Genome-wide association studies have helped to shed light on the genetic architecture of complex traits and diseases. In recent years, phenotyping of the samples has often gone beyond single traits and it has become common to record multi- to high-dimensional phenotypes for individuals. Whilst these rich datasets offer the potential to analyse complex trait structures and pleiotropic effects at a genome-wide level, multi-trait mapping of high-dimensional phenotypes poses computational challenges which have limited the analyses to a moderate number of traits.  

## LiMMBo
LiMMBo is a novel and computationally efficient approach for multivariate analysis of high-dimensional phenotypes based on linear mixed models with bootstrapping (LiMMBo). 
1. **LiMMBo is accessible as an open source, Python module**. The implementation, requirements and documentation can be found in [*limmbo*](https://github.com/HannahVMeyer/limmbo).

1. **Model performance**: LiMMBo's performance was tested in terms of covariance paramter estimation, calibration, run time and power based on simulated datasets. The data simulation strategy can be found in [*simulateData*](simulateData). Variance decompostion and genetic association analyses were conducted with files in [*setupLiMMBo*](setupLiMMBo) and evaluated with files from [*EvaluateLiMMBo*](EvaluateLiMMBo). The simulation studies demonstrate the statistical validity of LiMMBo and show that it yields results consistent with existing methods when applied to datasets with a moderate number of traits. In these cases, LiMMBo is faster than existing methods and can scale to hundreds of phenotypes. In a systematic power for genetic association studies on phenotypes with different trait architectures and sizes, multivariate genetic mapping via LiMMBo improves power relative to univariate approaches. The actual increase in power depends on phenotype dimension, pleiotropy and the proportion of phenotypic variance explained by genetics.

1. **Case study: Multi-trait mapping of 41 traits in yeast with LiMMBo**: To demonstrate that LiMMBo enables more powerful analyses than previously possible, we analysed a publicly available dataset with a large number of growth traits in yeast. We detect pleiotropic loci (i.e. loci that influence two or more phenotypic traits) for a number of known and novel correlated traits. The data processing and GWAS can be found in [*yeast*](yeast). 

## Publication
LiMMBo and associated analsyes have been described in this [PhD thesis](https://github.com/HannahVMeyer/Thesis) and publication (under review). Figures and Tables generated through code in this repository:
1. Thesis figures: 
  * Figure 4.2: [*EvaluateTime.R*](EvaluateLiMMBo/EvaluateTime.R)
  * Figure 4.3: [*EvaluateCovariance.R*](EvaluateLiMMBo/EvaluateCovariance.R)
  * Figure 4.4, 4.5: [*EvaluateCalibration.R*](EvaluateLiMMBo/EvaluateCalibration.R)
  * Figure 4.6, B7: [*EvaluatePower.R*](EvaluateLiMMBo/EvaluatePower.R)
  * Figure 5.1-5.4: [*phenotypes_yeast.R*](yeast/phenotypes/phenotypes_yeast.R)
  * Figure 5.5, 5.6, B2: [*GWAS_analysis.R*](yeast/GWAS/GWAS_analysis.R)
2. Publication figures:
  * Figure 2: [*EvaluateTime.R*](EvaluateLiMMBo/EvaluateTime.R)
  * Figure 3A: [*EvaluateCovariance.R*](EvaluateLiMMBo/EvaluateCovariance.R)
  * Figure 3B, Table S4: [*EvaluateCalibration.R*](EvaluateLiMMBo/EvaluateCalibration.R)
  * Figure 4, Table S7: [*EvaluatePower.R*](EvaluateLiMMBo/EvaluatePower.R)
  * Figure 5: [*GWAS_analysis.R*](yeast/GWAS/GWAS_analysis.R)
  * Figure S8-S10: [*phenotypes_yeast.R*](yeast/phenotypes/phenotypes_yeast.R)
