**Performance analyses of LiMMBo**

1. **Run time**: Compare empircal run timesfor variance decomposition with REML and LiMMBo via [*EvaluateTime.R*](EvaluateTime.R). Data for analysis generated with [*setupLiMMBo/runtime.sh*](setupLiMMBo/runtime.sh). Generates Figure 2 (publication) and Figure 4.2 (thesis).

1. **Covariance estimation**: Compare the covariance estimates of REML and LiMMBo to true covariance matrices (known from simulation) via [*EvaluateCovariance.R*](EvaluateCovariance.R). Data for analysis generated via [*setupLiMMBo/runtime.sh*](setupLiMMBo/runtime.sh). Generates Figure 3A (publication) and Figure 4.3 (thesis).

1. **Calibation**: Compare the calibration of genetic association studies with either simple LM mor ultivariate LMMs with REML and LiMMBo via [*EvaluateCalibration.R*](EvaluateCalibration.R). Data for analysis generated via [*setupLiMMBo/calibration_vd.sh*](setupLiMMBo/calibration_vd.sh) and [*setupLiMMBo/calibration_association.sh*](setupLiMMBo/calibration_association.sh). Generates Figure 3B, Table S4 (publication) and Figure 4.4, 4.5 (thesis).

1. **Power**: Compare power of genetic association studies for univariate LMM and multivariate LMMs with LiMMBo via [*EvaluatePower.R*](EvaluatePower.R). Data generated via [*setupLiMMBo/power_vd.sh*](setupLiMMBo/power_vd.sh) and [*setupLiMMBo/power_association.sh*](setupLiMMBo/power_association.sh). Generates Figure 4, S7 (publication) and Figure 4.6, B7 (thesis).

