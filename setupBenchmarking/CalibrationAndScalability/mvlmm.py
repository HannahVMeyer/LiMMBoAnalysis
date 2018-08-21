import numpy as np
import pandas as pd
import argparse
from pylmm import mvLMM

parser = argparse.ArgumentParser(
    prog="mvLMM",
    description=('Estimates the genetic and non-genetic trait'
                 'covariance matrix parameters of a linear mixed model '
                 'with random genetic and non-genetic effect via a '
                 'bootstrapping-based approach.'))
required = parser.add_argument_group('Basic required arguments')
required.add_argument(
    '-p',
    '--file_pheno',
    action="store",
    dest="file_pheno",
    required=False,
    help=('Path [string] to [(N+1) x (P+1)] .csv file of [P]'
          'phenotypes with [N] samples (first column: sample IDs, first '
          'row: phenotype IDs). Default: %(default)s'))
required.add_argument(
    '--pheno_delim',
    action="store",
    dest="pheno_delim",
    required=False,
    default=",",
    help=('Delimiter of phenotype file. Default: %(default)s'))
required.add_argument(
    '-k',
    '--file_kinship',
    action="store",
    dest="file_relatedness",
    required=False,
    default=None,
    help=('Path [string] to [N x (N+1)] file of kinship/relatedness '
          'matrix with [N] samples (first row: sample IDs). Required when '
          '--lmm/-lm. Default: '
          '%(default)s'))

options = parser.parse_args()
file_Y = options.file_pheno
file_K = options.file_kinship

if debug:
    file_Y='/homes/hannah/data/LiMMBo/simulateData/phenotypes/Calibration/Traits10_samples1000_NrSNP0_Cg0.3_modelnoiseFixedAndBggeneticBgOnly/seed201/Ysim.csv'
    file_K='/homes/hannah/data/LiMMBo/simulateData/genotypes/relatedEU_nopopstructure/N1000/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv'

Y = pd.io.parsers.read_csv(file_Y, index_col=0)
K = pd.io.parsers.read_csv(file_K)

# Instantiate a LMM object for the phentoype Y and fit the null model
L = mvLMM.mvLMM(Y,K)
R = L.getMax()

# The important results are then stored in M.mxCor
Cg= L.mxCor[0]
Cn = L.mxCor[1]
