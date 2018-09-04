import random

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog="chromosomeSampling",
    description=("Read snp count for sampled chromosomes"))

parser.add_argument(
        '-pf', '--prefix',
        action="store",
        dest='prefix',
        required=True,
        type=str,
        help=('Chromosome file prefix. Default: %(default)s'))

parser.add_argument(
        '-sf', '--suffix',
        action="store",
        dest='suffix',
        required=True,
        type=str,
        help=('Chromosome file suffix. Default: %(default)s'))

parser.add_argument(
        '-s', '--seed',
        action="store",
        dest='seed',
        required=True,
        type=int,
        help=('Seed for initiation of random number generator.'
        'Default: %(default)s'))

parser.add_argument(
        '-c', '--nrchr',
        action="store",
        dest='nrchr',
        required=False,
        type=int,
        default=5,
        help=('Number of chromosomes. Default: %(default)s'))

parser.add_argument(
        '-o', '--ofile',
        action="store",
        dest='ofile',
        required=True,
        type=str,
        help=('Output file name. Default: %(default)s'))

options = parser.parse_args()
seed = options.seed
prefix =  options.prefix
suffix = options.suffix
oname = options.ofile
nrChr = options.nrchr

random.seed(seed)
chr=random.sample(range(1,23), nrChr)
snps = []

for c in chr:
    fname = prefix + 'chr' + str(c) + suffix
    with open(fname) as f:
        content = f.readlines()
        snps.append(content[0].strip())

with open(oname, 'w') as f:
    f.write(",".join(map(str, chr)))
    f.write("\n")
    f.write(",".join(snps))
