import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog="GemmaVD-parser", description=("Parses "
    "Gemma log.txt files for Vg and Ve estimates under the null model and for
    time taken for variance estimation."))

parser.add_argument(
        '-f', '--file_log',
        action="store",
        dest='f',
        required=True,
        type=str,
        help=('Path to Gemma log.txt file. Default: %(default)s'))

parser.add_argument(
        '-o', '--outdir',
        action="store",
        dest='o',
        required=True,
        type=str,
        help=('Directory where Vg/Ve estimates and VD time will be saved. '
            'Default: %(default)s'))

vg_lines = False
vg_header = False
ve_lines = False
ve_header = False
vg_list = []
ve_list = []

options = parser.parse_args()
fname = options.f
directory = options.o

#fname = '/homes/hannah/tmp/phenotypeSimulator/GWAS_gemma_LMM_mt_chr21.log.txt'
with open(fname) as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    for line in content:
        if line.startswith('## REMLE estimate for Vg in the null model'):
            vg_lines = True
            vg_header = True
        if line.startswith('## se(Vg)'):
            vg_lines = False
        if line.startswith('## REMLE estimate for Ve in the null model'):
            ve_lines = True
            ve_header = True
        if line.startswith('## se(Ve)'):
            ve_lines = False
        if vg_lines and not vg_header:
            vg_list.append(line.split("\t"))
        if ve_lines and not ve_header:
            ve_list.append(line.split("\t"))
        if vg_header:
            vg_header = False
        if ve_header:
            ve_header = False
        if line.startswith('## total computation time'):
            time = float(line.split(" = ")[1].split(" ")[0])

vg = np.array(pd.DataFrame(vg_list).astype('float32'))
ve = np.array(pd.DataFrame(ve_list).astype('float32'))

vg[np.isnan(vg)] = 0
ve[np.isnan(ve)] = 0

vg_T = vg.copy().T
np.fill_diagonal(vg_T, 0)

ve_T = ve.copy().T
np.fill_diagonal(ve_T, 0)

vg = vg + vg_T
ve = ve + ve_T

np.savetxt("{}/Cg_gemma.csv".format(directory), vg, delimiter=",")
np.savetxt("{}/Cn_gemma.csv".format(directory), ve, delimiter=",")
np.savetxt("{}/VD_time_gemma.csv".format(directory), [time*60], delimiter=",")

