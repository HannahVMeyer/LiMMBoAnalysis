import h5py
import scipy as sp
import pandas as pd
import numpy as np
import os
import time
import pdb

class ratData:

    def computeMAF(self, X):
        bools = np.logical_and(X >= 0.5, X < 1.5)
        maf = 1 - (2*sp.sum(X>=1.5, axis=1) + sp.sum(bools, axis=1)) / \
                float(2.*X.shape[1])
        return maf

    def filterMAF(self, X, maf, id, mafThr):
        filter = np.logical_and(maf > mafThr, maf < 1 -mafThr)
        print('%d variants discarded' % (~filter).sum())
        X = X[filter,:]
        maf = maf[filter]
        id = id[filter]
        return X, maf, id

    def estimateKinship(self, X, maf):
        std = np.sqrt(2*maf*(1-maf))
        centered = X - 2*maf[:, np.newaxis]
        norm = centered/std[:, np.newaxis]
        K = sp.dot(norm.T,norm)
        return K

    def process(self, infile, outfile, estimateKin=False):
        f_chrom = h5py.File(infile,'r')
        chrom_group_in = f_chrom['imputed_genotypes']
        samples = chrom_group_in["row_header"]["rat"][:]
        nSamples = len(samples)
        snps_all = 0
        snps_maf = 0

        kinship = sp.zeros((nSamples, nSamples))
        
        for chrom in range(21):
            chrom_id = 'chr%s' % (chrom + 1)
            X = chrom_group_in[chrom_id]['array'][:].astype(float)
            pos = chrom_group_in[chrom_id]['col_header']['pos'][:].astype(int)
            id = np.array(["%s-%s-%s:%s" % ((chrom + 1), x, (chrom +1), x) 
                for x in pos])
            snps_all += len(id)
            pd.DataFrame(X, index=id, columns=samples).to_csv("%s_%s.csv" % (
                outfile, chrom_id), index=True, header=True)
            
            maf_all = self.computeMAF(X)
            pd.DataFrame(maf_all, index=id).to_csv(
                    "%s_%s_alleleFrequencies.csv" % (outfile, chrom_id),
                    index=True, header=False)
            X_maf, maf, id = self.filterMAF(X, maf=maf_all, id=id, 
                    mafThr=0.05)
            snps_maf += len(id)
            pd.DataFrame(X_maf, index=id, columns=samples).to_csv(
                    "%s_%s_maf0.05.csv" % (outfile, chrom_id), index=True, 
                    header=True)
            pd.DataFrame(maf, index=id).to_csv(
                    "%s_%s_alleleFrequencies0.05.csv" % (outfile, chrom_id),
                    index=True, header=False)
            
            if (estimateKin):
                kinship += self.estimateKinship(X_maf, maf)

        if (estimateKin):
            pd.DataFrame(kinship, columns=samples).to_csv("%s_kinship.csv" % 
                    outfile, sep=",", index=False, header=True)
            kinship = kinship/float(snps_maf)
            kinship = kinship + 10e-4*sp.eye(nSamples)
            pd.DataFrame(kinship, columns=samples).to_csv(
                    "%s_kinship_norm.csv" % outfile, sep=",", index=False,
                    header=True)
            
if __name__ == "__main__":

    in_file = os.path.join(snakemake.input.hf5file)
    out_file = os.path.join(snakemake.params.dir, snakemake.params.name)
    estimateKin = snakemake.params.estimateKin == "True"

    t1 = time.time()
    data = ratData()
    data.process(outfile=out_file, infile=in_file, estimateKin=estimateKin)
    t2 = time.time()
    print('... finished in %.2f seconds'%(t2-t1))

