import sys
import pdb
import h5py
import scipy as sp
import pandas as pd
import numpy as np
import os
import time
import gzip
import re
from optparse import OptionParser
# bsub -R 'select[gpfs]' python geno_simulator.py 15000 10 0


class GenoSimulator:

    def __init__(self, fn, nAncestors, nBlocks, maf=None, calc_covar=False,
                 popstructure=False, csv=False):
        """
        args:

        fn         : filename of the true genotypes
        nAncestors : number of parents for each new individual
        nBlock     : number of Snps that are jointly inherited by one parent
        maf        : minor allele frequency filter
        """
        self.fn = fn
        self.nAncestors = nAncestors
        self.nBlocks = nBlocks
        self.calc_covar = calc_covar
        self.popstructure = popstructure
        self.maf = maf
        self.csv = csv

    def generate_genos(self, nSamples, chroms, of, debug=False):
        """
        args:

        nSamples : number of genotypes to be generated
        of       : filename of the output file
        chroms   : chromosomes to use
        """

        f = h5py.File(self.fn, 'r')
        pop_ref = f['genotypes']['chrom1']['row_headers']['population'][:]
        sampleID_ref = f['genotypes']['chrom1']['row_headers']['sampleID'][:]
        f.close()

        if not self.popstructure:
            """ selecting ancestors """
            par_new = sp.zeros((nSamples, self.nAncestors), dtype='S7')
            sampleID_new = sp.array(['ID_%s' % (x + 1)
                                     for x in range(nSamples)], dtype='S7')
            for i in range(nSamples):
                par_new[i] = sp.random.choice(sampleID_ref, self.nAncestors,
                                              replace=False)
        else:
            """ selecting populations """
            pops = sp.unique(pop_ref)
            n_populations = pops.shape[0]
            Icv = sp.floor(n_populations * sp.arange(nSamples) / nSamples)
            pop_new = sp.zeros(nSamples, dtype='S7')
            for pop_i in range(n_populations):
                pop_new[Icv == pop_i] = pops[pop_i]
            """ selecting ancestors within populations """
            par_new = sp.zeros((nSamples, self.nAncestors), dtype='S7')
            sampleID_new = sp.array(['ID_%s' % (x + 1)
                                     for x in range(nSamples)])
            for i in range(nSamples):
                par_new[i] = sp.random.choice(
                        sampleID_ref[pop_ref == pop_new[i]], self.nAncestors, 
                        replace=False)
            idx_pop = sp.argsort(pop_new)
            pop_new = pop_new[idx_pop]
            par_new = par_new[idx_pop]

        if self.calc_covar:
            K = sp.zeros((nSamples, nSamples))

        for chrom_id in chroms:

            # reading in SNP information
            print('.. simulating chrom %d' % chrom_id)
            f = h5py.File(self.fn, 'r')
            pos = f['genotypes']['chrom%d' % chrom_id]['col_headers']['pos'][:]
            chrom = f['genotypes']['chrom%d' %
                                   chrom_id]['col_headers']['chrom'][:]
            rs = f['genotypes']['chrom%d' % chrom_id]['col_headers']['rs'][:]
            f.close()
            nSnps = len(pos)

            """ setting up file per chromosome """
            of_chr = "%s_chr%s.h5" % (of, chrom_id)
            f = h5py.File(of_chr, 'w')
            f_chrom = f.create_group('genotype')
            f_mat = f_chrom.create_dataset('matrix', (nSamples, nSnps),
                                           chunks=(nSamples, self.nBlocks))

            # simulate blocks
            idx_start = 0
            idx_stop = self.nBlocks
            while idx_start < nSnps:
                print(idx_start / float(nSnps))
                X = self.simulateBlock(idx_start, idx_stop, chrom_id, par_new)
                f_mat[:, idx_start:idx_stop] = X.T
                if self.calc_covar:
                    std = X.std(1)
                    I = std > 0
                    std = std[I]
                    Z = X[I, :]
                    Z -= Z.mean(1)[:, sp.newaxis]
                    Z /= std[:, sp.newaxis]
                    K += sp.dot(Z.T, Z)
                idx_start += self.nBlocks
                idx_stop += self.nBlocks
                if idx_stop > nSnps:
                    idx_stop = nSnps
            f_row = f_chrom.create_group('row_header')
            if self.popstructure:
                f_row['population'] = pop_new
            f_row['parents'] = par_new
            f_row['sample_ID'] = sampleID_new

            f_col = f_chrom.create_group('col_header')
            f_col['chrom'] = chrom
            f_col['pos'] = pos
            f_col['rs'] = rs

            f.close()

            if self.maf is not None:
                print("...filtering for maf %s" % self.maf)
                self.mafFilter(of=of, chrom_id=chrom_id)

            if self.csv:
                print("...creating .csv files")
                self.makeText(of=of, chrom_id=chrom_id)

        if self.calc_covar:
            pd.DataFrame(
                K,
                columns=sampleID_new).to_csv(
                "%s_kinship.csv" %
                of,
                sep=",",
                index=False,
                header=True)
            K = K / K.diagonal().mean()
            pd.DataFrame(
                K,
                columns=sampleID_new).to_csv(
                "%s_kinship_norm.csv" %
                of,
                sep=",",
                index=False,
                header=True)

    def getChroms(self, chromstring):
        """
        String of chromosomes to simulate

        Arguments:
            chromstring (string):
                comma-separated chr numbers (for single chr) or hyphen-
                separated chr numbers (for chr ranges) or combination of
                both for chr selection

        Returns:
            (numpy array)
                array containing list of chr IDs
        """
        print('Simulate chromosomes {}'.format(chrstring))
        search = re.compile('[^0-9,-]').search
        if bool(search(chrstring)):
            raise Exception(('chrstring can only contain integers '
                '(0-9), comma (,) and hyphen (-), but {} provided').format(
                    chrstring))
                chrlist = [x.split('-') for x in chrstring.split(',')]
                return np.array(chrlist)

    def simulateBlock(self, idx_start, idx_stop, chrom_id, parents):
        """ Read real genotype data """
        f = h5py.File(self.fn, 'r')
        G_ref = f['genotypes']['chrom%d' %
                               chrom_id]['matrix'][idx_start:idx_stop]
        sampleID_ref = f['genotypes']['chrom%d' %
                                      chrom_id]['row_headers']['sampleID'][:]
        f.close()

        nSamples = parents.shape[0]
        nSnps = G_ref.shape[0]
        G_new = sp.zeros((nSnps, nSamples))

        for i_sample in range(nSamples):
            idx_parent = sampleID_ref == sp.random.choice(parents[i_sample])
            G_new[:, i_sample] = G_ref[:, idx_parent][:, 0]

        return G_new

    def mafFilter(self, of, chrom_id, genotypes=None):
        """ setting up maf-filtered file per chromosome """

        f_chrom = h5py.File("%s_chr%s.h5" % (of, chrom_id), 'r')
        chrom_group_in = f_chrom['genotype']

        f_chrom_maf_name = "%s_chr%s_maf%s.h5" % (of, chrom_id, self.maf)
        f_chrom_maf = h5py.File(f_chrom_maf_name, 'w')
        chrom_group_out = f_chrom_maf.create_group('genotype')

        if genotypes is not None:
            X = genotypes[:]
        else:
            X = chrom_group_in['matrix'][:]

        MAF = (2 * sp.sum(X == 2, axis=0) + sp.sum(X == 1, axis=0)) /\
            float(2. * X.shape[0])
        invert = MAF > 0.5
        MAF[invert] = 1 - MAF[invert]
        filter = MAF > self.maf
        print('%d variants discarded' % (~filter).sum())
        X = X[:, filter]

        row_header = bytes("row_header", 'utf-8')
        chrom_group_out.create_dataset('matrix', data=X)
        h5py.h5o.copy(chrom_group_in.id, row_header,
                      chrom_group_out.id, row_header)
        chrom_group_out.move("/genotype/row_header/sample_ID",
                             "/genotype/row_header/sample_ID")
        colHin = chrom_group_in['col_header']
        colHout = chrom_group_out.create_group('col_header')
        for key in list(colHin.keys()):
            colHout.create_dataset(key, data=colHin[key][:][filter])
        f_chrom_maf.close()

    def makeText(self, of, chrom_id, pos=None, chrom=None, rs=None,
                 genotypes=None):

        f_chrom = h5py.File(
            "%s_chr%s_maf%s.h5" %
            (of, chrom_id, self.maf), 'r')
        chrom_group_in = f_chrom['genotype']

        if genotypes is not None:
            X = genotypes[:]
        else:
            X = chrom_group_in['matrix'][:].astype(int)
            pos = f_chrom['genotype']['col_header']['pos'][:].astype(int)
            chrom = f_chrom['genotype']['col_header']['chrom'][:].astype(int)
            rs = f_chrom['genotype']['col_header']['rs'][:]

        id = []
        for x in range(len(pos)):
            id.append("%s-%s-%s" % (chrom[x], pos[x], rs[x]))

        pd.DataFrame(np.vstack((np.array(id), X)).T).to_csv(
            "%s_chr%s_maf%s.csv" % (of, chrom_id, self.maf), index=False,
            header=False)

    def format(self, of, infile, chrom_id, row_header_in='row_headers',
               row_header_out='row_header', col_header_in='col_headers',
               col_header_out='col_header', sampleID="sampleID",
               genotypes=None):

        f_chrom_in = h5py.File("%s_chr%s.h5" % (infile, chrom_id), 'r')
        chrom_group_in = f_chrom_in['genotype']

        """ setting up maf-filtered file per chromosome """
        f_chrom_out = h5py.File("%s_chr%s.h5" % (of, chrom_id), 'w')
        chrom_group_out = f_chrom_out.create_group('genotype')

        X = chrom_group_in['matrix'][:].T.astype(float)

        chrom_group_out.create_dataset('matrix', data=X)
        h5py.h5o.copy(chrom_group_in.id, row_header_in,
                      chrom_group_out.id, row_header_out)
        chrom_group_out.move("/genotype/row_header/%s" % sampleID,
                             "/genotype/row_header/sample_ID")
        h5py.h5o.copy(chrom_group_in.id, col_header_in,
                      chrom_group_out.id, col_header_out)
        f_chrom_out.close()

    def estimateKinship(self, of, chroms, nSamples):
        K = sp.zeros((nSamples, nSamples))

        for chrom_id in chroms:
            f_chrom = h5py.File("%s_chr%s.h5" % (of, chrom_id), 'r')
            chrom_group_in = f_chrom['genotype']
            X = chrom_group_in['matrix'][:].T
            sample_ID = chrom_group_in['row_header']['sample_ID'][:]

            std = X.std(1)
            I = std > 0
            std = std[I]
            Z = X[I, :]
            Z -= Z.mean(1)[:, sp.newaxis]
            Z /= std[:, sp.newaxis]
            K += sp.dot(Z.T, Z)

        pd.DataFrame(
            K,
            columns=sample_ID).to_csv(
            "%s_kinship.csv" %
            of,
            sep=",",
            index=False,
            header=True)
        K = K / K.diagonal().mean()
        pd.DataFrame(
            K,
            columns=sample_ID).to_csv(
            "%s_kinship_norm.csv" %
            of,
            sep=",",
            index=False,
            header=True)


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--fn1000Genomes", dest='fn1000Genomes', type=str,
                      default=("/gpfs/nobackup/stegle/users/casale/mksum/",
                               "1000G/data/1000G_only4Eu.hdf5"))
    parser.add_option("--nSamples", dest='nSamples', type=int,
                      default=1000)
    parser.add_option("--seed", dest='seed', type=int, default=0)
    parser.add_option("--nAncestors", dest='nAncestors', type=int,
                      default=10)
    parser.add_option("--nBlocks", dest='nBlocks', type=int, default=1000)
    parser.add_option("--outdir", dest='outdir', type=str, default=None)
    parser.add_option("--fname", dest='fname', type=str, default=None)
    parser.add_option("--type", action="store", dest='type', default=None)
    parser.add_option("--chr", action="store", dest='chr', type=str, 
                      default="1-22")
    parser.add_option("--calc_cov", action="store_true", dest='calc_cov',
                      default=False)
    parser.add_option("--maf", action="store", dest='maf', type=float,
                      default=None)
    parser.add_option("--mafonly", action="store_true", dest='mafonly',
                      default=False)
    parser.add_option("--formatonly", action="store_true",
                      dest='formatonly', default=False)
    parser.add_option("--csv", action="store_true", dest='csv',
                      default=False)
    parser.add_option("--csvonly", action="store_true", dest='csvonly',
                      default=False)
    parser.add_option("--kinshiponly", action="store_true",
                      dest='kinshiponly', default=False)
    parser.add_option("--in_dir", dest='in_dir', type=str, default=None)
    parser.add_option("--rowHeaderGenotypes", dest='row_header_in',
                      type=str, default='row_headers')
    parser.add_option("--colHeaderGenotypes", dest='col_header_in',
                      type=str, default='col_headers')
    (opt, args) = parser.parse_args()

    fn1000Genomes = opt.fn1000Genomes
    nSamples = opt.nSamples
    nAncestors = opt.nAncestors
    nBlocks = opt.nBlocks
    seed = opt.seed
    type = opt.type
    maf = opt.maf
    chrstring = opt.chr
    calc_covar = opt.calc_cov
    csv = opt.csv
    if opt.fname is None:
        fname = '%dG_nAnc%d_nBlocks%d_seed%d_%s' % (
            nSamples, nAncestors, nBlocks, seed, opt.type)
    else:
        fname = opt.fname

    outdir = opt.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, fname)

    kinshiponly = opt.kinshiponly
    formatonly = opt.formatonly
    mafonly = opt.mafonly
    csvonly = opt.csvonly


    print('Number of Samples: %d' % nSamples)
    print('Number of Ancestors: %d' % nAncestors)
    print('... setting seed to %d' % seed)

    if type == "pop":
        popstructure = True
    else:
        popstructure = False

    sp.random.seed(seed)
    t1 = time.time()
    sim = GenoSimulator(fn1000G, nAncestors=nAncestors, nBlocks=nBlocks,
                        maf=maf, calc_covar=calc_covar, csv=csv,
                        popstructure=popstructure)
    chroms = sim.getChroms(chrstring)

    if not any([mafonly, formatonly, csvonly, kinshiponly]):
        sim.generate_genos(nSamples=nSamples, chroms=chroms,
                           of=outfile)
    elif mafonly:
        sim.mafFilter(of=out_file, chrom_id=chroms)
    elif csvonly:
        sim.makeText(of=out_file, chrom_id=chroms)
    elif kinshiponly:
        sim.estimateKinship(of=out_file, chroms=chroms,
                            nSamples=nSamples)
    else:
        in_file = os.path.join(opt.in_dir, opt.fname)
        sim.format(of=out_file, chrom_id=chroms, infile=in_file)
    t2 = time.time()
    print(('... finished in %.2f seconds' % (t2 - t1)))
