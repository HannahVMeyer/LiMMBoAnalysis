This folder contains files relating to this gemma-discussion thread:
https://groups.google.com/forum/#!topic/gemma-discussion/wnt0bY-v2fw

I am using gemma on redhat 7 system, XeonE526 processor.

1. info.txt contains:
  ```bash
  gemma0.96 > info.txt
  gemma0.97 >> info.txt
  ```
1. Test files
  1. Ysim_reg_gemma.txt: space-separated [N x P] matrix of N=1000 samples and P=10 phenotypes
  1. Kinship_eigenvec_gemma.txt: space-separated [N x N] matrix of eigenvectors of kinship matrix
  1. Kinship_eigenval_gemma.txt: space-separated [N x 1] matrix of eigenvalues of kinship matrix
  1. Genotypes_gemma.csv: comma-separated [(S + 3) x N] matrix of S genotypes for N samples, with columns 1-3 containing SNP information

1. test.sh
  Gemma commands tested

1. gemma0.97.err and gemma0.97.log error and log file for test run with gemma
1. gemma0.96.err and gemma0.96.log error and log file for test run with gemma0.96
1. gemma0.96_use_head.err and gemma0.96_use_head.log error and log file for test run with gemma0.96 and genotype size reduction via head
