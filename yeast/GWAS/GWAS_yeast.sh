###########################################################
###                                                     ###
### GWAS of 41 yeast quantitative growth traits         ###
###                                                     ###
###     * dataset from Bloom et al 2013                 ###
###     * phenotypes processed via                      ###
###         yeast/phenotypes/phenotypes_yeast.R         ###
###     * genotypes processed via                       ###
###         yeast/genotypes/genotypes_yeast.R           ###
###     * relationship estimated via                    ###
###         yeast/genotypes/relationship_yeast.R        ###
###                                                     ###
###     * GWAS via gwas.py								###
###         * univariate LMM                            ###
###         * multivariate LMMs with LiMMBo             ###
###                                                     ###
###                                                     ###
###########################################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/yeast
file_geno=$dir/inputdata/BYxRM.geno.limix.format.h5
file_pheno=$dir/inputdata/BYxRM.pheno.limix.format.h5
file_kinship=$dir/inputdata/BYxRM.grm.rel.csv
file_Cg=$dir/vd/subset_seed234/Cg_fit.csv
file_Cn=$dir/vd/subset_seed234/Cn_fit.csv
outdir=$dir/GWAS

seed=234
fdr=0.001
cpus=4
mem=10000

#######################
### Job submission  ###
#######################

# LiMMBo
logname=yeast_limmbo_seed$seed
bsub -q $q -g /permutation -sp 100 -J $logname \
-o $dir/vd/seedLiMMBo$seed/log/$logname.log \
-e $dir/vd/seedLiMMBo$seed/log/$logname.err \
-n $cpus \
-R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" \
-M $mem \
"runLiMMBo -of $dir/vd/seedLiMMBo$seed -lp 41 \
-sp 10 -seed $seed -pf $file_pheno -kf $file_kinship \
-tr gaussian -cpus $cpus -cache -v -t"


# multitrait
bsub -sp 100 -q $q -g /permutation -J multitrait_samples_lmm_yeast \
    -o $outdir/log/multitrait_lmm_yeast.log  \
    -e $outdir/log/multitrait_lmm_yeast.err \
	-w "ended($logname)" \
    -R "select[mem>$mem] rusage[mem=$mem]" \
    -M $mem "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
    -c genome  -pf $file_pheno -gf $file_geno -of $outdir \
    -kf $file_kinship -cgf $file_Cg -cnf $file_Cn -m multitrait \
    -set lmm -v -noPlot -fdr $fdr" >/dev/null

# singletrait
bsub -sp 100 -q $q -g /permutation -J singletrait_samples_lmm_yeast \
    -o $outdir/log/singletrait_lmm_yeast.log  \
    -e $outdir/log/singletrait_lmm_yeast.err \
    -R "select[mem>$mem] rusage[mem=$mem]" \
    -M $mem "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
    -c genome  -pf $file_pheno -gf $file_geno -of $outdir -kf $file_kinship \
    -cgf $file_Cg -cnf $file_Cn -m singletrait -set lmm -v \
    -noPlot -fdr $fdr" >/dev/null

