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
env='/homes/hannah/software/anaconda2/envs/limmbo-py2.7/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/yeast
#file_geno=$dir/inputdata/BYxRM.geno.limix.format.h5
file_geno=$dir/inputdata/BYxRM
#file_pheno=$dir/inputdata/BYxRM.pheno.limix.format.h5
pheno_prefix=$dir/inputdata/BYxRM_pheno
file_kinship=$dir/inputdata/BYxRM.3kb.grm.rel.csv
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

for method in mice phenix; do
    mkdir -p $dir/vd/seedLiMMBo$seed/$method/log
    mkdir -p $dir/vd/reml/$method/log
    logname=yeast_limmbo_seed${seed}_$method

    bsub -q $q -g /permutation -sp 100 -J $logname \
    -o $dir/vd/seedLiMMBao$seed/$method/log/$logname.log \
    -e $dir/vd/seedLiMMBol$seed/$method/log/$logname.err \
    -n $cpus \
    -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" \
    -M $mem \
    "runVarianceEstimation -o $dir/vd/seedLiMMBo$seed/$method \
    -sp 10 -seed $seed -p ${pheno_prefix}_$method.csv \
    -k $file_kinship --minCooccurrence 3 \
    -tr gaussian -cpus $cpus --limmbo -v -t"

    logname=yeast_reml_$method
    bsub -q $q -g /permutation -sp 100 -J $logname \
    -o $dir/vd/reml/$method/log/$logname.log \
    -e $dir/vd/reml/$method/log/$logname.err \
    -R "select[mem>$mem] rusage[mem=$mem]" \
    -M $mem \
    "runVarianceEstimation -o $dir/vd/reml/$method \
     -seed $seed -p ${pheno_prefix}_$method.csv \
    -k $file_kinship  \
    -tr gaussian --reml -v -t"
done


for method in mice phenix; do
    outdir=$dir/GWAS/$method
    mkdir -p $outdir/log
    logname=yeast_limmbo_$method
    cg=$dir/vd/seedLiMMBo$seed/$method/Cg_fit_seed$seed.csv
    cn=$dir/vd/seedLiMMBo$seed/$method/Cn_fit_seed$seed.csv

# multitrait
    bsub -sp 100 -q $q -g /permutation -J multitrait_samples_lmm_yeast \
        -o $outdir/log/multitrait_lmm_yeast.log  \
        -e $outdir/log/multitrait_lmm_yeast.err \
        #-w "ended($logname)" \
        -R "select[mem>$mem] rusage[mem=$mem]" \
        -M $mem "$env/python runAssociation -mt -lmm\
        -p ${pheno_prefix}_$method.csv -k $file_kinship \
        -g $file_geno -o $outdir \
        -cg $cg -cn $cn -tr gaussian \
        --plot -v -fdr $fdr" >/dev/null

# singletrait
bsub -sp 100 -q $q -g /permutation -J singletrait_samples_lmm_yeast \
    -o $outdir/log/singletrait_lmm_yeast.log  \
    -e $outdir/log/singletrait_lmm_yeast.err \
    -R "select[mem>$mem] rusage[mem=$mem]" \
    -M $mem "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
    -c genome  -pf $file_pheno -gf $file_geno -of $outdir -kf $file_kinship \
    -cgf $file_Cg -cnf $file_Cn -m singletrait -set lmm -v \
    -noPlot -fdr $fdr" >/dev/null
done


bsub -sp 100 -q $q -g /permutation -J multitrait_samples_lmm_yeast \
        -o $outdir/log/multitrait_lmm_yeast.log  \
        -e $outdir/log/multitrait_lmm_yeast.err \
        #-w "ended($logname)" \
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

