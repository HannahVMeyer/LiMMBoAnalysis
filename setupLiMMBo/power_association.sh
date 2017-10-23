###########################################
###                                     ###
###     Power analysis for LiMMBo:      ###
###                                     ###
###     # 1. variance decomposition in  ###
###         * setupLiMMBo/power_vd.sh   ###
###     # 2. Association analyses:      ###
###       	* univariate LMM            ###
###         * multivariate LMM          ###
###                                     ###
###########################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/Power
simulateGeno=~/data/simulateData/Power/genotypes
simulatePheno=~/data/simulateData/Power/phenotypes

kinship=~/data/simulateData/genotypes/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv
seedLiMMBo=29348
model=noiseFixedAndBggeneticFixedAndBg
NrSNP=20
fdr=0.001

#######################
### Job submission  ###
#######################

for seed in `seq 1 $NrSeed`; do
    for h2 in 0.2 0.5 0.8; do    
        for P in 10 50 100; do
            if [[ $P == 10 ]]; then
                p=5
                mem=500
            else
                mem=1000
                p=10
            fi
            for t in 0.01 0.04 0.1 0.2 0.4 0.6 0.8 1; do

                setup=samples${N}_traits${P}_NrSNP${NrSNP}_Cg${h2}_model$model
                analysisdir=$dir/$setup/seed$seed/TraitsAffected$t/estimateVD
                logname=Traits${P}_sampling${p}_seed${seed}_${h2}_${t}
                
				pheno=$simulatePheno/$setup/seed$seed/TraitsAffected$t/Ysim_TraitsAffected$t.csv
                cov=$simulatePheno/$setup/seed$seed/TraitsAffected$t/Covs_TraitsAffected$t.csv
                geno=$simulateGeno/$setup/seed$seed/TraitsAffected$t/SNP_NrSNP_TraitsAffected$t.csv
            	cgf=$analysisdir/seedLiMMBo$seedLiMMBo/Cg_fit_seed$seedLiMMBo.csv
            	cnf=$analysisdir/seedLiMMBo$seedLiMMBo/Cn_fit_seed$seedLiMMBo.csv
            
                if [[ ! -s $pheno ]]; then
                    continue
                fi
                
                if [[ ! -s $cgf ]]; then
                    continue
                fi
            	runningjobs=`bjobs -w | sed 's/  */ /g' | cut -d " " -f 7`
            
                if [[ $P -le 30 ]]; then
                    closedformdir=$dir/$setup/seed$seed/closedForm
                    mkdir -p $closedformdir/log
                    # linear mixed  model
                    if [[ ! $runningjobs =~ \
                        "multitrait_lmm_${P}traits_h2${h2}_closedForm_seed${seed}_$t" && \ 
                    ! -s $closedformdir/lmm_mt_pvalue_causalSNPs_closedForm.csv \ 
                    ]]; then
                        echo "multitrait_lmm_${P}traits_h2${h2}_closedForm_seed${seed}_$t"
                        bsub -sp 100 -q $q -g /permutation \ 
                        -J multitrait_lmm_${P}traits_h2${h2}_closedForm_seed${seed}_$t \
                        -o $closedformdir/log/multitrait_lmm_${P}traits.log \
                        -e  $closedformdir/log/multitrait_lmm_${P}traits.err \
                        -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \ 
                        "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                            -c causalSNPs \
                            -cf $cov -pf $pheno -gf $geno -of $closedformdir \
                            -kf $kinship -m multitrait -set lmm -v -fdr $fdr \
                            -tr gaussian -reg -empiricalP -noPlot" >/dev/null
                    fi
                fi 

                # linear mixed model, single-trait
                if [[ ! $runningjobs =~ "singletrait_lmm_${P}traits_h2${h2}_seed${seed}_$t" \
                    && ! -s $analysisdir/lmm_st_pempirical_causalSNPs$fdr.csv ]]; then
                    echo "singletrait_lmm_${P}traits_h2${h2}_seed${seed}_$t"
                    bsub -sp 100 -q $q -g /permutation \
                    -J singletrait_lmm_${P}traits_h2${h2}_seed${seed}_$t \
                    -o $analysisdir/log/singletrait_lmm_${P}traits.log \
                    -e $analysisdir/log/singletrait_lmm_${P}traits.err \
                    -R "select[mem>$mem] rusage[mem=$mem]" \
                    -M $mem \
                    "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                    -c causalSNPs -cf $cov -pf $pheno -gf $geno -of $analysisdir \
                    -kf $kinship -m singletrait -set lmm -fdr $fdr -tr gaussian \
                    -v -reg -noPlot -empiricalP " >/dev/null
                fi 

                # linear mixed model, multi-trait
                if [[ ! $runningjobs =~ "multitrait_lmm_${P}traits_h2${h2}_seed${seed}_$t" \
                    && ! -s $analysisdir/lmm_mt_pempirical_causalSNPs$fdr.csv \
                    && -s $cgf ]]; then
                    echo "multitrait_lmm_${P}traits_h2${h2}_seed${seed}_$t"
                    bsub -sp 100 -q $q -g /permutation \
                    -J multitrait_lmm_${P}traits_h2${h2}_seed${seed}_$t \ 
                    -o $analysisdir/log/multitrait_lmm_${P}traits.log \
                    -e $analysisdir/log/multitrait_lmm_${P}traits.err \
                    -R "select[mem>$mem] rusage[mem=$mem]" \
                    -M $mem \
                    "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                    -c causalSNPs -cf $cov -pf $pheno -gf $geno -of $analysisdir \
                    -cgf $cgf -cnf $cnf -kf $kinship -m multitrait -set lmm \
                    -fdr $fdr -tr gaussian -v -reg -empiricalP -noPlot"  >/dev/null
                fi
            done
		done
	done
done

