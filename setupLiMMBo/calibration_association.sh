###################################################
###                                             ###
###     Calibration analysis for LiMMBo         ###
###                                             ###
###     # 1. variance decomposition in          ###
###     setupLiMMBo/calibration_vd.sh           ###
###     # 2. GWAS with                          ###
###         * LMM REML                          ###
###         * LMM LiMMBo                        ###
###         * LM with PCs                       ###
###                                             ###
###     # repeat analysis for 10                ###
###       different seeds                       ###
###                                             ###
###################################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/CalibrationTmp
simulateGeno=~/data/simulateData/genotypes

model=noiseBgOnlygeneticBgOnly
NrSNP=20
seed=364
seedLiMMBo=29348
seedpermute=2040
model=noiseBgOnlygeneticBgOnlyi

declare -A
kinships=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_kinship_norm.csv"\
        ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv"\
            ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv")

declare -A
PCs=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_pcs.txt"\
        ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_pcs.txt"\
            ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_pcs.txt")

declare -A
genotypes=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_"\
        ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_"\
            ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_")

#######################
### Job submission  ###
#######################

runningjobs=`bjobs -w | sed 's/  */ /g' | cut -d " " -f 7`
for chr in `seq 1 22`; do
    mem=30000
    if [[ $chr -lt  11 ]]; then
        mem=50012
    fi
	for k in "${!kinships[@]}"; do
		kinship=${kinships[$k]}
		pcs=${PCs[$k]}
		geno=${genotypes[$k]}chr${chr}_maf0.02.h5
		
		for h2 in 0.2 0.5 0.8; do    

			for P in 10 20 30 40 50 60 70 80 90 100; do

				# extract number of required bootstrap runs from 
                # Bootstrap_sampling_scheme generated via BootstrapEstimates
				if [[ $P == 10 ]]; then
					p=5
				else 
					p=10
				fi
				
				setup=samples1000_NrSNP${NrSNP}_Cg${h2}_model$model
				pheno=$simulatePheno/$setup/$k/Ysim_samples1000_traits${P}_NrSNP${NrSNP}_Cg${h2}_model${model}.csv
				analysisdir=$dir/$setup/$k/nrtraits$P/estimateVD/nrtraits_samples$p
				
				if [[ $P -le 30 ]]; then
					closedformdir=$dir/$setup/$k/nrtraits$P/closedForm
					mkdir -p $closedformdir/log
					# linear mixed  model
					if [[ ! $runningjobs =~ "multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2}_closedForm" && ! -s $closedformdir/lmm_mt_pvalue_chr${chr}_${k}_closedForm.csv ]]; then
						echo "multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2}_closedForm"
						bsub -sp 100 -q $q -g /permutation \
                        -J multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2}_closedForm \
                        -o $closedformdir/log/multitrait_chr${chr}_lmm_${P}traits_${k}.log \
                        -e $closedformdir/log/multitrait_chr${chr}_lmm_${P}traits_${k}.err \
                        -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \
                        "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                        -c chr$chr -pf $pheno -gf $geno -of $closedformdir \
                        -kf $kinship -m multitrait -set lmm -v \
                        -fileend _${k}_closedForm" >/dev/null
					fi
				fi 

				# linear  model, with pcs
				if [[ ! $runningjobs =~ "multitrait_samples_chr${chr}_lm_pc_${P}traits_${k}_h2${h2}" && ! -s $analysisdir/lm_mt_pcs_pvalue_chr${chr}_${k}.csv ]]; then
					echo "multitrait_samples_chr${chr}_lm_pc_${P}traits_${k}_h2${h2}"
					bsub -sp 100 -q $q -g /permutation \
                    -J multitrait_samples_chr${chr}_lm_pc_${P}traits_$k_h2${h2} \
                    -o $analysisdir/log/multitrait_chr${chr}_lm_pc_${P}traits_${k}.log \
                    -e $analysisdir/log/multitrait_chr${chr}_lm_pc_${P}traits_${k}.err\
                    -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \
                    "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                    -c chr$chr -pf $pheno -gf $geno -of $analysisdir \
                    -m multitrait -set lm -v -pcsf $pcs -fileend _${k} \
                    -re" >/dev/null
				fi 

				# linear mixed  model
				if [[ ! $runningjobs =~"multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2}" && ! -s $analysisdir/lmm_mt_pvalue_chr${chr}_${k}.csv ]]; then
					echo "multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2}"
					bsub -sp 100 -q $q -g /permutation \
                        -J multitrait_samples_chr${chr}_lmm_${P}traits_${k}_h2${h2} \
                        -o $analysisdir/log/multitrait_chr${chr}_lmm_${P}traits_${k}.log \
                        -e $analysisdir/log/multitrait_chr${chr}_lmm_${P}traits_${k}.err \
                        -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \
                        "$env/python /homes/hannah/GWAS/analysis/GWAS/gwas.py \
                        -c chr$chr -pf $pheno -gf $geno -of $analysisdir \
                        -kf $kinship \
                        -cgf $analysisdir/seed$seed/Cg_fit_seed$seed.csv \
                        -cnf $analysisdir/seed$seed/Cn_fit_seed$seed.csv \
                        -m multitrait -set lmm -v -fileend _${k}" >/dev/null
				fi	
			done
		done
	done
done

