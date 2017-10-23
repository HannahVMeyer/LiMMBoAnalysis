###################################################
###                                             ###
###     Calibration analysis for LiMMBo         ###
###                                             ###
###     # 1. variance decomposition             ### 
###         with runLiMMBo and                  ###
###         runSimpleVD on:                     ###
###         * different levels of h2            ###
###         * different numbers of              ###
###           traits                            ###
###     # 2. Association study in:              ###
###     setupLiMMBo/calibration_association.sh  ###
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

declare -A kinships=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_kinship_norm.csv"\
    ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv"\
    ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv")

declare -A PCs=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_pcs.txt"\
    ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_pcs.txt"\
    ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_pcs.txt")

declare -A genotypes=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_"\
    ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_"\
    ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_")

#######################
### Job submission  ###
#######################

runningjobs=`bjobs -w | sed 's/  */ /g' | cut -d " " -f 7`
### 1. estimate cg/cn of simulated data for differing number of traits
for kinship in "${!kinships[@]}"; do
    kinshipfile=${kinships[$kinship]}
	for h2 in 0.2 0.5 0.8; do    
		for P in 10 30 50 100; do
			if [[ $P == 10 ]]; then
                p=5
			else 
				p=10
			fi

            if [[ $P -lt 30 ]]; then
                mem=10000
            elif [[ $P -lt 70 ]]; then
                mem=20000
            else
                mem=40000
            fi
			setup=samples1000_traits${P}_Cg${h2}_model$model
			pheno=$simulatePheno/$k/$setup/seed$seed/Ysim_seed${seed}.csv
		
			# create trait-size specific subdirectory
            logname=Traits${P}_sampling${p}_${k}_${h2}_cal
			analysisdir=$dir/$k/$setup/seed$seed/estimateVD/nrtraits_samples$p
			mkdir -p $analysisdir/log

            if [[ ! -e $analysisdir/Cg_all_bootstraps.p && \
 					! ${logname}_limmbo =~ $runningjobs ]]; then
				bsub -q $q -g /limmbo -sp 100 -J ${logname}_limmbo \ 
            	-o $analysisdir/log/${logname}_limmbo.log \ 
            	-e $analysisdir/log/${logname}_limmbo.err \
				-n 4 \
				-R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" \
            	-M $mem \
            	"runLiMMBo -of $analysisdir -lp $P -sp $p \
				-seed $seedLiMMBo -pf $pheno -kf $kinshipfile -cache -v -t"
			fi

            if [[ ! -e $analysisdir/Cn_mtSet.csv && $P < 50 ]]; then
                bsub -q $q  -sp 100 \
                -J Traits${P}_${k}_${h2}_mtSet \
                -o $analysisdir/log/Traits${P}_${k}_${h2}_mtSet.log \
                -e $analysisdir/log/Traits${P}_${k}_${h2}_mtSet.err \
                -R "select[mem>$mem] rusage[mem=$mem]" \
                -M $mem \
                "runSimpleVD -of $analysisdir -pf $pheno -kf $kinshipfile -cache \
                -v -t"
            fi
		done
	done
done
