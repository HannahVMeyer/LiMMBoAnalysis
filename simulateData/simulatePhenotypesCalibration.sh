
###########################################################
###                                                 	###
###     Simulate phenotypes for calibration analysis	###
###                                                 	###
###     Set up parameter combinations to be used    	###
###     with 'PhenotypeSimulator::simulatePhenotypes()'	###
###                                                 	###
###                                                 	###
###########################################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7
jg="/simulatePhenoLiMMBo"
mem=1000

# directories and analysis parameters
dir=~/data/simulateData
simulateGeno=$dir/genotypes

declare -A kinships=(["unrelatedEU_popstructure"]="$simulateGeno/unrelatedEU_popstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_kinship.csv"\
    ["unrelatedEU_nopopstructure"]="$simulateGeno/unrelatedEU_nopopstructure/1000G_nAnc10_nBlocks1000_seed256_only4Eu_noPop_kinship.csv"\
    ["relatedEU_nopopstructure"]="$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship.csv")

N=1000
model="noiseBgOnlygeneticBgOnly"
seed=364

#######################
### Job submission  ###
#######################

for k in "${!kinships[@]}"; do
    kinshipfile=${kinships[$k]}
    for P in 10 30 50 100; do
        for h2 in 0.2 0.5 0.8; do
            outstring=samples${N}_traits${P}_Cg${h2}_model$model
            analysis=samples${N}_traits${P}_Cg${h2}_model${model}
            subdir=$dir/Calibration/phenotypes/$k/$outstring
            
            # make subdirectory
            mkdir -p $subdir/log
            bsub -g $jg -q research-rh7 -J $analysis \
            -o $subdir/log/${analysis}_simulate.log \  
            -e $subdir/log/${analysis}_simulate.err \  
            -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \ 
            "Rscript -e 'PhenotypeSimulator::simulatePhenotypes()' \
                --args \
                --NrSamples=$N --NrPhenotypes=$P \
                --kinshipfile=$kinshipfile \
                --kinshipfileHasHeader \
                --seed=$seed \
				--cNrSNP=0 \
                --genVar=$h2 --h2bg=1 \
                --phi=1 \
                --directoryGeno=$dir/Calibration/genotypes/$k/$outstring \
                --directoryPheno=$dir/Calibration/phenotypes/$k/$outstring \
				--subdirectory=seed$seed \
				--norm \
                --saveTable \
				--saveRDS \
                --showProgress" 
        done
    done
done
