#######################################################
###                                                 ###  
###     Simulate phenotypes for power analysis      ###
###                                                 ###  
###		Set up parameter combinations to be used    ###
###     for simulateData/SimulatePhenotypesPower.sh ###
###                                                 ###
###                                                 ###
#######################################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7
jg="/simulatePhenoLiMMBo"
mem=3000

# directories and analysis parameters
dir=~/data/simulateData/Power
in=~/data/simulateData/genotypes/relatedEU_nopopstructure

kinshipfile=$in/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv
genoFilePrefix=$in/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_
genoFileSuffix=_maf0.02.csv

N=1000
totalSNPeffect=0.01
NrSNP=20
model="noiseFixedAndBggeneticFixedAndBg"

#######################
### Job submission  ###
#######################

for seed in `seq 1 50`; do for P in 10 20 30 40 50 60 70 80 80 100; do
        for h2 in 0.2 0.4 0.5 0.6 0.8; do
            h2s=`echo "$totalSNPeffect / $h2" | bc -l`
            outstring=samples${N}_traits${P}_NrSNP${NrSNP}_Cg${h2}_model$model
            analysis=${outstring}_seed$seed
            subdir=$dir/phenotypes/$outstring/seed$seed
            mkdir -p $subdir/log
            
            if [[ ! -e $subdir/TraitsAffected1/Ysim_TraitsAffected1.csv ]]; then
                bsub -g $jg -q research-rh7 -J $analysis \
                -o $subdir/log/${analysis}_simulate.log \  
                -e $subdir/log/${analysis}_simulate.err \  
                -R "select[mem>$mem] rusage[mem=$mem]" -M $mem \ 
                "Rscript ~/LiMMBo/simulateData/SimulatePhenotypesPower.R \
                    --args \
                    --NrSamples=$N --NrPhenotypes=$P \
                    --cNrSNP=20 \
                    --genoFilePrefix=$genoFilePrefix \
                    --genoFileSuffix=$genoFileSuffix \
                    --NrCausalChrom=5 \
                    --kinshipfile=$kinshipfile \
                    --seed=$seed \
                    --genVar=$h2 --h2s=$h2s \
                    --directoryGeno=$dir/genotypes/$outstring/seed$seed \
                    --directoryPheno=$dir/phenotypes/$outstring/seed$seed \
                    --norm"
            fi
        done
    done
done
