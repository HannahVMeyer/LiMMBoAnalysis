###################################################
###                                     		###
###     Power analysis for LiMMBo:      		###
###                                     		###
###     # 1. variance decomposition in  		###
###         * different numbers of traits		###
###         * different proportions of traits   ###
###			  affected							###
###         * different levels of genetics		###
###     # 2. Association analyses in      		###
###         * setupLiMMBo/power_association.sh  ###
###                                     		###
###################################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/Power
simulateGeno=~/data/simulateData/genotypes
simulatePheno=~/data/simulateData/Power/phenotypes

kinship=$simulateGeno/relatedEU_nopopstructure/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv
model=noiseFixedAndBggeneticFixedAndBg
NrSNP=20
N=1000
NrSeed=50
cpus=4

seedLiMMBo=29348

#######################
### Job submission  ###
#######################

runningjobs=`bjobs -w | sed 's/  */ /g' | cut -d " " -f 7`

### 1. estimate cg/cn of simulated data for differing number of traits
for seed in `seq 1 $NrSeed`; do
    for P in 10 20 30 40 50 60 70 80 80 100; do
        for h2 in 0.2 0.4 0.5 0.6 0.8; do
	# for h2 in 0.2 0.5 0.8; do    
	#	for P in 10 50 100; do
			if [[ $P < 60 ]]; then
                mem=10000
			else
                mem=20000
			fi
            if [[ $P == 10 ]]; then
				p=5
            else
                p=10
			fi
            for t in 0.01 0.04 0.1 0.2 0.4 0.6 0.8 1; do

                setup=samples${N}_traits${P}_NrSNP${NrSNP}_Cg${h2}_model$model
                pheno=$simulatePheno/$setup/seed$seed/TraitsAffected$t/Ysim_TraitsAffected$t.csv
                cov=$simulatePheno/$setup/seed$seed/TraitsAffected$t/Covs_TraitsAffected$t.csv
                
                if [[ ! -s $pheno ]]; then
                    continue
                fi

                # create trait-size specific subdirectory
                analysisdir=$dir/$setup/seed$seed/TraitsAffected$t/estimateVD
                mkdir -p $analysisdir/seedLiMMBo$seedLiMMBo/log
                mkdir -p $analysisdir/log

                logname=Traits${P}_sampling${p}_seed${seed}_${h2}_${t}
                if [[ ! -e $analysisdir/seedLiMMBo$seedLiMMBo/Cg_fit_seed$seedLiMMBo.csv ]]; then
                    if [[ ! $runningjobs =~ "$logname" ]]; then
                        bsub -q $q -g /limmbo -sp 100 -J $logname \ 
                        -o $analysisdir/seedLiMMBo$seedLiMMBo/log/$logname.log \ 
                        -e $analysisdir/seedLiMMBo$seedLiMMBo/log/$logname.err \
                        -n $cpus \
                        -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" \ 
                        -M $mem \
                        "runLiMMBo -of $analysisdir/seedLiMMBo$seedLiMMBo -lp $P \
                        -sp $p -seed $seedLiMMBo -pf $pheno -kf $kinship -cf $cov \
                        -tr gaussian -cpus $cpus -reg -cache -v -t"
                    fi
                fi
            done
		done
	done
done

