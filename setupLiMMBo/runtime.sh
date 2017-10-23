###########################################
###                                     ### 
###     Run time analysis for LiMMBo    ### 
###                                     ###
###     # runLiMMBo with:               ###
###         * different levels of h2    ###
###         * different numbers of      ###
###           traits                    ###
###     # repeat analysis for 10        ###
###       different seeds               ###
###                                     ###
###                                     ### 
###########################################

###################
### parameters  ###
###################

# submission parameters
env='/homes/hannah/software/anaconda2/envs/limix/bin'
q=research-rh7

# directories and analysis parameters
dir=~/data/LiMMBo/Runtime
simulateGeno=~/data/simulateData/genotypes
simulatePheno=~/data/simulateData/Calibration/phenotypes
simulateKinship=~/data/simulateData/genotypes/relatedEU_nopopstructure

seedLiMMBo=29348
model=noiseBgOnlygeneticBgOnly
NrSNP=20
k=relatedEU_nopopstructure
kinship=$simulateKinship/1000G_nAnc2_nBlocks1000_seed256_only4Eu_noPop_kinship_norm.csv

#######################
### Job submission  ###
#######################

for seed in `seq 1 10`; do
	for h2 in 0.2 0.5 0.8; do    
		for P in 10 20 30 40 50 60 70 80 90 100; do
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
			analysisdir=$dir/$k/$setup/seed$seed/estimateVD/nrtraits_samples$p
			mkdir -p $analysisdir/log

            runningjobs=`bjobs -w | sed 's/  */ /g' | cut -d " " -f 7`

            logname=Traits${P}_sampling${p}_seed${seed}_${k}_${h2}
            if [[ ! -e $analysisdir/Cg_all_bootstraps.p && \
 					! $runningjobs =~ ${logname}_limmbo ]]; then
				bsub -q $q -g /permutation -sp 100 -J ${logname}_limmbo \ 
            	-o $analysisdir/log/$logname.log \ 
            	-e $analysisdir/log/$logname.err \
            	-R "select[mem>$mem] rusage[mem=$mem]" \ 
            	-M $mem \
            	"runLiMMBo -of $analysisdir -lp $P -sp $p \
				-seed $seedLiMMBo -pf $pheno -kf $kinship -cache -v -t"
			fi

            if [[ ! -e $analysisdir/Cn_mtSet.csv && $P < 50 ]]; then
                bsub -q $q -g /permutation -sp 100 \
                -J Traits${P}_${k}_${h2}_mtSet \
                -o $analysisdir/log/Traits${P}_${k}_${h2}_mtSet.log \
                -e $analysisdir/log/Traits${P}_${k}_${h2}_mtSet.err \
                -R "select[mem>$mem] rusage[mem=$mem]" \
                -M $mem \
                "runSimpleVD -of $analysisdir -pf $pheno -kf $kinship -cache \
                -v -t"
            fi
		done
	done
done

