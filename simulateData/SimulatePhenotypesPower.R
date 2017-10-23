###################################################################
###                                                             ###
###     Simulate phenotypes for power analysis                  ###
###                                                             ###
###     Called by simulateData/runSimulatePhenotypesPower.sh    ###
###                                                             ###
###                                                             ###
###################################################################

###############################
### Libraries and Functions ###
###############################

library ("R.methodsS3")
library ("R.oo")
library ("R.utils")
library("PhenotypeSimulator")

simulate <- function(percentTraitsAffected, N, P, genoFilePrefix, 
                     genoFileSuffix, NrCausalSNPs, NrCausalChrom, kinshipfile, 
                     genVar, h2s, directoryGeno, directoryPheno, subdirectory, 
                     seed, phi, delta, rho, NrConfounders, 
                     pIndependentConfounders, pcorr) {
    verbose=TRUE
    Pgenetic <- round(percentTraitsAffected * P)
    if (Pgenetic < 1) {
        return(NULL)
    }
    message("Number of traits with fixed genetic effect: ", Pgenetic, " (", 
            percentTraitsAffected, "%)")
    # set genetic and noise models
    modelGenetic <- "geneticFixedAndBg"
    if (rho == 0) {
        modelNoise <- "noiseFixedAndBg"
    } else {
        modelNoise <- "noiseFixedAndBgAndCorrelated"
    }
    
    # kinship estimate based on standardised SNPs (as described in )
    kinship <- getKinship(kinshipfile = kinshipfile, norm=norm, 
                          verbose = verbose)
    
    # simulate fixed genetic effects (from non-standardised SNP genotypes)
    causalSNPs <- getCausalSNPs(NrCausalSNPs = NrCausalSNPs, 
                                NrChrCausal = NrCausalChrom,
                                genoFilePrefix = genoFilePrefix, 
                                genoFileSuffix = genoFileSuffix,  
                                standardise = FALSE, 
                                genoFileDelimiter = ",", verbose=verbose)
    
    genFixed <- geneticFixedEffects(N=N, P=Pgenetic, X_causal=causalSNPs)  
    if (!is.null(genFixed$shared)) {
        genFixed$shared <- cbind(genFixed$shared, matrix(0, ncol=P-Pgenetic, 
                                                     nrow=N))
    }
    if (!is.null(genFixed$independent)) {
        genFixed$independent <- cbind(genFixed$independent, 
                                  matrix(0, ncol=P-Pgenetic, nrow=N))
    }
    genFixed$cov_effect <- t(cbind(t(genFixed$cov_effect), 
                                   matrix(0, ncol=P-Pgenetic, 
                                          nrow=NrCausalSNPs)))
    
    # simulate random genetic effects
    genBg <- geneticBgEffects(kinship = kinship, P=P)
    
    # simulate fixed noise effects:
    noiseFixed <- noiseFixedEffects(N=N, P=P, NrConfounders=NrConfounders, 
                            pIndependentConfounders=pIndependentConfounders)
    # simulate correlated background if specified:
    if (rho != 0 ) {
        noiseCorrelated <- correlatedBgEffects(N, P, pcorr)
    } else {
        noiseCorrelated <- NULL
    }

    # simulate random noise effects
    noiseBg <- noiseBgEffects(N=N, P=P)
    
    # combine components into final phenotype with genetic variance component 
    # explaining 40% of total variance
    phenotype <- createPheno(N=N, P=P, noiseBg = noiseBg, 
                             noiseFixed = noiseFixed,
                             genFixed = genFixed, genBg = genBg, 
                             correlatedBg = noiseCorrelated,
                             modelNoise = modelNoise, 
                             modelGenetic = modelGenetic, 
                             genVar = genVar, h2s = h2s, phi = phi, 
                             delta = delta, rho = rho, pcorr = pcorr,  
                             theta = 1, verbose = verbose)
    
    # save phenotypes
    savePheno(phenotype, directoryGeno=directoryGeno, 
              directoryPheno=directoryPheno,
              outstring=paste("TraitsAffected", percentTraitsAffected, sep=""),
              saveAsTable=TRUE, saveAsRDS=TRUE, verbose=verbose)
}

############
### data ###
############

### read command line arguments
args <- commandArgs(asValue= TRUE, defaults=list(NrCausalChrom=5, seed=234,
                                                 phi=0.6, delta=0.4,
                                                 rho=0,
                                                 NrConfounders=10,
                                                 pIndependentConfounders=0.5,
                                                 pcorr=0.6))
cat(unlist(args), "\n")

N <- as.numeric(args$NrSamples)
P <- as.numeric(args$NrPhenotypes)
genoFilePrefix <- args$genoFilePrefix 
genoFileSuffix <- args$genoFileSuffix 
NrCausalSNPs <- as.numeric(args$cNrSNP)
NrCausalChrom <- as.numeric(args$NrCausalChrom)
kinshipfile <- args$kinshipfile
genVar <- as.numeric(args$h2)
h2s <- as.numeric(args$h2s)
directoryGeno <- args$directoryGeno
directoryPheno <- args$directoryPheno
seed <- as.numeric(args$seed)
norm <- as.logical(args$norm)
phi <- as.numeric(args$phi) 
delta <- as.numeric(args$delta)
rho <- as.numeric(args$rho)
pcorr <- as.numeric(args$pcorr)
NrConfounders <- as.numeric(args$NrConfounders)
pIndependentConfounders <- as.numeric(args$pIndependentConfounders)

################
### Analysis ###
################

traitModel <- c(0.01, 0.04, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
phenoSimulation <- sapply(traitModel, simulate, N, P, genoFilePrefix, 
       genoFileSuffix, NrCausalSNPs, NrCausalChrom, kinshipfile, 
       genVar, h2s, directoryGeno, directoryPheno, subdirectory, 
       seed, phi, delta, rho, NrConfounders, pIndependentConfounders, pcorr)


