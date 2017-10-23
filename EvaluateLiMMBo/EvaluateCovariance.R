###################################################################
###                                                             ###
### Evaluate covariance estimates of REML and LiMMBo compared   ###
### to true covariance matrices							        ###
###                                                             ###
### Data generated via setupLiMMBo/runtime.sh                   ###
###                                                             ###
### Generates Figure 3A, (publication)                          ###
###           Figure 4.3 (thesis)                               ###
###                                                             ###
###################################################################

###############################
### Libraries and Functions ###
###############################

library("data.table")
library("ggplot2")
library("wesanderson")
library("plyr")
library("dplyr")

library("reshape2") # melt, acast


SimilarityCovariance <- function(trueDirectory, fitDirectory,
                                 seedLiMMBo, seedSim, NAvalues=NA) {
    trueCg_file <- paste(trueDirectory, "/cov_Y_genBg_seed", seedSim,".csv", 
                         sep="")
    trueCn_file <- paste(trueDirectory, "/cov_Y_noiseBg_seed", seedSim, ".csv", 
                         sep="")
    fitCg_file <- paste(fitDirectory, "/Cg_fit_seed", seedLiMMBo, ".csv",
                        sep="")
    fitCn_file <- paste(fitDirectory, "/Cn_fit_seed", seedLiMMBo, ".csv", 
                        sep="")
    mtSetCg_file <- paste(fitDirectory, "/Cg_mtSet.csv", sep="")
    mtSetCn_file <- paste(fitDirectory, "/Cn_mtSet.csv", sep="")
    
    trueCg <- fread(trueCg_file, data.table=FALSE, stringsAsFactors=FALSE)
    trueCn <- fread(trueCn_file, data.table=FALSE, stringsAsFactors=FALSE)
    uniqueElements <- length(trueCg[lower.tri(as.matrix(trueCg), diag=TRUE)])

    if (file.exists(fitCg_file)) {
        
        fitCg <- fread(fitCg_file, data.table=FALSE, stringsAsFactors=FALSE)
        fitCn <- fread(fitCn_file, data.table=FALSE, stringsAsFactors=FALSE)
        rssCg_fit <- sum((trueCg[lower.tri(as.matrix(trueCg), diag=TRUE)] - 
                          fitCg[lower.tri(as.matrix(fitCg), diag=TRUE)])^2)
        rssCn_fit <-  sum((trueCn[lower.tri(as.matrix(trueCn), diag=TRUE)] -
                          fitCn[lower.tri(as.matrix(fitCn), diag=TRUE)])^2)
        rmseCg_fit <- sqrt(rssCg_fit/uniqueElements)
        rmseCn_fit <- sqrt(rssCn_fit/uniqueElements)
    } else {
        rssCg_fit <- NA
        rssCn_fit <- NA
        rmseCg_fit <- NA
        rmseCn_fit <- NA
    }
    if (file.exists(mtSetCg_file)) {
        mtSetCg <- fread(mtSetCg_file, data.table=FALSE, 
                         stringsAsFactors=FALSE)
        mtSetCn <- fread(mtSetCn_file, data.table=FALSE, 
                         stringsAsFactors=FALSE)
    	rssCg_mtSet <- sum((trueCg[lower.tri(as.matrix(trueCg), diag=TRUE)] - 
					  mtSetCg[lower.tri(as.matrix(mtSetCg), diag=TRUE)])^2)
    	rssCn_mtSet <-  sum((trueCn[lower.tri(as.matrix(trueCn), diag=TRUE)] -
                      mtSetCn[lower.tri(as.matrix(mtSetCn), diag=TRUE)])^2)
    	rmseCg_mtSet <- sqrt(rssCg_mtSet/uniqueElements)
    	rmseCn_mtSet <- sqrt(rssCn_mtSet/uniqueElements)
                      
    } else {
        rssCg_mtSet <- NAvalues
        rssCn_mtSet <- NAvalues
        rmseCg_mtSet <- NAvalues
        rmseCn_mtSet <- NAvalues
    }
    return(c(rssCg_fit, rssCn_fit, rssCg_mtSet, rssCn_mtSet, 
             rmseCg_fit, rmseCn_fit, rmseCg_mtSet, rmseCn_mtSet))
}

############
### data ### 
############

seedLiMMBo=29348
N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
model <- "noiseBgOnlygeneticBgOnly"
dirroot <- "~/GWAS/data/LiMMBo/Calibration"
simulationdir <- "~/GWAS/data/simulateData/Runtime/phenotypes"
Traits <- seq(10,100,10)
seed <- 1:10
alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)
K="relatedEU_nopopstructure"

# plotting parameters
axistext <- 12
axistitle <- 12
legendtext <- 12
legendtitle <- 12
color <- wes_palette(5, name="Darjeeling", type='continuous')[2:3]

################
### analysis ### 
################

# Compare covariance estimates from REML and LiMMBo to true (known from 
# simulation) covariance matrices

covarianceH2 <- lapply(genVariance, function(h2) {
    covarianceSeed <- lapply(seed, function(S) {
        covarianceTraits <- lapply(Traits, function(P) {
            if (P==10) sampling <- 5
            if (P>10) sampling <- 10
            
            analysisString <- paste(K, "/samples", N, "_traits", P, "_Cg", h2, 
                                    "_model", model, "/seed", S, sep="")
            setupString <- paste("samples", N, "_traits", P, "_NrSNP", NrSNPs,
                                 "_Cg", h2, "_model", model, sep="")
            
            fitDirectory <- paste(dirroot, "/", analysisString,
                                  "/estimateVD/nrtraits_samples", sampling, 
                                  sep="")
            trueDirectory <- paste(simulationdir,"/", analysisString, sep="")

            ss <- SimilarityCovariance(trueDirectory, fitDirectory, seedLiMMBo, 
                                        S, NAvalues=-1)
            
            summary <- data.frame(rssCg_LiMMBo =ss[1], rssCn_LiMMBo=ss[2],
                                  rssCg_RML = ss[3], rssCn_RML=ss[4], 
                                  rmseCg_LiMMBo =ss[5], rmseCn_LiMMBo=ss[6],
                                  rmseCg_RML=ss[7], rmseCn_RML=ss[8],
                                  P=P, h2=h2, model=model, kinship=K, seed=S,
                                  stringsAsFactors=FALSE)
        })
        covarianceTraits <- do.call(rbind, covarianceTraits)
    })
    covarianceSeed <- do.call(rbind, covarianceSeed)
})

covarianceSummary <- do.call(rbind, covarianceH2)
covarianceSummary <- melt(covarianceSummary, 
                          id.vars=c("kinship", "h2", "P", "model", "seed"),
                          measured.vars=c("rssCg_LiMMbo", "rssCn_LiMMBo", 
                                          "rssCg_RML", "rssCn_RML",
                                          "rmseCg_LiMMbo", "rmseCn_LiMMBo",
                                          "rmseCg_RML", "rmseCn_RML"),
                          value.name="SS",
                          variable.name="type")
covarianceSummary$component <- gsub("[rmse]{3,4}(C.*)_.*", "\\1", 
                                    covarianceSummary$type)
covarianceSummary$measure <- gsub("([rmse]{3,4})C.*_.*", "\\1", 
                                    covarianceSummary$type)
covarianceSummary$analysis <- gsub(".*_", "", covarianceSummary$type)
covarianceSummary$type <- gsub("[rmse]{3,4}", "", covarianceSummary$type)
covarianceSummary <- covarianceSummary[!is.na(covarianceSummary$SS),]
saveRDS(covarianceSummary, paste(dirroot, "/covarianceSummary.rds", sep=""))

## Plot root mean suqared deviation of covariance matrices 
## Figure 3A in publication, Figure 4.3 in thesis
p <- ggplot(filter(covarianceSummary, P %in% c(10, 20, 30, 50, 100),
                   measure == "rmse"), 
            aes(x=as.factor(P), y=SS, fill = as.factor(analysis)))
p  <- p + geom_boxplot() +
    scale_fill_manual(name="Method", values=color) + 

    labs(x = "Number of traits", y = "Root Mean Squared Deviation") +
    coord_cartesian(ylim = c(0.35, 0.75)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black",size=axistext,angle=0,
                                     hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=axistext,angle=0,
                                     hjust=0.5,vjust=0.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=axistitle,angle=0,
                                      hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=axistitle,angle=90,
                                      hjust=.5,vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=legendtext,angle=0,
                                     hjust=1,vjust=0,face="plain"),
          legend.title= element_text(colour="black",size=legendtitle,angle=0,
                                     hjust=1,vjust=0,face="plain"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom") 

ggsave(plot=p, file=paste(dirroot, "/covarianceSummary.pdf", sep=""), height=4, 
       width=5.2, units="in")
ggsave(plot=p, file=paste(dirroot, "/covarianceSummary.eps", sep=""), height=4, 
       width=5.2, units="in")
