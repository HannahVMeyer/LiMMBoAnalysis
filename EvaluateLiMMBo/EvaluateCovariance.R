###############################
### Libraries and Functions ###
###############################

library("data.table")
library("ggplot2")
library("wesanderson")
library("plyr")
library("dplyr")

library("reshape2") # melt, acast


SimilarityCovariance <- function(trueDirectory, fitDirectory, mtSetDirectory,
                                 seedLiMMBo, setupString) {
    trueCg_file <- paste(trueDirectory, "/cov_Y_genBg_", setupString,
                         ".csv", sep="")
    trueCn_file <- paste(trueDirectory, "/cov_Y_noiseBg_", setupString,
                                                      ".csv", sep="")
    fitCg_file <- paste(fitDirectory, "/Cg_fit_seed", seedLiMMBo, ".csv",
                        sep="")
    fitCn_file <- paste(fitDirectory, "/Cn_fit_seed", seedLiMMBo, ".csv", 
                        sep="")
    mtSetCg_file <- paste(mtSetDirectory, "/Cg_mtSet.csv", sep="")
    mtSetCn_file <- paste(mtSetDirectory, "/Cn_mtSet.csv", sep="")

    trueCg <- fread(trueCg_file, data.table=FALSE, stringsAsFactors=FALSE)
    trueCn <- fread(trueCn_file, data.table=FALSE, stringsAsFactors=FALSE)
    
    fitCg <- fread(fitCg_file, data.table=FALSE, stringsAsFactors=FALSE)
    fitCn <- fread(fitCn_file, data.table=FALSE, stringsAsFactors=FALSE)
    rssCg_fit <- sum((trueCg - fitCg)^2)
    rssCn_fit <- sum((trueCn - fitCn)^2)
    if (file.exists(mtSetCg_file)) {
        mtSetCg <- fread(mtSetCg_file, data.table=FALSE, 
                         stringsAsFactors=FALSE)
        mtSetCn <- fread(mtSetCn_file, data.table=FALSE, 
                         stringsAsFactors=FALSE)
        rssCg_mtSet <- sum((trueCg - mtSetCg)^2)
        rssCn_mtSet <- sum((trueCn - mtSetCn)^2)
    } else {
        rssCg_mtSet <- NA
        rssCn_mtSet <- NA
    }
    return(c(rssCg_fit, rssCn_fit, rssCg_mtSet, rssCn_mtSet))
}



############
### data ### 
############
seedLiMMBo=29348
N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
model <- "noiseBgOnlygeneticBgOnly"
kinship <- c("unrelatedEU_popstructure", "unrelatedEU_nopopstructure", 
             "relatedEU_nopopstructure")
dirroot <- "~/GWAS/data/LiMMBo/Calibration/samples1000_NrSNP20_Cg"
resultdir <- "~/GWAS/data/LiMMBo/Calibration"
simulationdir <- "~/GWAS/data/LiMMBo/simulateData/phenotypes"
Traits <- seq(10,100,10)

alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

################
### analysis ### 
################

covarianceH2 <- lapply(genVariance, function(h2) {
    covarianceKinship <- lapply(kinship, function(K) {
        covarianceTraits <- lapply(Traits, function(P) {
            if (P==10) sampling <- 5
            if (P>10) sampling <- 10

            fitDirectory <- paste(dirroot, h2, "_model", model, "/", K, 
                                  "/nrtraits", P, 
                                  "/estimateVD/nrtraits_samples", sampling, 
                                  "/seed", seedLiMMBo, sep="")
            trueDirectory <- paste(simulationdir, "/samples", N, "_NrSNP", 
                                   NrSNPs, "_Cg", h2, "_model", model, "/", K,
                                   sep="")
            mtSetDirectory <- paste(dirroot, h2, "_model", model, "/", K, 
                                  "/nrtraits", P, "/closedForm", sep="") 
            setupString <- paste("samples", N, "_traits", P, "_NrSNP", NrSNPs,
                                 "_Cg", h2, "_model", model, sep="")
            rss <- SimilarityCovariance(trueDirectory, fitDirectory, 
                                        mtSetDirectory, seedLiMMBo, 
                                        setupString)
            
            rss_summary <- data.frame(rssCg_fit =rss[1], rssCn_fit=rss[2], 
                                      rssCg_mtSet =rss[3], rssCn_mtSet=rss[4],
                                      P=P, h2=h2, model=model, kinship=K, 
                                      stringsAsFactors=FALSE)
        })
        covarianceTraits <- do.call(rbind, covarianceTraits)
    })
    covarianceKinship <- do.call(rbind, covarianceKinship)
})

covarianceSummary <- do.call(rbind, covarianceH2)
covarianceSummary <- melt(covarianceSummary, 
                          id.vars=c("kinship", "h2", "P", "model"),
                          measured.vars=c("rssCg_fit", "rssCn_fit", 
                                          "rssCg_mtSet", "rssCn_mtSet"),
                          value.name="RSS",
                          variable.name="type")
saveRDS(covarianceSummary, paste(resultdir, "/covarianceSummary.rds", sep=""))

axistext <- 18
axistitle <- 18
legendtext <- 18
legendtitle <- 18
striptext <- 20


png(file=paste(resultdir, "/covarianceSummary.png", sep=""), 
    height=1000, width=1000)
p <- ggplot(covarianceSummary, 
            aes(x=as.factor(P), y=RSS))
p + geom_boxplot(position=position_dodge(width=1)) +
    facet_grid(as.factor(h2) ~ as.factor(kinship)) +
    labs(x = "Number of traits", y = expression(-log[10](FDR))) +
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
          legend.key = element_rect(colour = NA)) 
dev.off()

