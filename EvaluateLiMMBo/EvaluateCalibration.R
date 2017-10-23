###################################################################
###                                                             ###
### Evaluate calibration of genetic association studies         ###
### simple LM and multivariate LMMs with REML andLiMMBo         ###
###                                                             ###
### Data generated via setupLiMMBo/calibration_vd.sh            ###
###                and setupLiMMBo/calibration_association.sh   ###
###                                                             ###
### Generates Figure 3B, Table S4 (publication)                 ###
###           Figure 4.4, 4.5 (thesis)                          ###
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

# remove Null list entries 
rmNulls <- function(compList) {
    nonNulls <- compList[!sapply(compList, is.null)]
}

# Assess calibration of p-values per FWER: ratio of SNPs significant for 
# threshold alpha divided by overall number of tests
calibration <- function (pvalues, alpha=c( 5e-5, 5e-6, 5e-7, 5e-8)) {
    snps <- length(pvalues)
    cal <- sapply(alpha, function(a, snps) {
            length(which(pvalues < a))/snps
        }, snps=snps)
    if (any(is.na(cal)))  {
        cal[is.na(cal)] <- NA
    }
    return(cal)
}

# wrapper to read in chromosome-wide GWAS results
readData <- function(gwas, directory, K, estimateString) {
	gwasdata <- lapply(1:22, function(chr, gwas) {
		missing <- NA
		data <- NA
		lm_file <- paste(directory, "/", gwas, "_pvalue_chr", chr,"_", K,
						 estimateString,  ".csv", sep="")

		if(file.exists(lm_file)) {
			data <- fread(lm_file, data.table=FALSE, stringsAsFactors=FALSE)$P
		} else {
			missing <- chr
		}
		return(list(data, missing))
	}, gwas=gwas)
    pvalues <- unlist(sapply(gwasdata, function(p) p[[1]]))
    pvalues <- pvalues[!is.na(pvalues)]
}

# format GWAS results: add expected pvalues and confidence intervals, transform
# into data.frame format suitable for ggplot
formatPValues <- function(cal) {
	analysis <- gsub("([RMLiBo]{2,6})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", "\\1",
					 names(cal))
	traits <- gsub("([RMLiBo]{2,6})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", "\\2",
				   names(cal))
	h2 <- gsub("([RMLiBo]{2,6})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", "\\3",
			   names(cal))
	kinship <- gsub("([RMLiBo]{2,6})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)",
                    "\\4_\\5",
					names(cal))

    # Generate expected p-values distribution and confidence intervals
	calOrdered <- apply(cal, 2, function(x) x[order(x)])
    S <- nrow(calOrdered)
    ci <- 0.95
	expected <- ppoints(S)
    clower   <- -log10(qbeta(ci,     1:S, S - 1:S + 1))
    cupper   <- -log10(qbeta(1 - ci, 1:S, S - 1:S + 1))

	cal <- melt(calOrdered, value.name="observed")[,-1]
	names(cal)[1] <- "type"
	cal$analysis <- rep(analysis, each=nrow(calOrdered))
	cal$traits <- rep(traits, each=nrow(calOrdered))
	cal$h2 <-rep(h2, each=nrow(calOrdered))
	cal$kinship <- rep(kinship, each=nrow(calOrdered))
	cal$expected <- expected
    cal$clower <- clower
    cal$cupper <- cupper
	
	cal$h2 <- factor(cal$h2, labels=paste("h[2]:", unique(cal$h2)))
	cal$traits <- factor(cal$traits, labels=paste("Traits:", 
												   unique(cal$traits)))
	return(cal)
}

############
### data ### 
############

N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
model <- "noiseBgOnlygeneticBgOnly"
kinship <- c("unrelatedEU_popstructure", "unrelatedEU_nopopstructure", 
             "relatedEU_nopopstructure")
dirroot <- "~/GWAS/data/LiMMBo/CalibrationOld/samples1000_NrSNP20_Cg"
resultdir <- "~/GWAS/data/LiMMBo/CalibrationOld"

alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

# plotting parameters
xlabel=expression(Expected~~-log[10](italic(p))) 
ylabel= expression(Observed~~-log[10](italic(p)))

axistext <- 8
axistitle <- 8
legendtext <- 8
legendtitle <- 8
striptext <- 8

################
### analysis ### 
################

# Calibration of LMM (with standard REML) for all three kinship structures and 
# moderate number of traits (10,20,30)
calibrationREML <- lapply(kinship, function(K) {
	calibrationH2 <- lapply(genVariance, function(h2) {
        calibrationTraits <- lapply(c(10, 20, 30), function(P) {
            
            directoryClosedForm=paste(dirroot, h2, "_model", model, "/", K, 
                                "/nrtraits", P, "/closedForm",  sep="")
            RML <- readData("lmm_mt", directoryClosedForm, K, 
											"_closedForm")
			if (length(RML) != 0) {
                tmp <- data.frame(RML=RML)
			} else {
				tmp <- NULL
			}
			if (!is.null(tmp)) {
            	names(tmp) <- paste(names(tmp), "_", P, "_", h2, 
                                        "_", K, sep="")
			}
			return(tmp)
        })
        cal_pvalues <- do.call(cbind, rmNulls(calibrationTraits))
    })
    cal_pvalues <- do.call(cbind, calibrationH2)	 
})

calibrationREML <- lapply(calibrationREML, formatPValues)
saveRDS(calibrationREML, paste(resultdir, 
									"/calibrationPvaluesREML.rds", sep=""))

# Calibration of LiMMBo for structured populations and low- to high trait
# numbers (10-100)
calibrationLiMMBo <- lapply("relatedEU_nopopstructure", function(K) {
    calibrationH2 <- lapply(genVariance, function(h2) {
        calibrationTraits <- lapply(c(10,20,30,50,100), function(P) {
            if (P == 10) sampling <- 5
            if (P > 10) sampling <- 10

            directory=paste(dirroot, h2, "_model", model, "/", K, "/nrtraits",
                            P, "/estimateVD/nrtraits_samples", sampling,
                            sep="")
            LiMMBo <- readData("lmm_mt", directory, K, "")
            if (length(LiMMBo) != 0) {
                tmp <- data.frame(LiMMBo=LiMMBo)
            } else {
                tmp <- NULL
            }
            if (!is.null(tmp)) {
                names(tmp) <- paste(names(tmp), "_", P, "_", h2,
                                        "_", K, sep="")
            }
            return(tmp)
        })
        cal_pvalues <- do.call(cbind, rmNulls(calibrationTraits))
    })
    cal_pvalues <- do.call(cbind, calibrationH2)
})

calibrationLiMMBo <- lapply(calibrationLiMMBo, formatPValues)
saveRDS(calibrationLiMMBo, paste(resultdir, 
									"/calibrationPvaluesLiMMBo.rds", sep=""))

# Calibration of simple linear model for all population structures and low- to
# high dimensional trait numbers
calibrationLM <- lapply(kinship, function(K) {
    calibrationH2 <- lapply(genVariance, function(h2) {
        calibrationTraits <- lapply(c(10, 50, 100), function(P) {
            if (P == 10) sampling <- 5
            if (P > 10) sampling <- 10
            
			directory=paste(dirroot, h2, "_model", model, "/", K, "/nrtraits",
                            P, "/estimateVD/nrtraits_samples", sampling,
                            sep="")

            LM <- readData("lm_mt", directory, K, "")
            if (length(LM) != 0) {
                tmp <- data.frame(LM=LM)
            } else {
                tmp <- NULL
            }
            if (!is.null(tmp)) {
                names(tmp) <- paste(names(tmp), "_", P, "_", h2,
                                        "_", K, sep="")
            }
            return(tmp)
        })
        cal_pvalues <- do.call(cbind, rmNulls(calibrationTraits))
    })
    cal_pvalues <- do.call(cbind, calibrationH2)
})

calibrationLM <- lapply(calibrationLM, formatPValues)
saveRDS(calibrationLM, paste(resultdir,
                                    "/calibrationPvaluesLM.rds", sep=""))

### Combine LMM pvalues from REML and LiMMBo and depict in qq-plots
### Figure 3B in publication, Figure 4.4 in thesis
calibrationPvalues <- rbind(calibrationLiMMBo[[1]], calibrationREML[[3]])
saveRDS(calibrationPvalues, paste(resultdir,
                    "/calibrationPvaluesLMMrelatedNoPopstructure.rds", sep=""))
calibrationPvalues <- filter(cal, observed < 0.001)

limits <- c(0, -log10(min(c(calibrationPvalues$expected, 
							calibrationPvalues$observed))))

png(file=paste(resultdir, "/calibrationSummaryQQAll.png", sep=""),  
    height=4, width=5.2, units="in", res=450)

p <- ggplot(data=calibrationPvalues, aes(x=-log10(expected), 
                                         y=-log10(observed)))
p + geom_segment(aes(x = 0, y = 0, xend = limits[2],
                         yend = limits[2]),
                         color="gray10") +
    geom_point(aes(color=as.factor(analysis)),
               size=1) +
	coord_fixed(ylim=limits, xlim=limits) +
    facet_grid(h2 ~ traits, labeller=label_parsed) + 
    scale_color_manual(values=wes_palette(n=5, 
                                          name="Darjeeling", 
                                          type='continuous')[c(2,3,5)],
                       name="Method") +
    xlab(xlabel) +
    ylab(ylabel) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=axistext, angle=0,
                                     hjust=.5, vjust=.5, face="plain"),
          axis.text.y = element_text(colour="black", size=axistext, angle=0,
                                     hjust=0.5, vjust=0.5, face="plain"),  
          axis.title.x = element_text(colour="black", size=axistitle, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          axis.title.y = element_text(colour="black", size=axistitle, angle=90,
                                      hjust=.5, vjust=.5, face="plain"),          
          strip.text = element_text(colour="black", size=striptext, angle=0,
                                    hjust=.5, vjust=.5, face="plain"),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(colour="black", size=legendtext, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.title= element_text(colour="black", size=legendtitle, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.key = element_rect(colour = NA)) 
dev.off()


### Combine pvalues from LM and LMM via LiMMBo and depict in qq-plots
### Table S4 in publication, Figure 4.5 in thesis
calibrationPvaluesCompare <- rbind(dplyr::filter(calibrationLiMMBo[[1]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 5),
									dplyr::filter(calibrationLM[[1]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3),
                                    dplyr::filter(calibrationLM[[2]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3),
                                    dplyr::filter(calibrationLM[[3]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3 ))
saveRDS(calibrationPvaluesCompare, paste(resultdir,
                                    "/calibrationPvaluesLMvsLMM.rds", sep=""))
calibrationPvaluesCompare$analysis <- factor(calibrationPvaluesCompare$analysis,
                                            levels = c("LM", "LiMMBo"))

limits <- c(0, -log10(min(c(calibrationPvaluesCompare$expected, 
							calibrationPvaluesCompare$observed))))


png(file=paste(resultdir, "/calibrationSummaryQQLMMvsLM.png", sep=""),
    height=4, width=5.2, units="in", res=450)

p <- ggplot(data=calibrationPvaluesCompare, aes(x=-log10(expected),
                                         y=-log10(observed)))
p + geom_segment(aes(x = 0, y = 0, xend = limits[2],
                         yend = limits[2]),
                         color="gray10") +
    geom_point(aes(color=analysis, shape=as.factor(kinship)),
               size=1) +
    coord_fixed(ylim=limits, xlim=limits) + 
    scale_shape_manual(values=c(0,1,2), name="Kinship") +
    scale_color_manual(values=wes_palette(n=5,
                                          name="Darjeeling",
                                          type='continuous')[c(4,2)],
                       name="Method") +
    xlab(xlabel) +
    ylab(ylabel) +
    coord_fixed() +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=axistext, angle=0,
                                     hjust=.5, vjust=.5, face="plain"),
          axis.text.y = element_text(colour="black", size=axistext, angle=0,
                                     hjust=0.5, vjust=0.5, face="plain"),
          axis.title.x = element_text(colour="black", size=axistitle, angle=0,
                                      hjust=.5, vjust=0 face="plain"),
          axis.title.y = element_text(colour="black", size=axistitle, angle=90,
                                      hjust=.5, vjust=.5, face="plain"),
          strip.text = element_text(colour="black", size=striptext, angle=0,
                                    hjust=.5, vjust=.5, face="plain"),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(colour="black", size=legendtext, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.title= element_text(colour="black", size=legendtitle, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.key = element_rect(colour = NA))
dev.off()

