###############################
### Libraries and Functions ###
###############################

library("data.table")
library("ggplot2")
library("wesanderson")
library("plyr")
library("dplyr")

library("reshape2") # melt, acast


rmNulls <- function(compList) {
                nonNulls <- compList[!sapply(compList, is.null)]
}
moonrise <- c(wes_palette(n=4, name="Moonrise1"), 
              wes_palette(n=4, name="Moonrise2"))

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

readData <- function(gwas, directory, K, estimateString) {
	gwasdata <- lapply(1:22, function(chr, gwas) {
        cat(chr, "\t", gwas, "\t",directory,"\t", K, "\n")
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




calibrationAcrossChromosomes <- function(gwas, directory, estimateString, K, 
                                         P) {
    gwasdata <- lapply(1:22, function(chr, gwas) {
        cat(chr, "\t", gwas, "\t",directory,"\t", K, "\n")
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
    missing <- unlist(sapply(gwasdata, function(p) p[[2]]))

    cal_pvalues <- calibration(pvalues)
    if (K == "relatedEU_nopopstructure" && !grepl("closed", estimateString)) {
        return(list(cal_pvalues, length(pvalues), missing, pvalues))
    } else {
        return(list(cal_pvalues, length(pvalues), missing, NULL))
    
    }
}

qq <- function(data, title=NULL, 
              xlabel=expression(Expected~~-log[10](italic(p))), 
              ylabel= expression(Observed~~-log[10](italic(p))), facet=FALSE) {
    p <- ggplot(data=data, aes(x=-log10(expected), y=-log10(observed)))
    p <- p + geom_point(aes(color=traits , shape=analysis)) +
        geom_abline(intercept=0,slope=1, col="black") +
        geom_hline(aes(yintercept=-log10(5e-8))) +
        xlim(c(0,max(data$expected))) +
        ylim(c(0,max(data$observed))) +
        xlab(xlabel) +
        ylab(ylabel) +
        #labs(title=title, x=xlabel, y=ylabel) + 
        theme_bw()
    if (facet) {
        p <- p + facet_grid( h2 ~ kinship)
    }
    p
}

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

	calOrdered <- apply(cal, 2, function(x) x[order(x)])
	expected <- ppoints(nrow(calOrdered))

	cal <- melt(calOrdered, value.name="observed")[,-1]
	names(cal)[1] <- "type"
	cal$analysis <- rep(analysis, each=nrow(calOrdered))
	cal$traits <- rep(traits, each=nrow(calOrdered))
	cal$h2 <-rep(h2, each=nrow(calOrdered))
	cal$kinship <- rep(kinship, each=nrow(calOrdered))
	cal$expected <- expected
	
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
#Traits <- seq(10,100,10)
Traits <- seq(10, 50, 10)

alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

################
### analysis ### 
################

calibrationRML <- lapply(kinship, function(K) {
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

calibrationRML <- lapply(calibrationRML, formatPValues)
saveRDS(calibrationRML, paste(resultdir, 
									"/calibrationPvaluesRML.rds", sep=""))


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

calibrationPvalues <- rbind(calibrationLiMMBo[[1]], calibrationRML[[3]])
saveRDS(calibrationPvalues, paste(resultdir,
                                    "/calibrationPvaluesLMMrelatedNoPopstructure.rds", sep=""))
calibrationPvalues <- filter(cal, observed < 0.001)

xlabel=expression(Expected~~-log[10](italic(p))) 
ylabel= expression(Observed~~-log[10](italic(p)))

axistext <- 8
axistitle <- 8
legendtext <- 8
legendtitle <- 8
striptext <- 8

limits <- c(0, -log10(min(calibrationPvalues$expected, 
							calibrationPvalues$observed)))

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

traitString <- paste("Traits: 10", "Traits: 50", "Traits: 100")

calibrationPvaluesCompare <- rbind(filter(calibrationLiMMBo[[1]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 5),
									filter(calibrationLM[[1]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3),
                                    filter(calibrationLM[[2]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3),
                                    filter(calibrationLM[[3]],
										  as.numeric(h2) == 3, 
                                          as.numeric(traits) == 3 ))
calibrationPvaluesCompare <- filter(calibrationPvaluesCompare, 
							as.numeric(h2) == 3, 	
							observed < 0.0001)


png(file=paste(resultdir, "/calibrationSummaryQQLMMvsLM.png", sep=""),
    height=4, width=5.2, units="in", res=450)

p <- ggplot(data=calibrationPvaluesCompare, aes(x=-log10(expected),
                                         y=-log10(observed)))
p + geom_segment(aes(x = 0, y = 0, xend = limits[2],
                         yend = limits[2]),
                         color="gray10") +
    geom_point(aes(color=as.factor(analysis), shape=as.factor(kinship)),
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

