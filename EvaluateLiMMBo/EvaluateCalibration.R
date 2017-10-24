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

calibrationAcrossChromosomes <- function(gwas, directory, estimateString, K, 
                                         P) {
    gwasdata <- lapply(1:22, function(chr, gwas) {
        cat(P, "\t", chr, "\t", gwas, "\t", estimateString,"\t", K, "\n")
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

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}
############
### data ### 
############

N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
model <- "noiseBgOnlygeneticBgOnly"
kinship <- c("unrelated_popstructure", "unrelated_nopopstructure", 
             "related_nopopstructure")
dirroot <- "~/GWAS/data/LiMMBo/CalibrationOld/samples1000_NrSNP20_Cg"
resultdir <- "~/GWAS/data/LiMMBo/Calibration"
Traits <- seq(10,100,10)

alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

################
### analysis ### 
################
if (! file.exist(paste(resultdir, "/calibrationSummary.rds", sep=""))) {
    calibrationH2 <- lapply(genVariance, function(h2) {
        calibrationKinship <- lapply(kinship, function(K) {
            calibrationTraits <- lapply(Traits, function(P) {
                if (P==10) sampling <- 5
                if (P>10) sampling <- 10
                
                alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)
    
                directory=paste(dirroot, h2, "_model", model, "/", K, 
                                "/nrtraits", P, "/estimateVD/nrtraits_samples", 
                                sampling, sep="")
                cal <- lapply(c("lm_mt_pcs", "lmm_mt"), 
                              calibrationAcrossChromosomes, directory=directory,  
                              estimateString="", K=K, P=P)
                missing <- sapply(cal, function(m) {
                    if (all(is.na(m[[3]]))) return(NA)
                    else paste(m[[3]][!is.na(m[[3]])], collapse=", ")
                })
                tmp_summary <- data.frame(alpha=alpha, lm=cal[[1]][[1]], 
                                  lmm=cal[[2]][[1]],
                                  P=rep(P, length(alpha)), 
                                  h2=rep(h2, length(alpha)),
                                  model=rep(model, length(alpha)),  
                                  kinship=rep(K, length(alpha)), 
                                  estimate=rep("LiMMBo", length(alpha)), 
                                  missing_lm=rep(missing[1], length(alpha)), 
                                  missing_lmm=rep(missing[2], length(alpha)),
                                  stringsAsFactors=FALSE)
    
                tmp_pvalues <- data.frame(lm=cal[[1]][[4]], lmm=cal[[2]][[4]])
                if (K == "relatedEU_nopopstructure") {
                    names(tmp_pvalues) <- paste(names(tmp_pvalues), "_", P, "_", 
                                                h2, "_", K, sep="")
                }
                
                if (P <= 30 ) {
                    directory=paste(dirroot, h2, "_model", model, "/", K, 
                                    "/nrtraits", P, "/closedForm",  sep="")
                    cal_closedForm <- 
                        calibrationAcrossChromosomes("lmm_mt", 
                                                     directory=directory,
                                                     estimateString=
                                                         "_closedForm", 
                                                     K=K, P=P)
                    if (all(is.na(cal_closedForm[[3]]))) {
                        missing <- NA
                    } else {
                        missing <- 
                            paste(cal_closedForm[[3]][!is.na(
                                cal_closedForm[[3]])],
                                  collapse=", ")
                    }
                    tmp_closedForm <- data.frame(alpha=alpha, 
                                                 lm=rep(NA, length(alpha)), 
                                                 lmm=cal_closedForm[[1]], 
                                                 P=rep(P, length(alpha)), 
                                                 h2=rep(h2, length(alpha)), 
                                                 model=rep(model,length(alpha)),  
                                                 kinship=rep(K, length(alpha)), 
                                                 estimate=rep("RML", 
                                                              length(alpha)), 
                                                 missing_lm=rep(missing[1], 
                                                                 length(alpha)), 
                                                 missing_lmm=rep(missing[3], 
                                                                 length(alpha)),
                                                stringsAsFactors=FALSE)
                    tmp_summary <- rbind(tmp_summary, tmp_closedForm)
                } 
               tmp <- list(tmp_summary, tmp_pvalues)
               return(tmp)
            })
            cal_summary <- ldply(calibrationTraits, function(x) x[[1]])
            cal_pvalues <- do.call(cbind, lapply(calibrationTraits, 
                                                 function(x) x[[2]]))
            tmp <- list(cal_summary, cal_pvalues)
            return(tmp)
        })
        cal_summary <- ldply(calibrationKinship, function(x) x[[1]])
        cal_pvalues <- lapply(calibrationKinship, function(x) x[[2]])
        cal_pvalues <- cal_pvalues[sapply(cal_pvalues, 
                                          function(x) dim(x)[1] != 0)]
        cal_pvalues <- do.call(cbind, cal_pvalues)
        tmp <- list(cal_summary, cal_pvalues)
        return(tmp)
    })

    calibrationSummary <- ldply(calibrationH2, function(x) x[[1]])
    calibrationSummary$alpha <- factor(calibrationSummary$alpha, 
                                      level=c(5e-5, 5e-6, 5e-7, 5e-8))
    calibrationSummary$kinship <- factor(gsub("EU", "",
                                              calibrationSummary$kinship), 
                                       level= kinship, 
                                       labels= gsub("_", "~", kinship))
    calibrationSummary$h2 <- factor(calibrationSummary$h2, 
                                    levels=unique(calibrationSummary$h2),
                                    labels=paste("h[2]:", 
                                                 unique(calibrationSummary$h2)))
    saveRDS(calibrationSummary, paste(resultdir, "/calibrationSummary.rds", 
                                      sep=""))
    
    calibrationPvalues <- do.call(cbind, lapply(calibrationH2, 
                                                function(x) x[[2]]))
    saveRDS(calibrationPvalues, paste(resultdir, "/calibrationPvalues.rds", 
                                      sep=""))
} else {
    calibrationSummary <- readRDS(paste(resultdir, "/calibrationSummary.rds", 
                                        sep=""))
    if (!file.exists(paste(resultdir, "/calibrationPvaluesOrderedMelt.rds", 
                           sep=""))) {
        if (!file.exists(paste(resultdir, "/calibrationPvaluesOrdered.rds", 
                               sep=""))) {
            calibrationPvalues <- readRDS(paste(resultdir, 
                                                "/calibrationPvalues.rds", 
                                                sep=""))
            analysis <- gsub("([lmm]{2,3})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", 
                             "\\1", names(calibrationPvalues))
            traits <- gsub("([lmm]{2,3})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", 
                           "\\2", names(calibrationPvalues))
            h2 <- gsub("([lmm]{2,3})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", "\\3", 
                       names(calibrationPvalues))
            kinship <- gsub("([lmm]{2,3})_(\\d{2,3})_(0\\.\\d)_(.*)EU_(.*)", 
                            "\\4_\\5", 
                            names(calibrationPvalues))
            
            calibrationPvaluesOrdered <- apply(calibrationPvalues, 2, 
                                               function(x) x[order(x)] )
            saveRDS(calibrationPvaluesOrdered, 
                    paste(resultdir, "/calibrationPvaluesOrdered.rds", sep=""))
        } else {
            calibrationPvaluesOrdered <- readRDS(
                paste(resultdir, "/calibrationPvaluesOrdered.rds", sep=""))
        }
        expected <- ppoints(nrow(calibrationPvaluesOrdered))
        subset=nrow(calibrationPvaluesOrdered)
        calibrationPvalues <- melt(calibrationPvaluesOrdered, 
                                   value.name="observed")
        calibrationPvalues$analysis <- rep(analysis, each=subset)
        calibrationPvalues$traits <- rep(traits, each=subset)
        calibrationPvalues$h2 <-rep(h2, each=subset)
        calibrationPvalues$kinship <- rep(kinship, each=subset)
        calibrationPvalues$expected <- expected[1:subset]
        calibrationPvalues$kinship <- factor(calibrationPvalues$kinship, 
                                             level= kinship, 
                                             labels= gsub("_", "~", kinship))
        calibrationPvalues$h2 <- factor(calibrationPvalues$h2, 
                                        levels=unique(calibrationPvalues$h2),
                                        labels=paste("h[2]:", 
                                                unique(calibrationPvalues$h2)))
        
        saveRDS(calibrationPvalues, paste(resultdir, 
                                          "/calibrationPvaluesOrderedMelt.rds", 
                                          sep=""))
    } else {
        calibrationPvalues <- readRDS(paste(resultdir, 
              "/calibrationPvaluesOrderedMelt.rds", 
              sep=""))
    }
}
calibrationSummary.m <- melt(calibrationSummary[,-c(6,9,10)], 
                             id.vars=c("alpha", "P","h2","kinship", "estimate"),
                             value.name="FWER",
                             variable.name="Model")
calibrationSummary.m <- filter(calibrationSummary.m, 
                               alpha %in% c("5e-05","5e-08"),
                               kinship == "related~nopopstructure",
                               P %in% c(10, 50, 100))
calibrationSummary.table <- summarySE(calibrationSummary.m, measurevar="FWER", 
                                      groupvars=c( "P", "Model", "alpha"),
                                      na.rm=TRUE)
calibrationSummary.table <-  dcast(calibrationSummary.table,  
                                   P + alpha ~ Model, value.var="FWER")

write.table(calibrationSummary.table, 
            paste(resultdir, "/calibrationSummaryTable.csv", sep=""), 
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")

# plotting parameters
xlabel=expression(Expected~~-log[10](italic(p))) 
ylabel= expression(Observed~~-log[10](italic(p)))

axistext <- 18
axistitle <- 18
legendtext <- 18
legendtitle <- 18
striptext <- 20

# Filter for subset of traits
calibrationPvalues <- filter(calibrationPvalues, traits %in% c(10,50,100))

# qq plot of lmm and lm for related samples
png(file=paste(resultdir, "/calibrationSummaryQQ.png", sep=""),  height=1000, 
    width=1000)

p <- ggplot(data=calibrationPvalues, aes(x=-log10(expected), 
                                         y=-log10(observed)))
p + geom_point(aes(color=factor(traits, levels=c(10, 50, 100)), 
                   shape=as.factor(analysis)),
               size=3) +
    facet_grid(h2 ~ kinship, labeller=label_parsed) + 
    scale_color_manual(values=wes_palette(n=5, 
                                          name="Darjeeling", 
                                          type='continuous')[c(2,3,5)],
                       name="Traits") +
    scale_shape(name="Model") +
    geom_abline(intercept=0, slope=1, col="black") +
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

# bar plot of FDR estimates for LM for realted samples, unrelated samples with
# and without population structure

png(file=paste(resultdir, "/calibrationSummaryPerModelLM.png", sep=""), 
    height=1000, width=1000)
p <- ggplot(filter(calibrationSummary, 
                   estimate == "LiMMBo", 
                   P %in% c(10,50,100)),
            aes(x=as.factor(P), y=-log10(lm)))
p +  geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5], 
                linetype="dashed") +
    geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6],  
               linetype="dashed") +
    geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7],  
               linetype="dashed") +
    geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8],  
               linetype="dashed") +
    geom_bar(stat = "identity", position="dodge", aes(fill=alpha), alpha=0.8) + 
    facet_grid( h2 ~ kinship, labeller=label_parsed) +
    scale_fill_manual(values = moonrise[5:8], name="FDR") +
    labs(x = "Number of traits", y = expression(-log[10](FDR))) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=axistext, angle=0,
                                     hjust=.5, vjust=.5, face="plain"),
          axis.text.y = element_text(colour="black", size=axistext, angle=0,
                                     hjust=0.5, vjust=0.5, face="plain"),  
          axis.title.x = element_text(colour="black", size=axistitle, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          axis.title.y = element_text(colour="black", size=axistitle, angle=90,
                                      hjust=.5, vjust=.5, face="plain"),          
          strip.text.x = element_text(colour="black", size=striptext.x, angle=0,
                                    hjust=.5, vjust=.5, face="plain"),
          strip.text.y = element_text(colour="black", size=striptext.y,
                                    hjust=.5, vjust=.5, face="plain"),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(colour="black", size=legendtext, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.title= element_text(colour="black", size=legendtitle, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          legend.key = element_rect(colour = NA))
dev.off()

# boxplot of FDR estimates for LMM with covariannce matrices derived from 
# RML or limmbo across all cohorts
pdf(file=paste(resultdir, "/calibrationSummary_closedForm.pdf", sep=""), 
    height=12, width=12)
p <- ggplot(filter(calibrationSummary, 
                   P <= 30,
                   alpha %in% c(5e-5, 5e-8)), aes(x=as.factor(P), 
                                                    y=-log10(lmm)))
p + geom_boxplot(position=position_dodge(width=1), 
                 aes(fill=alpha, alpha=estimate)) +
    scale_fill_manual(values = moonrise[c(5:8)], name="FDR threshold") +
	scale_alpha_manual(values = c(0.4, 0.8), name="VD method") + 
    labs(x = "Number of traits", y = expression(-log[10](FDR))) +
    geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
    geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[6]) +
    ylim(c(4,8)) +
    guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),
                                                         100, 0, 
                                                         alpha=c(0.4,0.8)),
                                                colour=NA))) +
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

