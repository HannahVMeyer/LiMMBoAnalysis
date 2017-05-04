###############################
### Libraries and Functions ###
###############################

library("data.table")
library("ggplot2")
library("wesanderson")
library("dplyr")
library("reshape2") # melt, acast

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
    return(list(cal_pvalues, length(pvalues), missing))
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
dirroot <- "~/GWAS/data/LiMMBo/Calibration/samples1000_NrSNP20_Cg"
resultdir <- "~/GWAS/data/LiMMBo/Calibration"
Traits <- seq(10,100,10)

alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

################
### analysis ### 
################

calibrationH2 <- lapply(genVariance, function(h2) {
    calibrationKinship <- lapply(kinship, function(K) {
        calibrationTraits <- lapply(Traits, function(P) {
            if (P==10) sampling <- 5
            if (P>10) sampling <- 10
            
            alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

            directory=paste(dirroot, h2, "_model", model, "/", K, "/nrtraits", 
                            P, "/estimateVD/nrtraits_samples", sampling, 
                            sep="")
            cal <- lapply(c("lm_mt_pcs", "lm_mt", "lmm_mt"), 
                          calibrationAcrossChromosomes, directory=directory,  
                          estimateString="", K=K, P=P)
            missing <- sapply(cal, function(m) {
                if (all(is.na(m[[3]]))) return(NA)
                else paste(m[[3]][!is.na(m[[3]])], collapse=", ")
            })
            tmp <- data.frame(alpha=alpha, lm_pc=cal[[1]][[1]], 
                              lm=cal[[2]][[1]], lmm=cal[[3]][[1]],
                              P=rep(P, length(alpha)), 
                              h2=rep(h2, length(alpha)),
                              model=rep(model, length(alpha)),  
                              kinship=rep(K, length(alpha)), 
                              estimate=rep("estimateVD", length(alpha)), 
                              missing_lm_pc=rep(missing[1], length(alpha)), 
                              missing_lm=rep(missing[2], length(alpha)), 
                              missing_lmm=rep(missing[3], length(alpha)))
            
            if (P <= 30 && FALSE) {
                directory=paste(dirroot, h2, "_model", model, "/", K, 
                                "/nrtraits", P, "/closedForm",  sep="")
                cal_closedForm <- 
                    calibrationAcrossChromosomes("lmm_mt", directory=directory,
                                                 estimateString="_closedForm", 
                                                 K=K)
                if (all(is.na(cal_closedForm[[3]]))) {
                    missing <- NA
                } else {
                    missing <- 
                        paste(cal_closedForm[[3]][!is.na(cal_closedForm[[3]])],
                              collapse=", ")
                }
                tmp_closedForm <- data.frame(alpha=alpha, 
                                             LM=rep(NA, length(alpha)), 
                                             sLM=rep(NA, length(alpha)), 
                                             LMM=cal_closedForm[[1]], 
                                             P=rep(P, length(alpha)), 
                                             h2=rep(h2, length(alpha)), 
                                             model=rep(model, length(alpha)),  
                                             kinship=rep(K, length(alpha)), 
                                             estimate=rep("closedForm", 
                                                          length(alpha)), 
                                             missing_lm_pc=rep(missing[1], 
                                                               length(alpha)), 
                                             missing_lm=rep(missing[2], 
                                                            length(alpha)), 
                                             missing_lmm=rep(missing[3], 
                                                             length(alpha)))
                tmp <- rbind(tmp, tmp_closedForm)
            } 
           return(tmp)
        })
        cal <- do.call(rbind, calibrationTraits)
    })
    cal <- do.call(rbind, calibrationKinship)
})

calibrationBgOnly <- do.call(rbind, calibrationH2)
calibrationBgOnly$alpha <- factor(calibrationBgOnly$alpha, 
                                  level=c(5e-5, 5e-6, 5e-7, 5e-8))
calibrationBgOnly$kinship <- factor(calibrationBgOnly$kinship, 
                                  level=kinship)
saveRDS(calibrationBgOnly, paste(resultdir, "/calibrationBGOnly_woRML.rds", sep=""))
calibrationBgOnly$LMM [calibrationBgOnly$LMM == 0 ] <- 5e-8 

axistext <- 16
axistitle <- 16
legendtext <- 16
legendtitle <- 16

pdf(file=paste(resultdir, "/calibrationBGOnly.pdf", sep=""),  height=12, 
    width=12)
p <- ggplot(filter(calibrationBgOnly, estimate == "LiMMBo",
                   alpha %in% c(5e-5, 5e-8)), 
            aes(x=as.factor(P), y=-log10(lm)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha), alpha=0.8) + 
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), title="LM") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[6]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))

p <- ggplot(filter(calibrationBgOnly, estimate =="LiMMBo",
                   alpha %in% c(5e-5, 5e-8)), 
            aes(x=as.factor(P), y=-log10(lm_pc)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha), alpha=0.8) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), 
         title="LM with PCs") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[6]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))

p <- ggplot(filter(calibrationBgOnly, estimate =="LiMMBo",
                   alpha %in% c(5e-5, 5e-8)), 
            aes(x=as.factor(P), y=-log10(lmm)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha), alpha=0.8) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), title="LMM") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[6]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()

calibrationBgOnly.m <- melt(calibrationBgOnly[,1:9], 
                            measure.vars=c("LM", "LMM", "sLM"),
                            value.name="fdr",
                            variable.name="analysis")
col <- wes_palette(n=8, name="Moonrise2", type = 'continuous')

pdf(file=paste(resultdir, "/calibrationBGOnly_combined.pdf", sep=""), 
    height=12, width=12)
p <- ggplot(filter(calibrationBgOnly.m, estimate == "LiMMBo",
                   alpha %in% c(5e-8),
                   analysis != "sLM"), 
            aes(x=as.factor(P), y=-log10(fdr)))

p + geom_bar(stat = "identity", position="dodge", aes(fill=analysis), 
             alpha=0.8) + 
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = col[c(1,5)], name="Model") +
    labs(x = "Number of traits", y = expression(-log[10](FDR))) +
    geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[6]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()

pdf(file=paste(resultdir, "/calibrationBGOnly_closedForm.pdf", sep=""), 
    height=12, width=12)
p <- ggplot(filter(calibrationBgOnly, P <= 30, !is.na(lmm), 
                   alpha %in% c(5e-5, 5e-8)), aes(x=as.factor(P), 
                                                    y=-log10(lmm)))
p + geom_boxplot(position=position_dodge(width=1), aes(fill=alpha, alpha=estimate)) +
    #facet_grid(kinship ~ h2) +
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
    theme(axis.text.x = element_text(colour="black",size=axistext,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=axistext,angle=0,hjust=0.5,vjust=0.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=axistitle,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=axistitle,angle=90,hjust=.5,vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=legendtext,angle=0,hjust=1,vjust=0,face="plain"),
          legend.title= element_text(colour="black",size=legendtitle,angle=0,hjust=1,vjust=0,face="plain"),
          legend.key = element_rect(colour = NA)) 
dev.off()

