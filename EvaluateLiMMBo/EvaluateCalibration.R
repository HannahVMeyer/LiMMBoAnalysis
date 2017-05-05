library("data.table")
library("ggplot2")
library("wesanderson")
library("dplyr")
library("plyr")


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

N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
model <- "noiseBgOnlygeneticBgOnly"
kinship <- c("unrelatedEU_popstructure", "unrelatedEU_nopopstructure", 
             "relatedEU_nopopstructure")
dirroot <- "/homes/hannah/GWAS/data/LiMMBo/Calibration/samples1000_NrSNP20_Cg"
resultdir <- "/homes/hannah/GWAS/data/LiMMBo/Calibration"
Traits <- seq(10,100,10)

calibrationH2 <- lapply(genVariance, function(h2) {
    calibrationKinship <- lapply(kinship, function(K) {
        calibrationTraits <- lapply(Traits, function(P) {
            if (P==10) sampling <- 5
            if (P>10) sampling <- 10
            
            alpha <- c(5e-5, 5e-6, 5e-7, 5e-8)

            directory=paste(dirroot, h2, "_model", model, "/", K, "/nrtraits", 
                            P, "/estimateVD/nrtraits_samples", sampling, 
                            sep="")
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
                names(tmp_pvalues) <- paste(names(tmp_pvalues), "_", P, "_", h2, 
                                        "_", K, sep="")
            }
            
            if (P <= 30 ) {
                directory=paste(dirroot, h2, "_model", model, "/", K, 
                                "/nrtraits", P, "/closedForm",  sep="")
                cal_closedForm <- 
                    calibrationAcrossChromosomes("lmm_mt", directory=directory,
                                                 estimateString="_closedForm", 
                                                 K=K, P=P)
                if (all(is.na(cal_closedForm[[3]]))) {
                    missing <- NA
                } else {
                    missing <- 
                        paste(cal_closedForm[[3]][!is.na(cal_closedForm[[3]])],
                              collapse=", ")
                }
                tmp_closedForm <- data.frame(alpha=alpha, 
                                             lm=rep(NA, length(alpha)), 
                                             lmm=cal_closedForm[[1]], 
                                             P=rep(P, length(alpha)), 
                                             h2=rep(h2, length(alpha)), 
                                             model=rep(model, length(alpha)),  
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
           #names(tmp) <- paste("Traits", Traits, sep="")
           return(tmp)
        })
        cal_summary <- ldply(calibrationTraits, function(x) x[[1]])
        cal_pvalues <- do.call(cbind, lapply(calibrationTraits, 
                                             function(x) x[[2]]))
        tmp <- list(cal_summary, cal_pvalues)
        #names(tmp) <- kinship
        return(tmp)
    })
    cal_summary <- ldply(calibrationKinship, function(x) x[[1]])
    cal_pvalues <- lapply(calibrationKinship, function(x) x[[2]])
    cal_pvalues <- cal_pvalues[sapply(cal_pvalues, function(x) dim(x)[1] != 0)]
    cal_pvalues <- do.call(cbind, cal_pvalues)
    tmp <- list(cal_summary, cal_pvalues)
    #names(tmp) <- paste("h2_", genVariance, sep="")
    return(tmp)
})

calibrationSummary <- ldply(calibrationH2, function(x) x[[1]])
calibrationPvalues <- do.call(cbind, lapply(calibrationH2, function(x) x[[2]]))
calibrationSummary$alpha <- factor(calibrationSummary$alpha, 
                                  level=c(5e-5, 5e-6, 5e-7, 5e-8))
saveRDS(calibrationSummary, paste(resultdir, "/calibrationSummary.rds", sep=""))
saveRDS(calibrationPvalues, paste(resultdir, "/calibrationPvalues.rds", sep=""))

pdf(file=paste(resultdir, "/calibrationSummary.pdf", sep=""),  height=12, 
    width=12)
p <- ggplot(filter(calibrationSummary, estimate == "LiMMBo"), 
            aes(x=as.factor(P), y=-log10(lm)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha)) + 
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), title="LM") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6]) +
 	geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))

p <- ggplot(filter(calibrationSummary, estimate =="LiMMBo"), 
            aes(x=as.factor(P), y=-log10(lm_pc)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha)) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), 
         title="LM with PCs") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6]) +
 	geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))

p <- ggplot(filter(calibrationSummary, estimate =="LiMMBo"), 
            aes(x=as.factor(P), y=-log10(lmm)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha)) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[5:8]) +
    labs(x = "Number of traits", y = expression(-log[10](FDR)), title="LMM") +
 	geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
 	geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6]) +
 	geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7]) +
 	geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8]) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))

dev.off()

pdf(file=paste(resultdir, "/calibrationSummary_closedForm.pdf", sep=""), 
    height=12, width=12)
p <- ggplot(filter(calibrationSummary, P <= 30), aes(x=as.factor(P), 
                                                    y=-log10(lmm)))
p + geom_bar(stat = "identity", position="dodge",
             aes(fill=alpha, alpha=estimate)) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[c(5:8)]) +
	scale_alpha_manual(values = c(0.4, 0.8)) + 
    labs(x = "Number of traits", y = expression(-log[11](FDR)), title="LMM") +
    geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
    geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6]) +
    geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7]) +
    geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8]) +
    ylim(c(0,10))+
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()
