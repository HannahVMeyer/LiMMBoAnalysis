###############################
### Libraries and Functions ###
###############################
library("data.table")
library("ggplot2")
library("wesanderson")
library("dplyr")


powerAnalysis <- function(gwas, directory,  fdr=0.001, 
                          alpha=c(0.1, 0.01, 0.001),trueSNPs) {
        p_file <- paste(directory, "/", gwas, "_pempirical_causalSNPs", fdr, 
                        ".csv", sep="")
        if(file.exists(p_file)) {
            p_data <- fread(p_file, data.table=FALSE, stringsAsFactors=FALSE, 
                            header=TRUE)
            detected <- sapply(alpha, function(thr) 
                               length(which(p_data < thr)))
        } else {
            detected <- rep(NA, length(alpha))
        }
}
################
### analysis ###
################

N=1000
NrSNPs=20
genVariance <- c(0.2, 0.5, 0.8)
Traits <- c(10, 50, 100)
seed <- 1:100
model <- "noiseFixedAndBggeneticFixedAndBg"
kinship <- "relatedEU_nopopstructure"
dirroot <- "/homes/hannah/GWAS/data/LiMMBo/Power"

dir=/homes/hannah/GWAS/data/LiMMBo/Power
seedLiMMBo=29348
model=noiseFixedAndBggeneticFixedAndBg

powerH2 <- lapply(genVariance, function(h2) {
    powerTraits <- lapply(Traits, function(P) {
        powerSeed <- lapply(seed, function(S) {
            
            if (P==10) sampling <- 5
            if (P>10) sampling <- 10

            directory=paste(dirroot, "/samples", N, "_traits", P, "_NrSNP", 
                            NrSNPs, "_Cg", h2, "_model", model, "/seed", S, 
                            "/estimateVD", sep="")

            power <- data.frame(sapply(c("lmm_mt", "lmm_st"), powerAnalysis, 
                                       directory=directory))
            power$alpha <- alpha
            power$Traits <- P
            power$H2 <- h2
            power$seed <- S
            
            return(power)
            })
        pow <- do.call(rbind, powerSeed)
    })
    pow <- do.call(rbind, powerTraits)
})

powerAll <- do.call(rbind, powerH2)
powerAll$alpha <- factor(powerAll$alpha, level=alpha)
saveRDS(powerAll, paste("~/GWAS/data/LiMMBo/Power/powerAll.rds", sep=""))

pdf(file="~/GWAS/data/LiMMBo/Power/powerAll.pdf",  height=12,width=12)
p <- ggplot(filter(powerAll.m, alpha==0.001),  aes(x=Traits, 
                                                   y=value/NrSNPs*100)) 
p + geom_bar(stat = "identity", position="dodge",aes(fill=variable)) + 
    facet_grid(kinship ~ H2) +
    scale_fill_manual(values = moonrise[3:8]) +
    labs(x = "Number of traits", y = "%detected true SNPs", title="0.001") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()

p <- ggplot(filter(calibrationBgOnly, estimate =="estimateVD"), 
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

p <- ggplot(filter(calibrationBgOnly, estimate =="estimateVD"), 
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

pdf(file="~/GWAS/data/LiMMBo/Calibration/calibrationBGOnly_closedForm.pdf", 
    height=12, width=12)
p <- ggplot(filter(calibrationBgOnly, P <= 30), aes(x=as.factor(P), 
                                                    y=-log10(lmm)))
p + geom_bar(stat = "identity", position="dodge",aes(fill=alpha, 
                                                     alpha=estimate)) +
    facet_grid(kinship ~ h2) +
    scale_fill_manual(values = moonrise[c(5:8)]) +
	scale_alpha_manual(values = c(0.4, 0.8)) + 
    labs(x = "Number of traits", y = expression(-log[10](FDR)), title="LMM") +
    geom_hline(yintercept = -log10(alpha[1]), colour = moonrise[5]) +
    geom_hline(yintercept = -log10(alpha[2]), colour = moonrise[6]) +
    geom_hline(yintercept = -log10(alpha[3]), colour = moonrise[7]) +
    geom_hline(yintercept = -log10(alpha[4]), colour = moonrise[8]) +
    ylim(c(0,10))+
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()
