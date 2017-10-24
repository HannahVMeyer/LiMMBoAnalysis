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
seedLiMMBo <- 29348
alpha <- c(0.1, 0.01, 0.001)

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
powerAll$alpha <- factor(powerAll$alpha, labels=paste("FDR: ", 
                         unique(powerAll$alpha)))
powerAll$H2 <- factor(powerAll$H2, level=,labels=paste("h[2]:", 
             unique(powerAll$H2)))
saveRDS(powerAll, paste("~/GWAS/data/LiMMBo/Power/powerAll.rds", sep=""))

powerAll.m <- melt(powerAll, 
                   id.vars=c("alpha", "Traits", "H2", "seed"),
                   measured.vars=c("lmm_mt", "lmm_st"), 
                   value.name="sigSNPs",
                   variable.name="Model")

powerAll.m <- powerAll.m[!is.na(powerAll.m$sigSNPs),]
color <- wes_palette(5, name="Darjeeling", type='continuous')[1:2]

pdf(file="~/GWAS/data/LiMMBo/Power/powerAll.pdf",  height=12,width=12)
p <- ggplot(powerAll.m,  aes(x=as.factor(Traits), 
                                                   y=sigSNPs/NrSNPs*100)) 
p + geom_boxplot(aes(fill=Model)) + 
    facet_grid(alpha ~ H2, labeller=label_parsed) +
    scale_fill_manual(values = color) +
    labs(x = "Number of traits", y = "%detected true SNPs") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
dev.off()