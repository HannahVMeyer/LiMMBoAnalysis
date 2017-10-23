###########################################################
###                                                     ###
### Evaluate power of genetic association studies for	###
### univariate LMM and multivariate LMMs with LiMMBo    ###
###                                                     ###
### Data generated with setupLiMMBo/power.sh            ###
###                                                     ###
### Generates Figure 4, S7 (publication)                ###
###           Figure 4.6 B1 (thesis)                    ###
###                                                     ###
###########################################################

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
seed <- 1:50
affected <- c(0.01, 0.04, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
model <- "noiseFixedAndBggeneticFixedAndBg"
kinship <- "relatedEU_nopopstructure"
dirroot <- "~/data/LiMMBo/Power"
alpha <- c(0.1, 0.01, 0.001)


# get association statistics for power analyses run through:
# ~/LiMMBo/setupLiMMBo/power.sh

powerH2 <- lapply(genVariance, function(h2) {
    powerTraits <- lapply(Traits, function(P) {
        powerSeed <- lapply(seed, function(S) {
           powerAffected <- lapply(affected, function(a) {
                if (P==10) sampling <- 5
                if (P>10) sampling <- 10

                directory=paste(dirroot, "/samples", N, "_traits", P, "_NrSNP", 
                                NrSNPs, "_Cg", h2, "_model", model, "/seed", S, 
                                "/TraitsAffected", a, "/estimateVD", sep="")

                power <- data.frame(sapply(c("lmm_mt", "lmm_st"), powerAnalysis, 
                                           directory=directory))
                power$affected <- a
                power$alpha <- alpha
                power$Traits <- P
                power$H2 <- h2
                power$seed <- S
                
                return(power)
            })
            pow <- do.call(rbind, powerAffected)
        })
        pow <- do.call(rbind, powerSeed)
    })
    pow <- do.call(rbind, powerTraits)
})

powerAll <- do.call(rbind, powerH2)
powerAll$alpha <- factor(powerAll$alpha, labels=paste("FDR:", 
                                                      unique(powerAll$alpha)))
powerAll$H2 <- factor(powerAll$H2, level=,labels=paste("h[2]:", 
                                                       unique(powerAll$H2)))
saveRDS(powerAll, paste(dirroot, "/powerTraitsAffectedNew.rds", 
						sep=""))

# Display the percentage of detected true SNPs for all parameter combinations
# (Supplementary Figure 7 in LiMMBo paper, figure B1 in thesis) 
powerAll.m <- melt(powerAll, 
                   id.vars=c("alpha", "Traits", "H2", "seed", "affected"),
                   measured.vars=c("lmm_mt", "lmm_st"), 
                   value.name="sigSNPs",
                   variable.name="Model")

powerAll.m <- powerAll.m[!is.na(powerAll.m$sigSNPs),]
color <- wes_palette(5, name="Darjeeling", type='continuous')[2:3]

p <- ggplot(filter(powerAll.m, alpha == "FDR: 0.01"),  aes(x=as.factor(Traits), 
                             y=sigSNPs/NrSNPs*100)) 
p <- p + geom_boxplot(aes(color=as.factor(Model)), outlier.colour = NA,
                 position=position_dodge(width=0.9)) + 
    geom_point(aes(color=as.factor(Model)), 
               position=position_jitterdodge(dodge.width=0.9),
               size=0.8) +
    facet_grid(as.factor(affected) ~ H2, labeller=label_parsed) +
    scale_fill_manual(values = color, name="Model") +
    scale_color_manual(values = color, name="Model") +
    labs(x = "Number of traits", y = "%detected true SNPs") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8),
          strip.background = element_rect(fill="white"))
ggsave(plot=p, file=paste(dirroot, "/powerAll.pdf", sep=""),
       height=12, width=12)
ggsave(plot=powerplot, file=paste(dirroot, "/powerAll.eps", sep=""), 
       height=10, width=10, units="in")

# Display the percentage of detected true SNPs for selected parameter
# combination (Figure 4 in LiMMBo paper; Figure 4.6 in thesis ) 
textsize <- 10

## ideal scenario, with all traits being affected
affected1_h02 <- ggplot(filter(powerAll.m, alpha == "FDR: 0.01", affected == 1, 
                         H2 == "h[2]: 0.2"),  
                  aes(x=as.factor(Traits), 
                      y=sigSNPs/NrSNPs*100))
affected1_h02 <- affected1_h02 + 
    geom_boxplot(aes(fill=Model)) + 
    scale_fill_manual(values = color) +
    labs(x =  "Number of traits", y = "%detected true SNPs") +
    ylim(c(0,52)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=textsize, angle=0,
                                     hjust=.5, vjust=1,face="plain"),
          axis.text.y = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),  
          axis.title.x = element_text(colour="black", size=textsize, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          axis.title.y = element_text(colour="black", size=textsize, angle=90,
                                      hjust=.5, vjust=.5, face="plain"),
          legend.text = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position = 'bottom')

## keep number of traits and genetic architeture constanst, 
## change number of traits affected
P50_h02 <- ggplot(filter(powerAll.m, alpha == "FDR: 0.01", Traits == 50, 
                    H2 == "h[2]: 0.2", affected > 0.1),  
            aes(x=as.factor((100*affected)), 
                y=sigSNPs/NrSNPs*100)) 
P50_h02 <- P50_h02 + 
    geom_boxplot(aes(fill=Model)) + 
    scale_fill_manual(values = color) +
    labs(x = "% affected traits", y = "") +
    ylim(c(0,52)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=textsize, angle=0,
                                     hjust=.5, vjust=1, face="plain"),
          axis.text.y = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),  
          axis.title.x = element_text(colour="black", size=textsize, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          legend.text = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position = 'bottom')

## keep number of traits and affected traits constanst, 
## change genetic background
P100_affected0.6 <- ggplot(filter(powerAll.m, alpha == "FDR: 0.01", 
                                  Traits == 100, affected == 0.6),  
                           aes(x=H2, y=sigSNPs/NrSNPs*100))
P100_affected0.6 <- P100_affected0.6 + 
    geom_boxplot(aes(fill=Model)) + 
    scale_fill_manual(values = color) +
    scale_x_discrete(labels=c(0.2, 0.5, 0.8)) +
    labs(x = expression(h[2]), y = "") +
    ylim(c(0,52)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=textsize, angle=0,
                                     hjust=.5, vjust=1, face="plain"),
          axis.text.y = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),  
          axis.title.x = element_text(colour="black", size=textsize, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          legend.text = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position = 'bottom')

legendModel <- get_legend(P100_affected0.6)
plots <-  plot_grid(affected1_h02  + theme(legend.position ='none'), 
                    P50_h02  + theme(legend.position ='none'), 
                    P100_affected0.6 + theme(legend.position ='none'),
                    align="h", nrow=1, labels=c("A", "B", "C"), 
                    label_size = 12, hjust=0)
powerplot <- plot_grid(plots, legendModel, nrow=2, rel_heights = c(5,1)) 
                          
ggsave(plot=powerplot, file=paste(dirroot, "/power.pdf", sep=""), 
       height=4, width=5.2, units="in")
ggsave(plot=powerplot, file=paste(dirroot, "/power.eps", sep=""), 
       height=4, width=5.2, units="in")


