###############################
### Libraries and Functions ###
###############################

library("wesanderson")
library("ggplot2")
library("reshape2")

library ("grid")
library("gridExtra")
library("gplots")
library("StatMatch")
library(dplyr)


switchPop <- function(x) {
    switch(EXPR=x, 
           "relatedEU_nopopstructure"="relatedNoPopStructure", 
           "unrelatedEU_popstructure" = "unrelatedPopStructure", 
           "unrelatedEU_nopopstructure" = "unrelatedNoPopStructure")
}

############
### data ### 
############

rootdir <- '~/GWAS/data/LiMMBo/Calibration'

# directories
h2 <- c(0.2, 0.5, 0.8)
kinship <- c('unrelatedEU_nopopstructure', 
             'unrelatedEU_popstructure', 
             'relatedEU_nopopstructure')

nrtraits_closedform <- c(10, 20)
nrtraits_bootstrap <- seq(10,100,10)

## closedform directories
closedformdir <- sapply(h2, function(h) {
    sapply(kinship, function(k) {
        sapply(nrtraits_closedform, function(n) {
            paste(rootdir, "/samples1000_NrSNP20_Cg", 
                  h, "_modelnoiseBgOnlygeneticBgOnly/", 
                  k, "/nrtraits", n, "/closedForm", sep="")
        })
    })
})

## boostrap directories
bootstrapdir <- sapply(h2, function(h) {
    sapply(kinship, function(k) {
        sapply(nrtraits_bootstrap, function(n) {
            if (n == 10) p=5
            if (n != 10) p=10
            paste(rootdir, "/samples1000_NrSNP20_Cg",
                  h, "_modelnoiseBgOnlygeneticBgOnly/", 
                  k, "/nrtraits", n, "/estimateVD/nrtraits_samples", 
                  p ,"/seed29348", sep="")
        })
    })
})

################
### analysis ### 
################

# average bootstrap run times
t_bs <- do.call(rbind, lapply(bootstrapdir, function(x) {
    tmp <- read.table(paste(x,"/process_time_bs.csv", sep=""), sep=",")
    tmp$analyses <- gsub(
      ".*(Cg0\\.\\d{1})_model.*/(.*)/nrtraits(\\d{2,3}).*"
      , "\\1_\\2_nrtraits\\3", x)
    tmp$traits <- gsub(
      ".*(Cg0\\.\\d{1})_model.*/(.*)/nrtraits(\\d{2,3}).*"
      , "\\3", x)
    tmp$population <- gsub(
      ".*(Cg0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*"
      , "\\2", x)
    tmp$h2 <- gsub(
      ".*Cg(0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*"
      , "\\1", x)
    colnames(tmp)[1] <- c("ProctimeBootstraps")
    tmp$ProctimeBootstraps <- as.numeric(tmp$ProctimeBootstraps)
    return(tmp)
    }))

t_bs$traits <- factor(t_bs$traits, levels=seq(10,100,10))
saveRDS(t_bs, paste(rootdir, "/t_bs.rds", sep=""))



ProctimeBootstraps_mean <- melt(acast(t_bs,  population ~ h2, 
                                      value.var="ProctimeBootstraps", mean))
ProctimeBootstraps_sd <- melt(acast(t_bs,  population ~ h2, 
                                    value.var="ProctimeBootstraps", sd))

ProctimeBootstraps_stats <- cbind(ProctimeBootstraps_mean, 
                                  ProctimeBootstraps_sd[,3])
colnames(ProctimeBootstraps_stats) <- c("popStructure", "h2", "mean", "sd")
ProctimeBootstraps_stats$popStructure <- 
    sapply(as.character(ProctimeBootstraps_stats$popStructure), switchPop)

saveRDS(ProctimeBootstraps_stats, paste(rootdir, "/ProctimeBootstrapsStats.rds",
                                      sep=""))

pdf(file=paste(rootdir, "/proctimeBootstraps_stats.pdf", sep=""), onefile=TRUE, 
    height=12,width=16, paper =   "a4r")
p <- ggplot(ProctimeBootstraps_stats, aes(x=as.factor(popStructure), y=mean, 
                                          fill=as.factor(h2)))
p + geom_bar(stat='identity', position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                  position=position_dodge(width=0.9), width=0.25) +
	scale_fill_manual(
	    values=wes_palette(4, name="Moonrise2", type='continuous')[1:3], 
	    name=expression(h^2)) +
	ylab("Process time [s]") +
 	xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black",size=16,angle=0,hjust=0.5,
                                     vjust=1,face="plain"),
          axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,
                                     vjust=0,face="plain"),  
          axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,
                                      vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,
                                      vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=16,angle=0,hjust=1,
                                     vjust=0,face="plain"),
          legend.title= element_text(colour="black",size=20,angle=0,hjust=1,
                                     vjust=0,face="plain"),
          legend.key = element_rect(colour = NA)) 
dev.off()


# bootstrap process times combined 
t_combined_bs <- do.call(rbind, lapply(bootstrapdir, function(x) {
    tmp <- read.table(paste(x,"/process_time_summary.csv", sep=""), sep=",")
    tmp <- data.frame(t(tmp[,-1]))
	tmp$Proctime <- sum(tmp)
    tmp$analyses <- gsub(".*(Cg0\\.\\d{1})_model.*/(.*)/nrtraits(\\d{2,3}).*",
                         "\\1_\\2_nrtraits\\3", x)
    tmp$traits <- gsub(".*(Cg0\\.\\d{1})_.*/(.*)/nrtraits(\\d{2,3}).*",
                       "\\3", x)
    tmp$popStructure <- gsub(
        ".*Cg(0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*",
        "\\2", x)
    tmp$h2 <- gsub(".*(Cg0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*",
                   "\\1", x)
    colnames(tmp)[1:2] <- c("ProctimeCombineBootstraps", 
                            "ProctimeSumBootstraps")
    return(tmp)
    }))

t_combined_bs$traits <- factor(t_combined_bs$traits, levels=seq(10,100,10))
t_combined_bs$popStructure <- sapply(as.character(t_combined_bs$popStructure), 
                                     switchPop)
t_combined_bs$setup <- "LiMMBo"
saveRDS(t_combined_bs, paste(rootdir, "/t_combined_bs.rds", sep=""))

# closed form process times
t_closedform <- do.call(rbind, lapply(closedformdir, function(x) {
	if (file.exists(paste(x,"/timeVarianceDecomposition_mtSet.csv", sep=""))) {
    	tmp <- read.table(paste(x, "/timeVarianceDecomposition_mtSet.csv", 
    	                        sep=""), sep=",")
    	tmp$analyses <- gsub(
    	    ".*(Cg0\\.\\d{1})_model.*/(.*)/nrtraits(\\d{2,3}).*", 
    	    "\\1_\\2_nrtraits\\3", x)
    	tmp$traits <- gsub(
    	    ".*(Cg0\\.\\d{1})_model.*/(.*)/nrtraits(\\d{2,3}).*",
    	    "\\3", x)
    	tmp$popStructure <- gsub(
    	    ".*Cg(0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*",
    	    "\\2", x)
    	tmp$h2 <- gsub(
    	    ".*(Cg0\\.\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*",
    	    "\\1", x)
    	colnames(tmp)[1] <- c("Proctime")
    	return(tmp)
	}
}))
t_closedform$traits <- factor(t_closedform$traits, levels=seq(10,100,10))
t_closedform$setup <- "RML"
t_closedform$popStructure <- sapply(as.character(t_closedform$popStructure), 
                                    switchPop)
saveRDS(t_closedform, paste(rootdir, "/t_closedform.rds", sep=""))

proctimeComparison <- rbind(t_combined_bs[,-c(1:2)], t_closedform)
saveRDS(proctimeComparison, paste(rootdir, "/proctimeComparison.rds", sep=""))

pdf(file=paste(rootdir, "/proctimeBootstrap_stats.pdf", sep=""), onefile=TRUE, 
    height=12,width=16, paper =   "a4r")
p <- ggplot(proctimeComparison, aes(x=traits, y=Proctime/3600, 
                                    fill=as.factor(h2), alpha=setup))
p + geom_bar(stat='identity', position=position_dodge()) +
	scale_fill_manual(values=
	                   wes_palette(4, name="Moonrise2", type='continuous')[1:3], 
	                  name=expression(h^2~ ":")) +
	scale_alpha_manual(values=c(0.6, 0.8), name="VD method") +
    geom_point(data=t_combined_bs, aes(x=traits, y=ProctimeSumBootstraps/3600, 
                                       shape="ProctimeSumBootstraps"), 
               position=position_dodge(width=0.9)) +
    geom_point(data=t_combined_bs, aes(x=traits, 
                                       y=ProctimeCombineBootstraps/3600, 
                                       shape="ProctimeCombineBootstraps"), 
               position=position_dodge(width=0.9)) +
	scale_shape_manual(values=c(0,18), name="Analysis step:") +
    facet_grid (popStructure ~ .) +
	ylab("Process time [h]") +
 	xlab("Traits") +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,
                                     vjust=1,face="plain"),
          axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,
                                     vjust=0,face="plain"),  
          axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,
                                      vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,
                                      vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=16,angle=0,hjust=1,
                                     vjust=0,face="plain"),
          legend.title= element_text(colour="black",size=16,angle=0,hjust=1,
                                     vjust=0,face="plain"),
          legend.key = element_rect(colour = NA),
          legend.position ='bottom') 
dev.off()




