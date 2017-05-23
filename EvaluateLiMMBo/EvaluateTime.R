###############################
### Libraries and Functions ###
###############################

library("wesanderson")
library("ggplot2")
library("reshape2") # melt, acast

library("dplyr") # filter
library("gplots")

switchPop <- function(x) {
    switch(EXPR=x, 
           "relatedEU_nopopstructure" = "relatedNoPopStructure", 
           "unrelatedEU_popstructure" = "unrelatedPopStructure", 
           "unrelatedEU_nopopstructure" = "unrelatedNoPopStructure")
}


complexityLiMMBo <- function(df, model) {
    if (model == "sumBS") {
        sum_bs <- df$bsruns *(df$sampled^2 + df$sampled^4)
        return(df$ProctimeSumBootstraps/3600 ~ sum_bs)
    }
    if (model == "combineBS") {
        combine_bs <- 0.5*(df$traits^2 + df$traits^4)
        return(df$ProctimeCombineBootstraps/3600 ~ combine_bs) 
    }
    if (model == "both") {
        sum_bs <- df$bsruns *(df$sampled^2 + df$sampled^4)
        combine_bs <- 0.5*(df$traits^2 + df$traits^4)
        return(df$Proctime/3600 ~ sum_bs + combine_bs) 
    }
}

complexityRML <- function(df) {
    Proctime <- df$Proctime
    traits <- df$traits^2 + df$traits^4
    return(Proctime/3600 ~ traits )
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
#t_bs <- readRDS("~/GWAS/data/LiMMBo/Calibration/t_bs.rds")

# include all runs where bootstrap size s=10 -> everything but traits=10
ProctimeBootstraps_mean <- melt(acast(filter(t_bs, traits != 10),  
                                      population ~ h2, 
                                      value.var="ProctimeBootstraps", mean))
ProctimeBootstraps_sd <- melt(acast(filter(t_bs, traits != 10),  
                                    population ~ h2, 
                                    value.var="ProctimeBootstraps", sd))

ProctimeBootstraps_stats <- cbind(ProctimeBootstraps_mean, 
                                  ProctimeBootstraps_sd[,3])
colnames(ProctimeBootstraps_stats) <- c("popStructure", "h2", "mean", "sd")
ProctimeBootstraps_stats$popStructure <- 
    sapply(as.character(ProctimeBootstraps_stats$popStructure), switchPop)

#ProctimeBootstraps_stats <- 
# readRDS("~/GWAS/data/LiMMBo/Calibration/ProctimeBootstrapsStats.rds")

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

p <- ggplot(filter(t_bs, traits != 10), aes(x=population, 
                                            y=ProctimeBootstraps/60, 
                                            fill=as.factor(h2)))
p + geom_boxplot(position=position_dodge()) +
    scale_fill_manual(
        values=wes_palette(4, name="Moonrise2", type='continuous')[1:3], 
        name=expression(h^2)) +
    ylab("Process time [min]") +
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
#t_combined_bs <- readRDS("~/GWAS/data/LiMMBo/Calibration/t_combined_bs.rds")

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
    	    ".*Cg0\\.(\\d{1})_model.*/(.*)/(nrtraits\\d{2,3}).*",
    	    "\\1", x)
    	colnames(tmp)[1] <- c("Proctime")
    	return(tmp)
	}
}))
#t_closedform <- readRDS("~/GWAS/data/LiMMBo/Calibration/t_closedform.rds")
t_closedform$traits <- factor(t_closedform$traits, levels=seq(10,100,10))
t_closedform$setup <- "RML"
t_closedform$popStructure <- sapply(as.character(t_closedform$popStructure), 
                                    switchPop)

proctimeComparison <- rbind(t_combined_bs[,-c(1:2)], t_closedform)
#proctimeComparison <- 
#readRDS("~/GWAS/data/LiMMBo/Calibration/proctimeComparison.rds")

bscount <- read.table(
    "~/GWAS/data/LiMMBo/Calibration/Boostrap_sampling_schemes.csv", 
    row.names=1, header=TRUE, sep=",")
t_combined_bs$sampled <- 10
t_combined_bs$sampled[t_combined_bs$traits == 10] <- 5

t_combined_bs$bsruns <- as.numeric(bscount[3,])
t_combined_bs$bsruns[t_combined_bs$traits == 10] <- bscount[2,1]

t_combined_bs$traits <- as.numeric(as.character(t_combined_bs$traits))
meanP10 <- mean(filter(t_bs, traits != 10)$ProctimeBootstraps)/60

sum_bs <- lm(complexityLiMMBo(t_combined_bs, model="sumBS"))
combine_bs <- lm(complexityLiMMBo(t_combined_bs, model="combineBS"))
both <- lm(complexityLiMMBo(t_combined_bs, model="both"))

t_closedform$traits <- as.numeric(as.character(t_closedform$traits))
rml <- lm(complexityRML(t_closedform))
predictedRML <- predict(rml, 
                        data.frame(traits=seq(10, 100, 10)^2 + 
                                       seq(10, 100, 10)^4 ))
color <- wes_palette(4, name="Moonrise2", type='continuous')[3:1]

pdf(file=paste(rootdir, "/proctimeBootstrap_stats.pdf", sep=""), onefile=TRUE, 
    height=12,width=16, paper =   "a4r")
p <- ggplot(t_combined_bs)
p + geom_boxplot( aes(x=traits, y=ProctimeSumBootstraps/3600, group=traits), 
                  position=position_dodge(), col=color[2]) +
    #geom_point(aes(x=traits, y=ProctimeSumBootstraps/3600), 
     #          position=position_dodge(width=0.9)) +
    geom_line(data=data.frame(y=sum_bs$fitted.values/3600, 
                              x =t_combined_bs$traits,
                              col="Sum bootstraps")
              , aes(x=x, y=y, group=1, color=col)) +
    geom_boxplot( aes(x=traits, y=ProctimeCombineBootstraps/3600, group=traits), 
                  position=position_dodge(), col=color[1]) +
    #geom_point(data=t_combined_bs, aes(x=traits, 
     #                                  y=ProctimeCombineBootstraps/3600), 
      #         position=position_dodge(width=0.9)) +
    geom_line(data=data.frame(y=combine_bs$fitted.values/3600, 
                              x =t_combined_bs$traits,
                              col="Combine bootstraps")
              , aes(x=x, y=y, group=1, color=col)) +  
    geom_point(aes(x=traits, y=Proctime/3600), 
              position=position_dodge(width=0.9), col=color[3]) +
    geom_line(data=data.frame(y=both$fitted.values, 
                              x =t_combined_bs$traits,
                              col="Total")
              , aes(x=x, y=y, group=1, color=col)) +
    scale_colour_manual(
        values=color, 
        name="Runtime components:") +
    #scale_y_log10() +
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

p <- ggplot(t_combined_bs)
p + geom_boxplot( aes(x=bsruns, y=ProctimeSumBootstraps/3600, group=bsruns),
                width = 100, col=color[2])  +
    geom_point(aes(x=bsruns, y=ProctimeSumBootstraps/3600))

p <- ggplot(t_combined_bs)
p + geom_boxplot( aes(x=traits, y=Proctime/3600, group=traits), 
                  position=position_dodge()) +
    geom_point(aes(x=traits, y=Proctime/3600), 
               position=position_dodge(width=0.9)) +
    geom_line(data=data.frame(y=both$fitted.values, 
                              x =t_combined_bs$traits)
              , aes(x=x, y=y, group=1)) +
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

    
    
p <- ggplot(proctimeComparison, aes(x=as.numeric(as.character(traits)), 
                                    y=Proctime/3600, 
                                    fill=as.factor(h2), alpha=as.factor(setup)))
p + geom_bar(stat='identity', position=position_dodge()) +
    scale_fill_manual(values= wes_palette(4, name="Moonrise2", 
                                      type='continuous')[1:3], 
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




