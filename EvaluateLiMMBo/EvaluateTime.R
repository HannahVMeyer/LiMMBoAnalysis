###############################
### Libraries and Functions ###
###############################
library("wesanderson")
library("ggplot2")
library("reshape2") # melt, acast

library("dplyr") # filter
library("plyr")
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

complexityREML <- function(df) {
    Proctime <- df$Proctime
    mtSet <- df$traits^4 + df$traits^5
    traits2 <- df$traits^2
    traits4 <- df$traits^4
    return(Proctime ~ mtSet)
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
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar))
    
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

rootdir <- '~/data/LiMMBo/Calibration'

# directories
h2 <- c(0.2, 0.5, 0.8)
kinship <- c('unrelatedEU_nopopstructure', 
             'unrelatedEU_popstructure', 
             'relatedEU_nopopstructure')

nrtraits_closedform <- seq(10,100,10)
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
            paste(rootdir, "Old/samples1000_NrSNP20_Cg",
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
#t_bs <- readRDS("~/GWAS/data/LiMMBo/CalibrationOld/t_bs.rds")

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


# bootstrap process times combined 
t_combined_bs <- do.call(rbind, lapply(bootstrapdir, function(x) {
    tmp <- read.table(paste(x,"/process_time_summary.csv", sep=""), sep=",")
    tmp <- data.frame(t(tmp[,-1]))
	tmp$ProctimeAll <- sum(tmp)
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

t_combined_bs <- t_combined_bs[!is.na(t_combined_bs$Proctime),]
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
t_closedform$setup <- "REML"
#t_closedform$popStructure <- sapply(as.character(t_closedform$popStructure), 
#                                    switchPop)
t_closedform$sampled <- NA
t_closedform$bsruns <- NA
t_closedform$traits <- as.numeric(as.character(t_closedform$traits))
t_closedform$Component <- "REML"
t_closedform <- t_closedform[,c(2:ncol(t_closedform),1)]

bscount <- read.table(
    "~/GWAS/data/LiMMBo/Calibration/Boostrap_sampling_schemes.csv", 
    row.names=1, header=TRUE, sep=",")
t_combined_bs$sampled <- 10
t_combined_bs$sampled[t_combined_bs$traits == 10] <- 5
t_combined_bs$traits <- as.numeric(as.character(t_combined_bs$traits))
t_combined_bs$bsruns <- sapply(t_combined_bs$traits, function(x) {
    if (x == 10 ) {
        bscount[2,as.numeric(gsub("NrTraits", "", colnames(bscount))) == x]
    } else {
        bscount[3,as.numeric(gsub("NrTraits", "", colnames(bscount))) == x]
    }})

t_combined_bs.m <- melt(t_combined_bs, id.vars = c( "analyses", "traits",
                                                  "popStructure", "h2",
                                                  "setup", "sampled", "seed",
                                                  "bsruns" ),
                      value.name = "Proctime",
                      variable.name = "Component")

proctimeComparison <- rbind(t_combined_bs.m, t_closedform)
proctimeComparisonSummary <- summarySE(proctimeComparison, measurevar="Proctime", 
                groupvars=c("traits", "setup", "Component"))

toFitLiMMBo <- dplyr::filter(t_combined_bs, popStructure == "relatedNoPopStructure")
toFitRML <- dplyr::filter(t_closedform, popStructure == "relatedNoPopStructure")

sum_bs <- lm(complexityLiMMBo(toFitLiMMBo, model="sumBS"))
combine_bs <- lm(complexityLiMMBo(toFitLiMMBo, model="combineBS"))
both <- lm(complexityLiMMBo(toFitLiMMBo, model="both"))
reml <- lm(complexityREML(toFitRML))

predictedREML <- predict(reml, 
                        data.frame(REML=1000*(seq(10, 100, 10)^4 +
                                                  seq(10, 100, 10)^4) + 
                                       seq(10, 100, 10)^5 ),
                        se.fit=TRUE)


color <- wes_palette(5, name="Darjeeling", type='continuous')[c(5, 3, 2, 1)]
pd <- position_dodge() # move them .05 to the left and right
textsize <- 12

p <- ggplot(proctimeComparison)
p <- p + geom_point(data=proctimeComparisonSummary,
                 aes(x=as.factor(traits), y=Proctime/3600, colour=Component)) +
    geom_errorbar(data=proctimeComparisonSummary, 
                  aes(x=as.factor(traits), ymin=(Proctime-sd)/3600, 
                      ymax=(Proctime+sd)/3600,  colour=Component), 
                  width=.1) +
    geom_line(data=data.frame(y=sum_bs$fitted.values, 
                              x = as.factor(toFitLiMMBo$traits),
                              Component="ProctimeSumBootstraps")
              , aes(x=x, y=y, group=1, colour=Component)) +
    geom_line(data=data.frame(y=combine_bs$fitted.values, 
                              x =as.factor(toFitLiMMBo$traits),
                              Component="ProctimeCombineBootstraps")
              , aes(x=x, y=y, group=1, colour=Component)) +
    geom_line(data=data.frame(y=both$fitted.values, 
                              x = as.factor(toFitLiMMBo$traits),
                              Component="Proctime")
              , aes(x=x, y=y, group=1, colour=Component)) +
    geom_line(data=data.frame(y=predictedREML$fit, 
                              x = as.factor(seq(10,100,10)),
                              #y=predictedRML$fit/3600, 
                              #x = as.factor(seq(10,100,10)),
                              Component="REML")
              , aes(x=x, y=y, group=1, colour=Component)) +
    scale_colour_manual(
        values=color, 
        name="Runtimes:",
        labels=c("Combine bootstraps",  "Sum bootstraps", "LiMMBo (total)", "REML"),
        guide=guide_legend(nrow=2)) +
    ylim(c(-5,110)) +
    ylab("Process time [h]") +
    xlab("Traits") +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black",size=textsize,angle=0,hjust=.5,
                                     vjust=1,face="plain"),
          axis.text.y = element_text(colour="black",size=textsize,angle=0,hjust=1,
                                     vjust=0,face="plain"),  
          axis.title.x = element_text(colour="black",size=textsize,angle=0,hjust=.5,
                                      vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=textsize,angle=90,hjust=.5,
                                      vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=textsize,angle=0,hjust=1,
                                     vjust=0,face="plain"),
          legend.title = element_blank(),
          #legend.title= element_text(colour="black",size=16,angle=0,hjust=1,
          #                           vjust=0,face="plain"),
          legend.key = element_rect(colour = NA),
          legend.position ='bottom') 

ggsave(plot=p, file="~/data/LiMMBo/Calibration/proctime_new.pdf", height=4, 
       width=5.2, units="in")

#p + geom_boxplot(data=filter(proctimeComparison, 
#                             popStructure == "relatedNoPopStructure"),
#                 aes(x=as.factor(traits), y=Proctime/3600, colour=Component),
#                 position="identity") +

