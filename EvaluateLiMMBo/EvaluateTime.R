###########################################################
###                                                     ###
### Evaluate run times for variance decomposition with  ###
### REML and LiMMBo                                     ###
###                                                     ###
### Data generated with setupLiMMBo/runtime.sh          ###
###                                                     ###
### Generates Figure 2 (publication)                    ###
###           Figure 4.2 (thesis)                       ###
###                                                     ###
###########################################################

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


complexityLiMMBo <- function(df, model, N) {
    if (model == "sumBS") {
        sum_bs <- df$bsruns * (N * df$sampled^4 + 
					df$sampled^5)
        return(df$ProctimeSumBootstraps/3600 ~ sum_bs)
    }
    if (model == "combineBS") {
        combine_bs <- 0.5 * df$traits^2
        return(df$ProctimeCombineBootstraps/3600 ~ combine_bs) 
    }
    if (model == "both") {
        sum_bs <- df$bsruns * (N * df$sampled^4 + 
					df$sampled^5)
        combine_bs <- df$traits^2 
        return(df$Proctime/3600 ~ sum_bs + combine_bs) 
    }
}

complexityREML <- function(df, N) {
    Proctime <- df$Proctime
    traits <- N * df$traits^4 + df$traits^5
    return(Proctime/3600 ~ traits )
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

# parameters
h2 <- c(0.2, 0.5, 0.8)
kinship <- 'relatedEU_nopopstructure'
nrtraits <- seq(10, 100, 10)
seed <- 1:10

# directories
rootdir <- '~/data/LiMMBo/Calibration'

resultsdir <- sapply(h2, function(h) {
    sapply(seed, function(s) {
        sapply(nrtraits, function(P) {
            if (P == 10) p <- 5
            if (P != 10) p <- 10
            paste(rootdir, "/relatedEU_nopopstructure/samples1000_traits", P, 
                  "_Cg", h, "_modelnoiseBgOnlygeneticBgOnly/seed", s, 
                  "/estimateVD/nrtraits_samples", p, sep="")
        })
    })
})


## get information of number of bootstrap per run
bscount <- read.table(paste(rootdir, "/Boostrap_sampling_schemes.csv", sep=""),
                        row.names=1, header=TRUE, sep=",")
################
### analysis ### 
################

# get run times for LiMMBo runs generated via  ~/LiMMBo/setupLiMMBo/runtime.sh

## run times of variance decomposition for all bootstraps
t_bs <- do.call(rbind, lapply(resultsdir, function(x) {
    if (file.exists(paste(x,"/process_time_bs.csv", sep=""))) {
        tmp <- read.table(paste(x,"/process_time_bs.csv", sep=""), sep=",")
        tmp$analyses <- gsub(
          ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
          , "\\2_seed\\3_nrtraits\\1", x)
        tmp$traits <- gsub(
          ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
          , "\\1", x)
        tmp$population <- "relatedEU_nopopstructure" 
        tmp$h2 <- gsub(
          ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
          , "\\2", x)
        tmp$seed <- as.numeric(gsub(
          ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
          , "\\3", x))
        colnames(tmp)[1] <- c("ProctimeBootstraps")
        tmp$ProctimeBootstraps <- as.numeric(tmp$ProctimeBootstraps)
    } else {
        tmp <- NA
    }
    return(tmp)
}))

t_bs <- t_bs[!is.na(t_bs$ProctimeBootstraps),]
t_bs$traits <- factor(t_bs$traits, levels=seq(10,100,10))

# include all runs where bootstrap size s=10 -> everything but traits=10
ProctimeBootstraps_mean <- melt(acast(dplyr::filter(t_bs, traits != 10),  
                                      population ~ h2, 
                                      value.var="ProctimeBootstraps", mean))
ProctimeBootstraps_sd <- melt(acast(dplyr::filter(t_bs, traits != 10),  
                                    population ~ h2, 
                                    value.var="ProctimeBootstraps", sd))

ProctimeBootstraps_stats <- cbind(ProctimeBootstraps_mean, 
                                  ProctimeBootstraps_sd[,3])
colnames(ProctimeBootstraps_stats) <- c("popStructure", "h2", "mean", "sd")
write.table(ProctimeBootstraps_stats, paste(rootdir, 
                                            "/ProctimeBootstrapsStats.csv",
                                            sep=""), sep=",",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

## run times for variance decomposition via REML 
t_closedform <- do.call(rbind, lapply(resultsdir, function(x) {
	if (file.exists(paste(x,"/process_time_mtSet.csv", sep=""))) {
    	tmp <- read.table(paste(x, "/process_time_mtSet.csv", 
    	                        sep=""), sep=",")
		tmp$analyses <- gsub(
		  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
		  , "\\2_seed\\3_nrtraits\\1", x)
		tmp$traits <- gsub(
		  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
		  , "\\1", x)
		tmp$popStructure <- "relatedEU_nopopstructure"
		tmp$h2 <- gsub(
		  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
		  , "\\2", x)
        tmp$seed <- as.numeric(gsub(
          ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
          , "\\3", x))
    	colnames(tmp)[1] <- c("Proctime")
    	return(tmp)
	}
}))
t_closedform$setup <- "REML"
t_closedform$popStructure <- sapply(as.character(t_closedform$popStructure), 
                                    switchPop)
t_closedform$sampled <- NA
t_closedform$bsruns <- NA
t_closedform$traits <- as.numeric(as.character(t_closedform$traits))
t_closedform$Component <- "REML"
t_closedform <- t_closedform[,c(2:ncol(t_closedform),1)]

## run times for combined bootstrap processing
t_combined_bs <- do.call(rbind, lapply(resultsdir, function(x) {
    if (file.exists(paste(x,"/process_time_summary.csv", sep=""))) {
		tmp <- read.table(paste(x,"/process_time_summary.csv", sep=""), 
                          sep=",")
		tmp <- data.frame(t(tmp[,-1]))
		tmp$Proctime <- sum(tmp)
	} else {
		tmp <- data.frame(t(c(NA, NA)))
		tmp$Proctime <- NA
	}
	tmp$analyses <- gsub(
	  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
	  , "\\2_seed\\3_nrtraits\\1", x)
	tmp$traits <- gsub(
	  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
	  , "\\1", x)
	tmp$popStructure <- "relatedEU_nopopstructure"
	tmp$h2 <- gsub(
	  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
	  , "\\2", x)
	tmp$seed <- as.numeric(gsub(
	  ".*traits(\\d{2,3})_(Cg0\\.\\d{1})_model.*/seed(\\d*)/.*"
	  , "\\3", x))
	colnames(tmp)[1:2] <- c("ProctimeCombineBootstraps", 
								"ProctimeSumBootstraps")
    return(tmp)
    }))

t_combined_bs <- t_combined_bs[!is.na(t_combined_bs$Proctime),]
t_combined_bs$popStructure <- sapply(as.character(t_combined_bs$popStructure), 
                                     switchPop)
t_combined_bs$setup <- "LiMMBo"
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

## Combine REML and LiMMBo results
proctimeComparison <- rbind(t_combined_bs.m, t_closedform)
proctimeComparisonSummary <- summarySE(proctimeComparison, 
									   measurevar="Proctime", 
                					   groupvars=c("traits", "setup", 
													"Component"))

## Fit lines to run times (and predict run times for REML) case
sum_bs <- lm(complexityLiMMBo(t_combined_bs, model="sumBS", N=1000))
combine_bs <- lm(complexityLiMMBo(t_combined_bs, model="combineBS", N=1000))
both <- lm(complexityLiMMBo(t_combined_bs, model="both", N=1000))
reml <- lm(complexityREML(t_closedform, N=1000))
predictedREML <- predict(reml, 
                        data.frame(traits=1000*seq(10, 100, 10)^4 + 
                                       seq(10, 100, 10)^5 ),
                        se.fit=TRUE)

## Plot results 
## Figure 2 LiMMBo paper, Figure 4.2 thesis

### Plotting parameters
color <- wes_palette(5, name="Darjeeling", type='continuous')[c(5, 3, 2, 1)]
textsize <- 12

### Mean, std and lines for four run time components
p <- ggplot()
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
                              Component="REML")
              , aes(x=x, y=y, group=1, colour=Component)) +
    scale_colour_manual(
        values=color, 
        name="Runtimes:",
        labels=c("Combine bootstraps",  "Sum bootstraps", "LiMMBo (total)", 
				"REML"),
        guide=guide_legend(nrow=2)) +
    ylim(c(-5,110)) +
    ylab("Process time [h]") +
    xlab("Traits") +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black", size=textsize, angle=0,
									hjust=.5, vjust=1, face="plain"),
          axis.text.y = element_text(colour="black", size=textsize, angle=0,
                                     hjust=1, vjust=0, face="plain"),  
          axis.title.x = element_text(colour="black", size=textsize, angle=0,
                                      hjust=.5, vjust=0, face="plain"),
          axis.title.y = element_text(colour="black", size=textsize, angle=90,
                                      hjust=.5, vjust=.5, face="plain"),
          legend.text = element_text(colour="black", size=textsize, angle=0, 
                                     hjust=1, vjust=0, face="plain"),
          legend.title = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position ='bottom') 

ggsave(plot=p, file=paste(rootdir, "/proctime.pdf", sep=""), height=4, 
       width=5.2, units="in")

ggsave(plot=p, file=paste(rootdir, "/proctime.eps", sep=""), height=4, 
       width=5.2, units="in")
