###############################
### Libraries and Functions ###
###############################

library("ggplot2")
#library ("grid")
library("gridExtra")
library("gplots")
library("reshape2")
library("ggrepel")

ttcooccurence <- function(nrtraits, nrtraitssampled) {
    1/nrtraits * 1/(nrtraits-1) * nrtraitssampled * (nrtraitssampled-1)
}

buildSamplingMatrix <- function(P, s, minCooccurrence=3, countOnly=FALSE) {
    # determine number of predefinded tuples per sampling matrix
    counts <- matrix(0, ncol=P, nrow=P)
    bootstrap <- NULL
    while (min(counts) < minCooccurrence) {
        it <- sample(1:P, s, replace=FALSE)
        bootstrap <- rbind(bootstrap, it)
        counts[it, it] <- counts[it, it] + 1
    }
    if (countOnly) {
        return(nrow(bootstrap))
    } else {
        return(list(bootstrap, counts))
    }
}



################
### analysis ###
################

### via combinatorics ### 
nrtraits <- seq(10, 100, 10)
nrsamples <- 10

cooccurence <- sapply(nrtraits, function(t) {
    sapply(nrsamples, function(s,t){
        if (t != 10) {
            return(ttcooccurence(t,s))
        } else {
            return(ttcooccurence(t,5))
        }
    }, t=t)
})

colnames(cooccurence) <- paste("NrTraits", nrtraits, sep="") 
rownames(cooccurence) <- paste("NrTraitsSampled:", nrsamples)

# set cooccurence to fixed value: 3
cooccurence_fixed <- 3/cooccurence

# Number of cooccurences in 1000 bootstraps for 100 traits and 10 sampling traits
bs2000_NrTraits100_samples10 <- cooccurence[rownames(cooccurence) == "NrTraitsSampled: 10", colnames(cooccurence) == "NrTraits100"]*2000
#bs2000_NrTraits100_samples20 <- cooccurence[rownames(cooccurence) == "NrTraitsSampled: 10", colnames(cooccurence) == "NrTraits100"]*2000

bs_for_count_eq_18 <- apply(cooccurence, c(1,2) ,function(co) bs2000_NrTraits100_samples10/co)
write.table(round(bs_for_count_eq_18), "~/GWAS/data/LiMMBo/Calibration/Boostrap_sampling_schemes.csv", quote=FALSE, col.names=NA, row.names=TRUE, sep=",")

### via simulation ###
set.seed(100)
it=1000

counts <- lapply(nrtraits, function(t, nrsamples ) {
    if (t != 10) {
        tmp_bs <- sapply(rep(t, it), buildSamplingMatrix,  s= nrsamples, countOnly=TRUE)
    } else {
        tmp_bs <- sapply(rep(t, it), buildSamplingMatrix,  s=5, countOnly=TRUE)
    }
    return(tmp_bs)
}, nrsamples=10)

summary_sampling <- do.call(rbind, counts)  
colnames(summary_sampling) <- 1:it
rownames(summary_sampling) <-  nrtraits
summary_sampling.m <- melt(summary_sampling)
colnames(summary_sampling.m) <- c("NrTraits", "Run", "Count")

bs_count <- data.frame(NrTraits=seq(10, 100, 10), 
                       NrTraitsSampled=as.factor(c(5, rep(10,9))), 
                       bs=c(as.numeric(cooccurence_fixed[2,1]),
                            as.numeric(cooccurence_fixed[3,-1])))

axistext <- 16
axistitle <- 16
legendtext <- 16
legendtitle <- 16

p <- ggplot(summary_sampling.m, aes(x=NrTraits, y=Count))
p + geom_boxplot(aes(group=as.factor(NrTraits))) + 
    #geom_point(data=bs_count, aes(x=NrTraits, y=bs),pch=19,size=2) +
    #geom_text_repel(data=bs_count, aes(x=NrTraits, y=bs, label=round(bs)), size=6) +
    #geom_hline(yintercept=countsNrTraits100_samples10[1,9]) +
    ylab("Number of bootstraps") +
    xlab("Number of traits") +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black",size=axistext,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=axistext,angle=0,hjust=0.5,vjust=0.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=axistitle,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=axistitle,angle=90,hjust=.5,vjust=.5,face="plain"),
          legend.text = element_text(colour="black",size=legendtext,angle=0,hjust=1,vjust=0,face="plain"),
          legend.title= element_text(colour="black",size=legendtitle,angle=0,hjust=1,vjust=0,face="plain"),
          legend.key = element_rect(colour = NA)) 


  






