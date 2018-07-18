###############################
### Libraries and Functions ###
###############################

library("ggplot2")
library("gridExtra")
library("gplots")
library("reshape2")
library("ggrepel")

ttcooccurence <- function(nrtraits, nrtraitssampled) {
    1/nrtraits * 1/(nrtraits-1) * nrtraitssampled * (nrtraitssampled-1)
}

estimateCounts <- function( nrbootstraps, nrtraits, nrsamples) {
    traits <- seq(1,nrtraits, 1)
    sampling <-  sapply(1:nrbootstraps, function(dummy) 
		sample(traits, nrsamples, replace=FALSE))
    sampling_matrix <- matrix(nc=nrtraits, nr=nrtraits, 0)
    for (i in 1:nrbootstraps) {
        sampling_matrix[sampling[,i], sampling[,i]] <- 
			sampling_matrix[sampling[,i], sampling[,i]] + 1
    }
    diagonal <- diag(sampling_matrix)
    diag(sampling_matrix) <- NA
    return(list(diagonal=diagonal, 
				all= as.vector(sampling_matrix)[!is.na(as.vector(
														sampling_matrix))]))
}
############
### Data ###
############

rootdir <- "~/data/LiMMBo/Calibration"

nrtraits <- seq(10, 100, 10)
nrsamples <- c(3, 5, 10)
it=1000

################
### analysis ###
################

### via combinatorics ### 
cooccurence <- sapply(nrtraits, function(t) {
    sapply(nrsamples, function(s,t){
        if (t ==  s) {
            return(NA)
        } else {
            return(ttcooccurence(t,s))
        }
    }, t=t)
})

colnames(cooccurence) <- paste("NrTraits", nrtraits, sep="") 
rownames(cooccurence) <- paste("NrTraitsSampled:", nrsamples)

# set cooccurence to fixed value: 3
cooccurence_fixed <- 3/cooccurence

# Number of bootstraps fixed to 2000:
bs2000_100_10 <- cooccurence[rownames(cooccurence) == 
                                            "NrTraitsSampled: 10", 
                                            colnames(cooccurence) == 
                                                "NrTraits100"] * 2000

bs2000 <- apply(cooccurence, c(1,2) ,function(co) bs2000_100_10/co)
write.table(round(bs2000), paste(rootdir, "/Boostrap_sampling_schemes.csv", 
                                sep=""),
            quote=FALSE, col.names=NA, row.names=TRUE, sep=",")


### via simulation ###
set.seed(100)

# number of bootstraps
bs=c(20, 40, 60, 80, 100, 150, 300, 500, 1000,1500,2000)

counts <- lapply(nrtraits, function(t) {
    tmp <- lapply(nrsamples, function(s,t){
        if (s < t) {
            tmp_bs <- sapply(bs, estimateCounts, nrtraits=t, nrsamples=s)
            colnames(tmp_bs) <- paste("Nr_bootstraps", bs, sep="")
        } else {
            tmp_bs <- matrix(NA, nc=length(bs), nr=4)
        }
        return(tmp_bs)
    }, t=t)
    names(tmp) <- rownames(cooccurence)
    counts_mean <- do.call(c, lapply(1:length(tmp), function(s) tmp[[s]][1,]))
    counts_min <- do.call(c, lapply(1:length(tmp), function(s) tmp[[s]][2,]))
    counts_max <- do.call(c, lapply(1:length(tmp), function(s) tmp[[s]][3,]))
    counts_sd <- do.call(c, lapply(1:length(tmp), function(s) tmp[[s]][4,]))
    tmp_t <- data.frame(counts = counts_mean, maxcount = counts_max, 
						mincount = counts_min, sdcount = counts_sd, 
						NrTraitsSampled = factor(rep(names(tmp_s), 
													each=length(bs)),  
										levels=paste("Nr of traits sampled:",
													c(3,5,10))), 
						NrBootstraps = bs, 
						NrTraits = as.factor(rep(t,length(bs))))
    return(tmp_t)
})
   

counts <- lapply(nrtraits, function(t) {
    tmp_s <- lapply(nrsamples, function(s,t){
        if (s < t) {
            tmp_bs <- lapply(bs, estimateCounts, nrtraits=t, nrsamples=s)
            colnames(tmp_bs) <- paste("Nr_bootstraps", bs, sep="")
        } else {
            tmp_bs <- matrix(NA, nc=length(bs), nr=4)
        }
        return(tmp_bs)
    }, t=t)
    names(tmp_s) <- rownames(cooccurence)
    counts_diagonal <- do.call(c, lapply(1:length(tmp_s), function(s) tmp_s[[s]][1,]))
    counts_all <- do.call(c, lapply(1:length(tmp_s), function(s) tmp_s[[s]][2,]))
     tmp_t <- data.frame(counts = counts_mean, maxcount = counts_max, mincount = counts_min, sdcount = counts_sd, NrTraitsSampled = factor(rep(names(tmp_s),each=length(bs)),  levels=paste("Nr of traits sampled:",c(3,5,10))), NrBootstraps = bs, NrTraits = as.factor(rep(t,length(bs))))
    return(tmp_t)
})

 
summary_sampling <- do.call(rbind, counts)    



counts <- lapply(nrtraits, function(t, nrsamples ) {
    if (t != 10) {
        tmp_bs <- sapply(rep(t, it), buildSamplingMatrix,  s= nrsamples, 
                         countOnly=TRUE)
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


universe <- function(v,k,t){
    choose(k,t) * v^t
}

blocks <- function(v,k) { v^k}
P=c(50, 100, 200, 300, 400, 500,600, 700, 800, 900, 1000)
U=sapply(P, universe, k=10, t=2)
S=sapply(P, blocks, k=10)

plot(P, U, pch=20)
plot(P, S, pch=20)



