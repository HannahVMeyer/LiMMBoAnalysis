###############################
### Libraries and Functions ###
###############################

library("gplots")
require("gtable")
require("ggplot2")
require("grid")
library("gridExtra")
library("cowplot")
library("wesanderson")
library("scales")

sizeLabel=20
sizeText=20

theme_hm <- function() {
    theme(#axis.line = element_blank(),
          axis.title=element_text(size=sizeLabel),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.text=element_text(size=sizeText),
          legend.title=element_text(size=sizeLabel),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid=element_blank())
}

numberTicks <- function(n) {function(limits) pretty(limits, n)}
scaleTicks <- function(x) sprintf("%.2f", x)

############
### data ###
############
#directory="~/data/LiMMBo/simulateData/genotypes"

infile <- snamemake@input
outfile <- snakemake@output

################
### Analysis ###
################

### compute single value decomposition/eigendecomposition of kinship matrix
kinship_svd <- svd(kinship)

### save first 20 PCs as covariates for LM
write.table(data.frame(colnames(kinship), kinship_svd$u[,1:20]),
            paste(gsub(".csv", "_pc.csv", infile), sep=""), 
            sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)

### Depict structure of kinships as split by first two eignevectors
### set color schemes
col_hm <- colorRampPalette(col=c("white", 
                                 wes_palette(16, name="Moonrise2", 
                                             type='continuous')[c(1)]))(80)
col_pop <-  wes_palette(16, name="Moonrise2", type='continuous')[c(13,11, 9,7)]
col_NA <- wes_palette(16, name="Moonrise2", type='continuous')[1]
col_sp <- wes_palette(16, name="Moonrise2", type='continuous')[13]

### create heatmap ggplot objects
kinship_clustering <- hclust(dist(kinship))$order
kinship_4hm <- kinship[kinship_clustering, kinship_clustering]
diag(kinship_4hm) <- NA
kinship_4hm.m <- reshape2::melt(kinship_4hm)

kinship_hm <- ggplot(kinship_4hm.m, aes(Var1, Var2)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradientn(colours = col_hm, na.value=col_NA, 
                         limits=c(min(kinship_4hm, na.rm=TRUE),
                                  max(kinship_4hm, na.rm=TRUE))) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    guides(fill=FALSE) +
    labs(x="Samples", y="Samples") +
    theme_hm()


### create scatterplot ggplot objects
kinship_4sp.m <- data.frame(PC1=kinship_svd$u[, 1], PC2=kinship_svd$u[, 2])
kinship_sp <- ggplot(kinship_4sp.m, aes(PC1, PC2)) + 
    geom_point(color=col_sp, size=0.6)  +
    theme_bw() +
    theme(axis.title.x=element_text(size=sizeLabel),
          axis.title.y=element_text(size=sizeLabel),
          axis.text.x=element_text(size=sizeText),
          axis.text.y=element_text(size=sizeText))


### combine plots
png(paste(directory, "/simulatedCovarianceMatrices_kinship.png", sep=""), 
    width=1500, height=1800)
plot_grid(kinship_hm, kinship_sp, nrow=2, align='hv')
dev.off()




