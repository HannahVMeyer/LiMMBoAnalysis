## Libraries ####
library(data.table)
library(VIM)
library(Hmisc)
library(corrplot)

## directories ####
rawdir <- "~/data/LiMMBo/rat/rawdata"
directory <- "~/data/LiMMBo/rat/processeddata"

## data ####
phenotypes_normal <- read.table(paste(rawdir, "/phenotypes_normal.csv", sep=""), 
                            header=TRUE, row.names=1, stringsAsFactors=FALSE, 
                            sep=",")
phenotypes_notNormal <- read.table(paste(rawdir, "/phenotypes_notNormal.csv", 
                                        sep=""),
                              header=TRUE, row.names=1, stringsAsFactors=FALSE, 
                              sep=",")
col <- c('#fc8d62','#8da0cb')

## 1. pattern of missing data ####
pdf(paste(directory, "/missing_data_pattern.pdf", sep=""), height=30, width=15)
aggr_plot <- aggr(phenotypes_normal, col=col, prop=TRUE, numbers=TRUE, 
                  sortVars=TRUE, 
                  labels=FALSE, cex.axis=.7, gap=3, 
                  ylab=c("Frequency","Pattern"), combined=TRUE, border=NA,
                  bars=FALSE, cex.numbers=0.9, oma = c(10,4,4,2) + 0.1,
                  only.miss=FALSE, ylim=c(0, 0.42))
dev.off()

frequency_missingness <- data.frame(missing=
                                        apply(phenotypes_normal, 2, function(x)
                                            length(which(is.na(x)))/length(x)))
frequency_missingness$complete <- 1 - frequency_missingness$missing

per_sample_missingness <- data.frame(missing=
                                         apply(phenotypes_normal, 1, function(x) 
                                             length(which(is.na(x)))/length(x)))
per_sample_missingness$complete <- 1 - per_sample_missingness$missing

Samples2Keep <- per_sample_missingness$missing <= 0.20
Traits2Keep <- frequency_missingness$missing <= 0.20

pheno <- pheno[Samples2Keep, Traits2Keep]

pdf(paste(directory, "/missing_data_pattern.pdf", sep=""), height=30, width=15)
aggr_plot <- aggr(pheno, col=col, prop=TRUE, numbers=TRUE, 
                  sortVars=TRUE, 
                  labels=FALSE, cex.axis=.7, gap=3, 
                  ylab=c("Frequency","Pattern"), combined=TRUE, border=NA,
                  bars=FALSE, cex.numbers=0.9, oma = c(10,4,4,2) + 0.1,
                  only.miss=FALSE, ylim=c(0, 0.42))
dev.off()

# 2. dataset with no missing values
noNA_samples <- !apply(pheno, 1, function(x) any(is.na(x)))
pheno_noNA <- pheno[noNA_samples,]

## a) correlation between phenotypes 
pheno_noNA_cor <- rcorr(as.matrix(pheno_noNA), type="spearman")$r
pheno_noNA_p <- rcorr(as.matrix(pheno_noNA), type="spearman")$P
pheno_noNA_n <- diag(rcorr(as.matrix(pheno_noNA), type="spearman")$n)
pheno_noNA_padjust <- apply(pheno_noNA_p, 1, p.adjust)

pdf(paste(directory, "/", strftime(Sys.time(), "%Y%m%d"), 
          "_correlation_pheno_noNA.pdf", sep=""), height=20, width=15)
col_corr <- colorRampPalette(col=c('#f1a340','#f7f7f7','#998ec3'))(100)
corrplot(pheno_noNA_cor, tl.col='black', method ="ellipse", col=col_corr,
         order="hclust", insig="blank", p.mat=pheno_noNA_padjust, addrect=7, 
         tl.cex=1.2, cl.cex=1.2,tl.offset=0.2, cl.offset=0.2, cl.align.text='l')
dev.off() 

# 3. generate matrix with artificial missingness
## a) all samples with at least one and less then 20% missing phenotypes
pheno_NA <- pheno[which(!noNA_samples),]

## b) generate matrix of size of pheno_noNA and insert NAs with similar sample/
## pheno missingness distribution
set.seed(3422)
pheno_small <- pheno[sample(1:nrow(pheno), nrow(pheno_noNA)),]

## c) introduce missingness in the fully-phenotyped dataset
pheno_addNA <- pheno_noNA
pheno_addNA[is.na(pheno_small)] <- NA 

frequency_missingness_sample <- data.frame(
    missing=apply(pheno_addNA, 2, function(x) 
        length(which(is.na(x)))/length(x)))
frequency_missingness_original <- frequency_missingness[
    rownames(frequency_missingness) %in% rownames(frequency_missingness_sample),]

lm_freq <- lm(frequency_missingness_sample$missing ~ 
                  frequency_missingness_original$missing)
plot(frequency_missingness_original$missing, 
     frequency_missingness_sample$missing, pch=20, 
     ylab=expression(-log[10](missingness[sample])), 
     xlab= expression(-log[10](missingness[all])))

## d) plot dataset
pdf(paste(directory, "/", strftime(Sys.time(), "%Y%m%d"),
          "_missing_data_pattern_simulated.pdf", sep=""), height=20, width=15)
aggr_plot <- aggr(pheno_addNA, col=col, prop=TRUE,numbers=TRUE, sortVars=TRUE, 
                  labels=names(pheno), cex.axis=.7, gap=3, 
                  ylab=c("Frequency","Pattern"), combined=FALSE, border=NA,
                  bars=TRUE, cex.numbers=0.9, oma = c(10,4,4,2) + 0.1,
                  only.miss=FALSE,ylim=c(0, 0.5))
dev.off()

