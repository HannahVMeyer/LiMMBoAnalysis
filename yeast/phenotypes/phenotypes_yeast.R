###########################################################
###                                                     ###
### Understand, impute and format  yeast quantitative   ###
### growth traits (Bloom et al (2013))      	        ###
###                                                     ###
###     * missing traits, missingness patterns,         ###
###         trait correlations					        ###
###     * testing for missingness mechanism             ###
###         ** Missing completely at random				###
### 			 (Little's test for MCAR)          		###
###    		** Missing at random by prediction patterns	###
###     * imputation via MICE				            ###
###                                                     ###
### Generates Figure  S8-S10 (publication)              ###
###           Figure 5.1-5.4 (thesis)             		###
###                                                     ###
###########################################################

###############################
### Libraries and Functions ###
###############################

library(corrplot)
library(Hmisc)
library(mice)
library(VIM)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(BaylorEdPsych) #LittleMCAR
source("~/projects/utils/corrplotHM.R")

assignPrior <- function(cor_matrix, threshold) {
        corr_info <- cor_matrix
        corr_info[abs(corr_info) < threshold] <- 0
        corr_info[abs(corr_info) >= threshold] <- 1   
        return(corr_info)
}

corrThreshold <- function(threshold, cor_matrix) {
    corr_info <- assignPrior(cor_matrix, threshold)
    pheno_corr_thr <- apply(corr_info, 1,function(x) length(which(x ==1)))
    pheno_no_corr <- names( pheno_corr_thr)[which( pheno_corr_thr == 1)]
    corr_info[which(rownames(corr_info) %in% pheno_no_corr),] <- 1
    diag(corr_info) <- 0
    return(corr_info)
}

cor2matrices <- function(mat1, mat2) {
    cor_r <- sapply(1:ncol(mat1), function(x) rcorr(mat1[,x], mat2[,x])$r[1,2])
    cor_p <- sapply(1:ncol(mat1), function(x) rcorr(mat1[,x], mat2[,x])$P[1,2])
    cor_padjust <-  p.adjust(cor_p)
    return(data.frame(cor_r=cor_r, cor_p=cor_p, cor_padjust= cor_padjust))
}


combineImpute <- function(imputelist, npheno) { 
    mat <- complete(imputelist, 'r')
    nimpute <- ncol(mat)/npheno
    sets <- seq(1, ncol(mat), nimpute)
    imputemedian <- sapply(sets, function(s, mat) {
        apply(mat[, s:(s+nimpute-1)], 1, median)
    }, mat=mat)
    return(imputemedian)
}

Teff <- function(test) {
    # 1. get correlation matrix
    corr_matrix <-  cor(test,  method="spearman")
    # 2. Get eigenvalues of correlation matrix:
    eigenval <- eigen(corr_matrix, only.value=TRUE, symmetric=TRUE)$values
    # 3. Determine effective number of tests:
    t <- sum(sqrt(eigenval))^2/sum(eigenval)
    return(t)
}


############
### data ###
############

directory <- "~/data/LiMMBo/yeast/inputdata"

col=wes_palette(5, name="Darjeeling", type='continuous')[c(4,5)]
col_corr <- colorRampPalette(col=c(col[1], "white",col[2]))(100)

geno <- read.table(paste(directory, "/BYxRM_GenoData.txt", sep=""), header=TRUE)
load(paste(directory, "/cross.Rdata", sep=""))
pheno <- cross$pheno
genoinfo <- cross$geno

################
### Analysis ###
################

# 1. pattern of missing data
## a) distribution: Figure S8.A (publication), Figure 5.1A (thesis)
pdf(paste(directory, "/missing_data_pattern.pdf", sep=""), height=5, width=7)
aggr_plot <- aggr(pheno, col=col, prop=TRUE, numbers=FALSE, sortVars=TRUE, 
                  labels=names(pheno), cex.axis=.7, gap=3, 
                  ylab=c("Combinations"), combined=TRUE, border=NA,
                  bars=TRUE, cex.numbers=0.9, oma = c(8,4,4,2) + 0.1,
                  only.miss=FALSE, ylim=c(0, 0.42))
dev.off()

frequency_missingness <- data.frame(missing=apply(pheno, 2, function(x) 
    length(which(is.na(x)))/length(x)))
frequency_missingness$complete <- 1 - frequency_missingness$missing

per_sample_missingness <- data.frame(missing=apply(pheno, 1, function(x) 
    length(which(is.na(x)))/length(x)))
per_sample_missingness$complete <- 1 - per_sample_missingness$missing

# b) Missing data mechanism
### MCAR
Samples2Keep <- per_sample_missingness$missing <= 0.20
MCARtest <- LittleMCAR(pheno[Samples2Keep,])
pMCAR <- pchisq(MCARtest$chi.square, MCARtest$df, lower.tail=FALSE)

### MAR: Figure 5.2 (thesis)
pheno_missing <- as.data.frame(abs(is.na(pheno)))
corrMiss_P <- t(sapply(1:ncol(pheno), function(y) {
    sapply(1:ncol(pheno), function(x) {
        rcorr(missing[!is.na(pheno[,y]),x], pheno[!is.na(pheno[,y]),y])$P[1,2]
    })
}))
corrMiss_R <- t(sapply(1:ncol(pheno), function(y) {
    sapply(1:ncol(pheno), function(x) {
        rcorr(missing[!is.na(pheno[,y]),x], pheno[!is.na(pheno[,y]),y])$r[1,2]
    })
}))
corrMiss_Padjust <- apply(corrMiss_P, 2, p.adjust)
colnames(corrMiss_R) <- colnames(pheno)
rownames(corrMiss_R) <- colnames(pheno)

pdf(paste(directory, "/correlation_missingness.pdf", sep=""), 
    height=7, width=7)
par(xpd=TRUE)
corrplotHM(corrMiss_R, tl.col='black',method ="color", col=col_corr,
           cl.lim = c(-max(abs(corrMiss_R), na.rm=TRUE),
                      max(abs(corrMiss_R), na.rm=TRUE)), cl.cex=0.7, tl.cex=0.7,
           tl.offset=0.2, cl.offset=0.2, cl.align.text='l', is.corr=FALSE,
           na.label="square", na.label.col="grey",
           mar = c(4,6,7,0), addgrid.col = 'grey')
mtext(text="Phenotype", side=2, line=3, cex=1.1, adj=0.55)
mtext(text="Missingness", side=1, line=2, cex=1.1, adj=0.5)
dev.off()

# 2. dataset with no missing values
noNA_samples <- !apply(pheno, 1, function(x) any(is.na(x)))
pheno_noNA <- pheno[noNA_samples,]

## a) correlation between phenotypes: Figure S9 (publication),
## Figure 5.3 (thesis)
pheno_noNA_cor <- rcorr(as.matrix(pheno_noNA), type="pearson")$r
pheno_noNA_p <- rcorr(as.matrix(pheno_noNA), type="pearson")$P
pheno_noNA_n <- diag(rcorr(as.matrix(pheno_noNA), type="pearson")$n)
pheno_noNA_padjust <- apply(pheno_noNA_p, 1, p.adjust)

pdf(paste(directory, "/correlation_pheno_noNA.pdf", sep=""), 
    height=7, width=7)
corrplot(pheno_noNA_cor, tl.col='black', method ="ellipse", col=col_corr,
         order="hclust", insig="blank", p.mat=pheno_noNA_padjust, addrect=7, 
         tl.cex=0.7, cl.cex=0.7, tl.offset=0.2, cl.offset=0.2, 
         cl.align.text='l')
dev.off() 


# 3. generate matrix with artificial missingness
## a) all samples with at least one and less then 20% missing phenotypes
soi <- intersect(which(per_sample_missingness$missing < 0.20), 
                 which(!noNA_samples))
pheno_NA <- pheno[soi,]

## b) generate matrix of size of pheno_noNA and insert NAs with similar sample/
## pheno missingness distribution
set.seed(34221)
pheno_small <- pheno[sample(1:nrow(pheno), nrow(pheno_noNA)),]

## c) introduce missingness in the fully-phenotyped dataset
pheno_addNA <- pheno_noNA
pheno_addNA[is.na(pheno_small)] <- NA 

frequency_missingness_sample <- data.frame(
    missing=apply(pheno_addNA, 2, function(x) 
        length(which(is.na(x)))/length(x)))
lm_freq <- lm(log(frequency_missingness_sample$missing) ~ 
                  log(frequency_missingness$missing))
plot(-log(frequency_missingness$missing), 
     -log(frequency_missingness_sample$missing), pch=20, 
     ylab=expression(-log[10](missingness[sample])), 
     xlab= expression(-log[10](missingness[all])))

## d) plot dataset: Figure S8.B (publication), Figure 5.1B (thesis)
pdf(paste(directory, "/missing_data_pattern_simulated.pdf", sep=""), 
    height=5, width=7)
aggr_plot <- aggr(pheno_addNA, col=col, prop=TRUE,numbers=FALSE, sortVars=TRUE, 
                  labels=names(pheno), cex.axis=.7, gap=3, 
                  ylab=c("Combinations"), combined=TRUE, border=NA,
                  bars=TRUE, cex.numbers=0.9, oma = c(8,4,4,2) + 0.1,
                  only.miss=FALSE,ylim=c(0, 0.5))
dev.off()

# 4. impute artifically created missing data 
## a) create predictor matrices based in correlations of phenotypes: 
## if trait-trait corr > threshold, use as predictor
## if no correlation greater then threshold use all traits
## design: rows are targets, 0/1 in columns specify whether trait is used as 
## predictor or not

corr_info0.0 <- quickpred(pheno_addNA, mincor=0.0, minpuc=0.2)
corr_info0.1 <- quickpred(pheno_addNA, mincor=0.1, minpuc=0.2)
corr_info0.2 <- quickpred(pheno_addNA, mincor=0.2, minpuc=0.2)
corr_info0.3 <- quickpred(pheno_addNA,mincor=0.3, minpuc=0.2)


## b) impute with different predictor matrix schemes
imputeData_Corr0.0 <- mice(pheno_addNA,m=20, predictorMatrix=corr_info0.1, 
                           maxit=30, meth='pmm', seed=500)
imputeData_Corr0.1 <- mice(pheno_addNA, m=20, predictorMatrix=corr_info0.1, 
                           maxit=30,meth='pmm',seed=500) 
imputeData_Corr0.2 <- mice(pheno_addNA, m=20, predictorMatrix=corr_info0.2, 
                           maxit=30,meth='pmm',seed=500) 
imputeData_Corr0.3 <- mice(pheno_addNA, m=20, predictorMatrix=corr_info0.3, 
                           maxit=30,meth='pmm',seed=500) 

complete_Corr0.0 <-combineImpute(imputeData_Corr0.0,npheno=ncol(pheno_noNA))
complete_Corr0.1 <-combineImpute(imputeData_Corr0.1,npheno=ncol(pheno_noNA))
complete_Corr0.2 <-combineImpute(imputeData_Corr0.2,npheno=ncol(pheno_noNA))
complete_Corr0.3 <-combineImpute(imputeData_Corr0.3,npheno=ncol(pheno_noNA))

cor_Corr0.0_noNA_cor <- cor2matrices(complete_Corr0.0, pheno_noNA)
cor_Corr0.1_noNA_cor <- cor2matrices(complete_Corr0.1, pheno_noNA)
cor_Corr0.2_noNA_cor <- cor2matrices(complete_Corr0.2, pheno_noNA)
cor_Corr0.3_noNA_cor <- cor2matrices(complete_Corr0.3, pheno_noNA)

meanCorr0.0 =mean(cor_Corr0.0_noNA_cor$cor_r)
meanCorr0.1 = mean(cor_Corr0.1_noNA_cor$cor_r)
meanCorr0.2 = mean(cor_Corr0.2_noNA_cor$cor_r)
meanCorr0.3 = mean(cor_Corr0.3_noNA_cor$cor_r)

medianCorr0.0 =median(cor_Corr0.0_noNA_cor$cor_r)
medianCorr0.1 = median(cor_Corr0.1_noNA_cor$cor_r)
medianCorr0.2 = median(cor_Corr0.2_noNA_cor$cor_r)
medianCorr0.3 = median(cor_Corr0.3_noNA_cor$cor_r)

cor_setups <- data.frame(pheno = factor(colnames(pheno_noNA)), 
                         corrAll=cor_Corr0.0_noNA_cor$cor_r, 
                         corrCorr0.1= cor_Corr0.1_noNA_cor$cor_r, 
                         corrCorr0.2=cor_Corr0.2_noNA_cor$cor_r, 
                         corrCorr0.3= cor_Corr0.3_noNA_cor$cor_r)
cor_setups_melt <- melt(cor_setups, variable.name="setup")
cor_setups_melt$type <- gsub("corr(.*)", "\\1", cor_setups_melt$setup)
cor_setups_melt$x <- 1

# 5. Impute full data set: Figure S10 (publication),
## Figure 5.4(thesis)
## a) get best predictors for each trait 
predictors <- colnames(cor_setups)[2:5][apply(cor_setups[,2:5], 1, which.max)]

## b) Filter phenos that cannot be imputed
cutoff<- 0.95
Traits2Keep <- sapply(1:length(predictors), function(x) {
    cor_setups[x,colnames(cor_setups) == predictors[x]] > cutoff
})


xaxis_color <- rep('black', ncol(pheno_noNA))
xaxis_color[!Traits2Keep] <- 'darkred'
text_size <- 12
rect_left <- seq(0.5, 45.5, 2)
rectangles <- data.frame(xmin = rect_left, 
                         xmax = rect_left + 1, 
                         ymin = 0.84, 
                         ymax = 1.003)

p <- ggplot()
p <- p + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
              fill='gray80', alpha=0.8) +
    geom_point(data=cor_setups_melt, aes(x=as.numeric(as.factor(pheno)), 
                                         y=value, color=type), size=1, 
               position=position_dodge(width=0.5)) +
    scale_color_manual(values=wes_palette(4, name="Darjeeling", 
                                          type='continuous'), name='Predictors') +
    theme_bw() +
    scale_y_continuous(limits=c(0.84, 1.003), expand = c(0,0) ) +  
    scale_x_continuous(breaks = seq(1, 46, 1), limits=c(0.5, 46.5),
                       minor_breaks=seq(0.5, 46.5, 1),
                       labels = colnames(pheno_addNA), expand = c(0,0) ) +
    theme(axis.title.y = element_text(size=text_size + 2),
          axis.title.x = element_text(size=text_size + 2),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12, 
                                     color=xaxis_color),
          axis.text.y = element_text(size=text_size),
          legend.title = element_text(size=text_size + 2),
          legend.text = element_text(size=text_size),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(colour = 'grey', size = 0.5)) +
    labs(x="Phenotype", y="Pearson Correlation") +
    geom_hline(yintercept=cutoff)

ggsave(plot=p, file=paste(directory,  
                          "/imputation_correlation_median_imputationvalue.pdf", 
                          sep=""), height=5, width=8, units="in")

Samples2Keep <- per_sample_missingness$missing <= 0.20
pheno_filtered <- pheno[Samples2Keep, Traits2Keep]

## c) desgin predictor matrix for Traits2Kepp
predictors2Keep <- predictors[Traits2Keep]
corr_info0.0_2Keep <- corr_info0.0[Traits2Keep, Traits2Keep]
corr_info0.1_2Keep <- corr_info0.1[Traits2Keep, Traits2Keep]
corr_info0.2_2Keep <- corr_info0.2[Traits2Keep, Traits2Keep]
corr_info0.3_2Keep <- corr_info0.3[Traits2Keep, Traits2Keep]

predictorMatrix <- do.call(rbind, lapply(1:ncol(pheno_filtered), function(x) {
    if (predictors2Keep[x] == "corrAll") {
        tmp <- rep(1, ncol(pheno_filtered))
        tmp[x] <- 0
        return(tmp)
    }
    if (predictors2Keep[x] == "corrCorr0.1") {
        return(corr_info0.1_2Keep[x,])
    }
    if (predictors2Keep[x] == "corrCorr0.2") {
        return(corr_info0.2_2Keep[x,])
    }
    if (predictors2Keep[x] == "corrCorr0.3") {
        return(corr_info0.3_2Keep[x,])
    }
}))

## d) impute missing data
imputeData_filtered <- mice(pheno_filtered, m=20, maxit=30, 
                            predictorMatrix=predictorMatrix, 
                            meth='pmm', seed=500) 

complete_filtered <-combineImpute(imputeData_filtered, 
                                  npheno=ncol(pheno_filtered))
colnames(complete_filtered) <- colnames(pheno_filtered)

write.table(complete_filtered, paste(directory, 
                                     "/BYxRM_pheno_format.txt", sep=""), 
            sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

# 6. compute effective number of tests
effective_number_ofTests <-  Teff(as.matrix(complete_filtered))
write.table(c("effective_number_ofTests:", effective_number_ofTest), 
            paste(directory, "/BYxRM_pheno_effectiveNumberofTests.txt", sep=""),
            sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

save.image(paste(directory, "/phenotypes_yeast.Rdata", sep=""))
