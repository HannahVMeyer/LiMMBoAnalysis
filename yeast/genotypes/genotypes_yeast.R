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


############
### data ###
############

directory <- "~/data/LiMMBo/yeast/inputdata"

geno <- read.table(paste(directory, "/BYxRM_GenoData.txt", sep=""), header=TRUE)
load(paste(directory, "/cross.Rdata", sep=""))
genoinfo <- cross$geno

pheno_filtered <- fread(paste(directory, "/BYxRM_pheno_format.txt", sep=""), 
                        sep="\t", 
                        header=TRUE,
                        data.table=FALSE)

################
### Analysis ###
################

### 1. genotypes ###
# a) Filter genotypes for samples that passed phenotype filter
geno_filtered <- geno[,c(TRUE, 
                         (colnames(geno) %in% pheno_filtered[,1])[-1])]

# b) make plink formated files
## i) map file: chr, SNP ID, Centimorgan, bp position
### 1) get centimorgan positions
genomap <- unlist(sapply(genoinfo, function(x) x$map))
names(genomap) <- gsub("\\d{1,2}\\.", "", names(genomap))
genomap <- genomap[match(geno_filtered[,1], names(genomap))]

### 2) get chromosomes
chr <- gsub("(\\d*)_chr0?(\\d{1,2})_(\\d*)_\\w*_\\w*", "\\2", names(genomap))

### 3) get bp position
bp <- gsub("(\\d*)_chr0?(\\d{1,2})_(\\d*)_(\\w*)_(\\w*)", "\\3", names(genomap))

### 4) alleles
a1 <- gsub("(\\d*)_chr0?(\\d{1,2})_(\\d*)_(\\w*)_(\\w*)", "\\4", names(genomap))
a2 <- gsub("(\\d*)_chr0?(\\d{1,2})_(\\d*)_(\\w*)_(\\w*)", "\\5", names(genomap))

map <- data.frame(chr=as.numeric(chr), 
                  snpid=names(genomap), 
                  cm=as.numeric(genomap), 
                  bp=as.numeric(bp))
write.table(map, file=paste(directory, "/BYxRM.map", sep=""), 
            sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

## ii) ped file: FID, IID, ID father, ID mother, sex, pheno, allele calls
### 1) convert B/R to respective alleles
BR2a1a2diploid <- do.call(cbind, lapply(1:length(a1), function(SNP) {
    tmp <- lapply(geno_filtered[SNP,-1], function(call) {
        if (call == "B") return(c(a1[SNP], a1[SNP]))
        if (call == "R") return(c(a2[SNP], a2[SNP]))
        else return(c(NA,NA))
    })
    tmp <- do.call(rbind, tmp)
    return(tmp)
}))

ped <- data.frame(FID=colnames(geno_filtered)[-1], 
                  IID=colnames(geno_filtered)[-1], 
                  father=rep(0, nrow(BR2a1a2diploid)), 
                  mother=rep(0, nrow(BR2a1a2diploid)), 
                  sex=rep(0, nrow(BR2a1a2diploid)), 
                  pheno= rep(0,  nrow(BR2a1a2diploid)),
                  BR2a1a2diploid )
write.table(ped, file=paste(directory, "/BYxRM.ped", sep=""), sep=" ", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)


