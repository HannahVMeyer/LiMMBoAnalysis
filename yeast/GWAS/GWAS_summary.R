###############################
### Libraries and Functions ###
###############################
library("data.table")
library("ggplot2")
library("wesanderson")
library("cowplot")

source("~/LiMMBo/manhattanplot.R")

# estimate number of independent tests (according to Galwey, Genetic 
# Epidemiology 2009)
Teff <- function(test) {
    # 1. get correlation matrix
    corr_matrix <-  cor(test,  method="spearman")
    # 2. Get eigenvalues of correlation matrix:
    eigenval <- eigen(corr_matrix, only.value=TRUE, symmetric=TRUE)$values
    eigenval[eigenval < 0 ] <- 0
    # 3. Determine effective number of tests:
    t <- sum(sqrt(eigenval))^2/sum(eigenval)
}

# compute distance based on 'angle between vectors'
dotproddist <- function(x) {
    tmp <- sapply(1:nrow(x), function(l) {
        sapply(1:nrow(x), function(r) {
            x[l,] %*% x[r,]
        })
    })
    return(as.dist(tmp))
}

# filter SNPs with lowest p-value from LD block
filterLD <- function(fac, snp, ld, sigSNP, p_original, multitrait=FALSE) {
    if (multitrait) {
        ldtrait  <- ld[[fac]]
        ptrait <-  p_original[[fac]]
    } else {
        ldtrait  <- ld
        ptrait <-  p_original
    }
    snp_vec <- c(unlist(strsplit(ldtrait[snp,"TAGS"], split="|", fixed=TRUE)),  
                 ldtrait[snp,"SNP"])
    if (any(snp_vec %in%  "NONE")) {
        min_ld_snp <- ptrait[which(ptrait$SNP %in% ldtrait[snp,"SNP"]),]
    }
    if (! any(snp_vec %in%  "NONE")) {
        ld_snps <- ptrait[which(ptrait$SNP %in% snp_vec),]
        ### if nrow(ld_snps) > 1): snps which are in LD with SNPS in question 
        ### are significantly associated with phenotype as well
        ### find minimun of all ld SNPs and only keep SNP with minimum p-value
        if (nrow(ld_snps) > 1) {
            min_ld_snp_tmp <- ld_snps[which.min(ld_snps$P),]
            ### make sure that this SNP has pmin out of all SNPs it's in LD with
            snp_vec2 <- c(unlist(strsplit(
                ldtrait[which(ldtrait[,"SNP"] == min_ld_snp_tmp$SNP),"TAGS"], 
                split="|", fixed=TRUE)), 
                ldtrait[which(ldtrait[,"SNP"] == min_ld_snp_tmp$SNP), "SNP"])
            ld_snps2 <- ptrait[which(ptrait$SNP %in% snp_vec2),]
            min_ld_snp <- ld_snps2[which.min(ld_snps2$P),]
        } else {
            min_ld_snp <- ld_snps
        }
    }
    return( min_ld_snp)
}


############
### data ###
############

directory='~/GWAS/data/LiMMBo/feasabilityBootstrap/yeast'

effective_number_ofTests=32.76497

# estimated in 
fdr_multi=1.99537104453e-05
fdr_single=1.30007477173e-05

genotypes <-  fread(paste(directory, "/inputdata/BYxRM.gen", sep=""),sep=" ", 
                    data.table=FALSE)
genotypes <- data.frame(genotypes[, 1:5], sapply(
    seq(1, ncol(genotypes[, -c(1:5)]),3), function(pos) 
        genotypes[,pos+6] + 2*genotypes[,pos+7]))

psingle <- fread(paste(directory, "/GWAS/pvalues_lmm_single_min_adjust.csv", 
                       sep=""), sep=",", data.table=FALSE, header=TRUE)
colnames(psingle) <- c("ID", "P")

psingle$SNP <- ld$SNP
psingle$CHR <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\1", psingle$ID))
psingle$BP <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\2", psingle$ID))
psingle$TYPE <- 'singletrait'

bsingle <- fread(paste(directory, "/GWAS/betas_lmm_single.csv", sep=""), sep=","
                 , data.table=FALSE)
bsingle$SNP <- ld$SNP

pany <- fread(paste(directory, "/GWAS/pvalues_lmm_any.csv", sep=""), sep=",", 
              data.table=FALSE)
colnames(pany) <- c("ID", "P")
pany$SNP <- ld$SNP
pany$CHR <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\1", pany$ID))
pany$BP <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\2", pany$ID))
pany$TYPE <- 'multitrait'

bany <- fread(paste(directory, "/GWAS/betas_lmm_any.csv", sep=""), sep=",", 
              data.table=FALSE)
bany$SNP <- ld$SNP

################
### analysis ###
################

# Filter SNPs based on pruned dataset
kb <- c(3, 5, 8, 10, 12, 15, 20, 30, 50, 100, 0)
sig_tags <- sapply(kb, function(k) {
    if (k != 0) {
        pruned <- fread(paste(directory, "/inputdata/BYxRM.", k,"kb.pruned.bim", 
                          sep=""), data.table=FALSE)
        psingle_pruned <- psingle[which(psingle$SNP %in% pruned[,2]),]
        pany_pruned <- pany[which(pany$SNP %in% pruned[,2]),]

        psingle_pruned_sig <- psingle_pruned[psingle_pruned$P < fdr_single,]
        pany_pruned_sig <- pany_pruned[pany_pruned$P < fdr_single,]
        return(c(NrSNPs=nrow(pruned),
                 multitrait=nrow(pany_pruned_sig), 
                 singletrait=nrow(psingle_pruned_sig)))
    } else {
        return(c(NrSNPs=nrow(pany),
             multitrait=nrow(pany[pany$P< fdr_single, ]), 
             singletrait=nrow(psingle[psingle$P< fdr_single, ])))
    }
})

sig_tags <- rbind(sig_tags, sig_tags[2,]/sig_tags[3,])
colnames(sig_tags) <- paste("ld_pruned", kb, "kb", sep="")
colnames(sig_tags)[which(kb == 0)] <- "All SNPs"
rownames(sig_tags)[4] <- "multitrait:singletrait"

# chose 15kb filter to display marker SNPs
pruned <- fread(paste(directory, "/inputdata/BYxRM.", 15, "kb.pruned.bim", 
                      sep=""), data.table=FALSE)
psingle_pruned <- psingle[which(psingle$SNP %in% pruned[,2]),]
pany_pruned <- pany[which(pany$SNP %in% pruned[,2]),]
psingle_pruned_sig <- psingle_pruned[psingle_pruned$P < fdr_single,]
pany_pruned_sig <- pany_pruned[pany_pruned$P < fdr_single,]

p_pruned_sig <- rbind(psingle_pruned_sig, pany_pruned_sig)
p_pruned_sig$MARKER <- "tag"

# filter pvalues based on LD - retain SNP with lowest pvalue from ld block (3kb)
# LD window estimate from Liti et al (2009) Population genomics of domestic and 
# wild yeasts
ld <- fread(paste(directory, "/inputdata/BYxRM.3kb.tags.list", sep=""),sep=" ", 
data.table=FALSE)  
ld$ID <- gsub("\\d*_chr0?(\\d{1,2})_(\\d*)_.*", "chr\\1:\\2",ld$SNP)

pany_ldfiltered <- do.call(rbind, lapply(seq_along(1:nrow(pany)), function(snp) 
    filterLD(snp=snp, ld=ld, p_original=pany)))
pany_ldfiltered <- pany_ldfiltered[!duplicated(pany_ldfiltered$SNP),]

psingle_ldfiltered <- do.call(rbind, 
                              lapply(seq_along(1:nrow(psingle)), function(snp) 
                                  filterLD(snp=snp, ld=ld, p_original=psingle)))
psingle_ldfiltered <- psingle_ldfiltered[!duplicated(psingle_ldfiltered$SNP),]

pany_ldfiltered <- do.call(rbind, lapply(seq_along(1:nrow(pany)), function(snp) 
    filterLD(snp=snp, ld=ld, p_original=pany)))
pany_ldfiltered <- pany_ldfiltered[!duplicated(pany_ldfiltered$SNP),]

psingle_ldfiltered <- do.call(rbind, 
                              lapply(seq_along(1:nrow(psingle)), function(snp) 
                                  filterLD(snp=snp, ld=ld, p_original=psingle)))
psingle_ldfiltered <- psingle_ldfiltered[!duplicated(psingle_ldfiltered$SNP),]

# get unique SNPs from LD filtered, most significant
ldfilteredSNPs <- unique(c(pany_ldfiltered$SNP, psingle_ldfiltered$SNP))


# Combine pvalues for manhattan plot  
pany_commonfiltered <- pany[pany$SNP %in% ldfilteredSNPs,]
pany_commonfiltered$MARKER <- ""
psingle_commonfiltered <- psingle[psingle$SNP %in% ldfilteredSNPs,]
psingle_commonfiltered$MARKER <- ""
p_ldfiltered <- rbind(psingle_commonfiltered, pany_commonfiltered, p_pruned_sig)
p_ldfiltered$TYPE <- as.factor(p_ldfiltered$TYPE)
p_ldfiltered$MARKER <- as.factor(p_ldfiltered$MARKER)

# Get corresponding beta values and set betas of non-significant SNPs to 0
bany_ldfiltered <- bany[which(bany$SNP %in% ldfilteredSNPs),]
bany_ldfiltered[which(pany_commonfiltered$P >= fdr_single),2:42] <- 0
write.table(bany_ldfiltered, 
            paste(directory, "/GWAS/",  strftime(Sys.time(), "%Y%m%d"), 
                  "bany_ldfiltered_nonSigBetasEqZero.txt", sep=""), 
            sep="\t",col.names=TRUE, row.names=FALSE) 

# Cluster traits based on dotproduct (angle between beta-vectors)
bany_ldfiltered_dist <- dotproddist(t(as.matrix(bany_ldfiltered [,-c(1,43)])))
bany_ldfiltered_clustered <- hclust(bany_ldfiltered_dist)
bany_ldfiltered <- data.frame(SNP=bany_ldfiltered$SNP, 
                bany_ldfiltered[,-c(1,43)][,bany_ldfiltered_clustered$order])
bany_ldfiltered.m <- melt(bany_ldfiltered) 

# Set up plotting parameters
colManhattan <- wes_palette(20, name="Moonrise1", type='continuous')[c(16,4)]
colEffectSizes=colorRampPalette(c(
    wes_palette(16, name="Moonrise2", type='continuous')[1], "white", 
    wes_palette(16, name="Moonrise2", type='continuous')[6]))((5))

# manhattan plot
pman <- manhattan(p_ldfiltered, genomewideline=-log10(fdr_single), 
                  size.y.labels=10, size.x.labels=10, xscale=TRUE, mtvsst=TRUE, 
                  cols=colManhattan, a=0.8, 
                  colGenomewideline=wes_palette(20, name="Moonrise1", 
                                                type='continuous')[19]) + 
        theme(panel.border = element_blank(),
          legend.text=element_text(size=10),
          legend.title=element_text(size=12)) +
        labs(colour="GWAS set-up")

# get SNp position for effect size heatmap
bany_ldfiltered.m$pos <- sapply(bany_ldfiltered.m$SNP, function(s) 
    unique(pman$data$pos[pman$data$SNP==s]) )

# heatmap of effect sizes
phm <- ggplot(bany_ldfiltered.m, aes(pos, variable)) + 
        geom_point(aes(colour=value), shape='|', size=5) +
        scale_color_gradientn(colours=colEffectSizes, 
                              values=c(0, 0.2, 
                                       1-range(bany_ldfiltered.m$value)[2]/
                                        sum(abs(range(bany_ldfiltered.m$value)))
                                       , 0.8, 1)) +
        scale_y_discrete(position = "right", expand=c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        labs(colour="Effect size") +
        theme(axis.line = element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.length= unit(0, "lines"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10),
          legend.text=element_text(size=10),
          legend.title=element_text(size=12),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())

# plot via cowplot funcion plot_grid
pdf(paste(directory, "/GWAS/",  strftime(Sys.time(), "%Y%m%d"), 
          "ManhattanEffectsizes.pdf", sep=""), width=20, height=12)
plot_grid(pman, phm, ncol=1, nrow=2, align='v', labels=c('(a)', '(b)'))
dev.off()

png(paste(directory, "/GWAS/",  strftime(Sys.time(), "%Y%m%d"), 
          "ManhattanEffectsizes.png", sep=""), width=1200, height=1200)
plot_grid(pman, phm, ncol=1, nrow=2, align='v', labels=c('(a)', '(b)'))
dev.off()





