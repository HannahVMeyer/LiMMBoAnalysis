###########################################################
###                                                     ###
### Analyse results from GWAS of 41 yeast quantitative  ###
### growth traits via yeast/GWAS/GWAS_yeast.sh          ###
###                                                     ###
###     * dataset from Bloom et al 2013                 ###
###     * phenotypes processed via                      ###
###         yeast/phenotypes/phenotypes_yeast.R         ###
###     * genotypes processed via                       ###
###         yeast/genotypes/genotypes_yeast.R           ###
###     * relationship estimated via                    ###
###         yeast/genotypes/relationship_yeast.R        ###
###     * GWAS via  yeast/GWAS/GWAS_yeast.sh            ###
###         * univariate LMM                            ###
###         * multivariate LMMs with LiMMBo             ###
###                                                     ###
### Generates Figure 5 (publication)                    ###
###           Figure 5.5, 5.6 (thesis)                  ###
###                                                     ###
###########################################################



###############################
### Libraries and Functions ###
###############################
library("data.table")
library("GenomicRanges")
library("plyr")
library("pvclust")
library("splitstackshape")
library("ggplot2")
library("cowplot")
library('ggdendro')
library('dendextend')
library("wesanderson")
source("~/GWAS/analysis/general/pvclustHM.R")
source("~/LiMMBo/yeast/GWAS/manhattanplot.R")

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


chrArab2chrRom <- function(x, direction="arab2rom") {
        if (direction == "arab2rom") {
            switch(EXPR=x, "1"="I", "2"="II", "3"="III", "4"="IV", "5"="V", 
                   "6"="VI", "7"="VII", "8"="VIII", "9"="IX", "10"="X", 
                   "11"="XI", "12"="XII", "13"="XIII", "14"="XIV", "15"="XV", 
                   "16"="XVI", "Mito"="Mito")
        } else {
            switch(EXPR=x, "I"=1, "II"=2, "III"=3, "IV"=4, "V"=5, "VI"=6, 
                   "VII"=7, "VIII"=8, "IX"=9, "X"=10, "XI"=11, "XII"=12, 
                   "XIII"=13, "XIV"=14, "XV"=15, "XVI"=16, "Mito"=17, "M"=17)
        }
}

############
### data ###
############
directory='~/data/LiMMBo/feasabilityBootstrap/yeast'

# fdr thresholds for significcance in multi-variate uni-variate LMM via
# permutation obtained from GWAS_yeast.sh calling gwas.py 
fdr_multi=1.150566e-05
fdr_single=8.636346e-06

# tag SNPs within 3kb window with r2 > 0.8 according to Liti et al (2009) Nature
ld <- fread(paste(directory, "/inputdata/BYxRM.3kb.tags.list", sep=""),sep=" ", 
            data.table=FALSE)  
ld$ID <- gsub("\\d*_chr0?(\\d{1,2})_(\\d*)_.*", "chr\\1:\\2",ld$SNP)

# p values of multi-variate, any effect LMM 
pany <- fread(paste(directory, "/GWAS/pvalues_lmm_any.csv", sep=""), sep=",", 
              data.table=FALSE)
colnames(pany) <- c("ID", "P")
pany$SNP <- ld$SNP
pany$CHR <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\1", pany$ID))
pany$BP <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\2", pany$ID))
pany$TYPE <- 'multitrait'

# p values of single-variate LMM 
psingle <- fread(paste(directory, "/GWAS/pvalues_lmm_single_min_adjust.csv", 
                       sep=""), sep=",", data.table=FALSE, header=TRUE)
colnames(psingle) <- c("ID", "P")
psingle$SNP <- ld$SNP
psingle$CHR <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\1", psingle$ID))
psingle$BP <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\2", psingle$ID))
psingle$TYPE <- 'singletrait'

# effect sizes
bany <- fread(paste(directory, "/GWAS/betas_lmm_any.csv", sep=""), sep=",", 
data.table=FALSE)
colnames(bany)[1] <- "SNP"
bany$ID <- ld$SNP

# yeast chromosome sizes generated via 
# fetchChromSizes sacCer3 sacCer3.chrom.sizes > sacCer3.chrom.sizes.txt
chrom_size <- read.table(paste(directory, "/inputdata/sacCer3.chrom.sizes.txt",
                               sep=""), stringsAsFactors=FALSE)
chrom_size$chr <- unlist(sapply(gsub("chr", "", chrom_size$V1), 
                         chrArab2chrRom, direction="rom2arab"))
chrom_size <- chrom_size[order(chrom_size$chr),]
colnames(chrom_size)[1:2] <- c("CHR", "length")

################
### analysis ###
################

####################
### SNP analysis ###
####################

# Filter pvalues based on LD (retain SNP with lowest pvalue from ld block)
pany_ldfiltered <- do.call(rbind, lapply(seq_along(1:nrow(pany)), function(snp) 
    filterLD(snp=snp, ld=ld, p_original=pany)))
pany_ldfiltered <- pany_ldfiltered[!duplicated(pany_ldfiltered$SNP),]

psingle_ldfiltered <- do.call(rbind, 
                              lapply(seq_along(1:nrow(psingle)), function(snp) 
                                  filterLD(snp=snp, ld=ld, p_original=psingle)))
psingle_ldfiltered <- psingle_ldfiltered[!duplicated(psingle_ldfiltered$SNP),]

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
write.table(sig_tags, file=
                paste(directory, "/GWAS/significant_tagSNPs.txt",sep=""),
            col.names=NA, row.names=TRUE, quote=FALSE, sep="\t")


# chose 15kb for displaying marker SNPs
pruned <- fread(paste(directory, "/inputdata/BYxRM.", 15, "kb.pruned.bim", 
                      sep=""), data.table=FALSE)
psingle_pruned <- psingle[which(psingle$SNP %in% pruned[,2]),]
pany_pruned <- pany[which(pany$SNP %in% pruned[,2]),]
psingle_pruned_sig <- psingle_pruned[psingle_pruned$P < fdr_single,]
pany_pruned_sig <- pany_pruned[pany_pruned$P < fdr_single,]

p_pruned_sig <- rbind(psingle_pruned_sig, pany_pruned_sig)
p_pruned_sig$MARKER <- "tag"

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

# get yeast annotation and store as GRanges and gtf
yeastgenes <- fread(paste(directory, 
                          '/inputdata/20170207_ScerevisaeR64-1-1.txt', sep=""),
                    data.table=FALSE, stringsAsFactors = FALSE)
colnames(yeastgenes) <- c("geneID", "transcriptID", "chr","start", "end",
                          "strand", "geneName" )
yeastgenes$CHRROM <- unlist(sapply(yeastgenes$chr, chrArab2chrRom, 
                                   direction="rom2arab"))
yeastgenes <- yeastgenes[order(yeastgenes$CHRROM, yeastgenes$start),]
grYeast = with(yeastgenes, GRanges(chr, 
                                   IRanges(start=start, end=end, names=geneID), 
                                   strand=strand, geneID=geneID, 
                                   transcriptID=transcriptID,
                                   geneName=geneName))

# map snps to genes
snpsGWAS <- data.frame(chrArab= pany$CHR, start=pany$BP, end= pany$BP+1, 
                            SNP=pany$SNP, ID=pany$ID)
snpsGWAS$chr <- sapply(snpsGWAS$chrArab, chrArab2chrRom)
grSnpsGWAS = with(snpsGWAS, GRanges(chr, IRanges(start=start, end=end,  
                                                           names=ID), SNPID=ID))

snpsGWASGenes <- subsetByOverlaps(grYeast, grSnpsGWAS)
overlaps <- findOverlaps(grYeast,grSnpsGWAS)
rsid <- CharacterList(split(names(grSnpsGWAS)[subjectHits(overlaps)], 
                            queryHits(overlaps)))
mcols(snpsGWASGenes) <- DataFrame(mcols(snpsGWASGenes), rsid)
snpsGWASGenes.df <- data.frame(snpsGWASGenes, stringsAsFactors=FALSE)
snpsGWASGenes.df$rsid <- sapply(snpsGWASGenes.df$rsid, paste, collapse=",")
snpsGWASGenes.df$seqnames <- as.character(snpsGWASGenes.df$seqnames)
snpsGWASGenes.df$strand <- as.character(snpsGWASGenes.df$strand)

write.table(snpsGWASGenes.df, file=
                paste(directory, "/GWAS/ReferenceGWASGenes.txt",sep=""),
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# convert from genes with multiple SNPs, to single SNPs mapping to genes
ID2Genes <- apply(as.matrix(snpsGWASGenes.df), 1, function(x) {
    ids <- unlist(strsplit( x[9], ","))
    df <- expandRows(as.data.frame(t(x)), length(ids), count.is.col = FALSE)
    df$SNP <- ids
    return(df)
})
ID2Genes <- do.call(rbind, ID2Genes)
beta_map <- merge(ID2Genes, bany, by.x=10, by.y=1)
pbeta_map <- merge(beta_map[,c(1:2, 7, 9, 11:52)], pany[,1:2], by=1)

# some SNPs are mapped to more than one gene; get these SNPs + paste gene column
duplicateSNPs <- pbeta_map[which(pbeta_map$SNP %in% 
                                     pbeta_map$SNP[duplicated(pbeta_map$SNP)]),]
mergedSNPs <- lapply(unique(duplicateSNPs$SNP), function(s) {
    geneIDs <- paste(as.character(duplicateSNPs$geneID[duplicateSNPs$SNP == s]), 
                     collapse=",")
    geneNames <- paste(
        as.character(duplicateSNPs$geneName[duplicateSNPs$SNP == s]), 
        collapse=",")
    data.frame(SNP=s, seqnames=duplicateSNPs[duplicateSNPs$SNP == s, ][1,2],
               geneID=geneIDs, geneName=geneNames,
               duplicateSNPs[duplicateSNPs$SNP == s,][1, 5:ncol(duplicateSNPs)])
})
mergedSNPs <- do.call(rbind, mergedSNPs)
pbeta_map <- pbeta_map[ -which(pbeta_map$SNP %in% 
                       pbeta_map$SNP[duplicated(pbeta_map$SNP)]),]
pbeta_map <- rbind(pbeta_map, mergedSNPs)

############################
### effect size analysis ###
############################

# Remove all non-significant SNPs
pbeta_map_ld <- pbeta_map[pbeta_map$SNP %in% pany_ldfiltered$ID,]
pbeta_map_sig <- pbeta_map_ld[pbeta_map_ld$P < fdr_single,]
rownames(pbeta_map_sig) <- pbeta_map_sig$SNP
saveRDS(pbeta_map_sig, paste(directory, "/GWAS/pbeta_map_ld_sig.rds", sep=""))

# cluster and determine stability of clusters with pvclust
method.dist <- "correlation"

## for traits
pbeta_map_traits_pv <- pvclust(pbeta_map_sig[,5:45], nboot=50000, iseed=10,
                          method.dist=method.dist)
saveRDS(pbeta_map_traits_pv, paste(directory, "/GWAS/",  
                              strftime(Sys.time(), "%Y%m%d"), 
                              "pbeta_map_traits_pv_50k_fdr5e-5", 
                              method.dist, ".rds", sep=""))
# pbeta_map_traits_pv <- readRDS( paste(directory, "/GWAS/20170626pbeta_map_traits_pv_50k_fdr5e-5correlation.rds", sep=""))
pvclustPlot(pbeta_map_traits_pv, print.bp=FALSE, print.edge=TRUE,
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(pbeta_map_traits_pv, max.only=TRUE, highestEdge = 20)
pbeta_map_traits_clusters <- pvclustPick(pbeta_map_traits_pv, 
                                         max.only=TRUE
                                         , highestEdge = 20)$clusters
pbeta_map_traits_clusters <- lapply(seq_along(pbeta_map_traits_clusters),
                                    function(x) {
                                        pbeta_map_traits_clusters[[x]] <- 
                                        cbind(pbeta_map_traits_clusters[[x]], 
                                              paste("cluster", x))
                                    })
pbeta_map_traits_clusters <- data.frame(do.call(rbind, 
                                                pbeta_map_traits_clusters),
                                        stringsAsFactors=FALSE)
colnames(pbeta_map_traits_clusters) <- c("ID", "cluster")

## for SNPs
pbeta_map_snps_pv <- pvclust(t(pbeta_map_sig[,5:44]), nboot=10000, iseed=10, 
                               method.dist=method.dist)
saveRDS(pbeta_map_snps_pv, paste(directory, "/GWAS/",  
                                   strftime(Sys.time(), "%Y%m%d"), 
                                   "pbeta_map_snps_pv_", 
                                   method.dist, ".rds", sep=""))
#pbeta_map_snps_pv <- readRDS("~/GWAS/data/LiMMBo/feasabilityBootstrap/yeast/GWAS/20170526pbeta_map_snps_pv_fdr5e-5correlation.rds")

pvclustPlot(pbeta_map_snps_pv, print.bp=FALSE, print.edge=TRUE, labels=FALSE, 
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(pbeta_map_snps_pv, highestEdge = 180)
pbeta_map_snps_clusters <- pvclustPick(pbeta_map_snps_pv, 
                                       highestEdge = 180)$clusters
pbeta_map_snps_clusters <- lapply(seq_along(pbeta_map_snps_clusters),
                                    function(x) {
                                        pbeta_map_snps_clusters[[x]] <- 
                                            cbind(pbeta_map_snps_clusters[[x]], 
                                                  paste("cluster", x))
                                    })
pbeta_map_snps_clusters <- data.frame(do.call(rbind, 
                                                pbeta_map_snps_clusters),
                                        stringsAsFactors=FALSE)
colnames(pbeta_map_snps_clusters) <- c("ID", "cluster")

## order according to cluster order
snp <- pbeta_map_sig$SNP
snp <- factor(snp[pbeta_map_snps_pv$hclust$order], 
              levels=unique(snp[pbeta_map_snps_pv$hclust$order])) 

pbeta_map_sig_order <- pbeta_map_sig[,5:45][pbeta_map_snps_pv$hclust$order,
                                    pbeta_map_traits_pv$hclust$order]
pbeta_map_sig_order.m <- melt(data.frame(SNP=snp, pbeta_map_sig_order), 
                              value.name = "beta",
                              variable.name="trait")

chr <- gsub("chr(\\d{1,2}):.*", "\\1", rownames(pbeta_map_sig_order))
chr <- factor(chr, levels=1:16)

#############
### plots ###
#############


##  rectangles for significant SNP clusters across different chromosomes
cluster_cols <-  data.frame(ID=rownames(pbeta_map_sig_order), 
                            stringsAsFactors = FALSE)
cluster_cols$chr <- gsub("chr(\\d{1,2}).*", "\\1", cluster_cols$ID)
cluster_cols$cluster <- sapply(cluster_cols$ID, function(x) {
    if(any(pbeta_map_snps_clusters$ID == x)) {
        pbeta_map_snps_clusters$cluster[pbeta_map_snps_clusters$ID == x]
    } else {
        "no_cluster"
    }})
cluster_cols$geneID <- pbeta_map_sig$geneID[pbeta_map_snps_pv$hclust$order]
cluster_cols$geneName <- pbeta_map_sig$geneName[pbeta_map_snps_pv$hclust$order]

cluster2chromosomes <- sapply(unique(cluster_cols$cluster) , function(x) {
    tmp <- cluster_cols[which(cluster_cols$cluster == x),]
    if (length(unique(tmp$chr)) == 1) {
        return("single chromosome")
    } else {
        return(paste(length(unique(tmp$chr)), "chromosomes"))
    }
})
multipleChr <- cluster2chromosomes[! grepl("single", cluster2chromosomes)]
multipleChr <- multipleChr[!grepl("no", names(multipleChr))]

cluster_cols$clustertype <- sapply(cluster_cols$cluster, function(x) {
    cluster2chromosomes[which(names(cluster2chromosomes) == x)]
})

write.table(cluster_cols, paste(directory, "/GWAS/", strftime(Sys.time(), 
                                                              "%Y%m%d"), 
                  "SNP_clusters.csv", sep=""), sep=",", quote=FALSE, 
            col.names=TRUE)

cluster_cols_start <- which(!duplicated(cluster_cols$cluster))
cluster_cols_start <- cluster_cols_start[-which(
    grepl("no",unique(cluster_cols$cluster)))]
cluster_cols_end <- sapply(cluster_cols_start, function(x) 
    max(which(cluster_cols$cluster==cluster_cols$cluster[x])))
cluster_cols_rect <- data.frame(xmin=cluster_cols_start, 
                                xmax=cluster_cols_end, 
                                ymin=rep(0, length(cluster_cols_start)), 
                                ymax=rep(42, length(cluster_cols_start)), 
                                colour=unique(cluster_cols$cluster)[-which(
                                    grepl("no",unique(cluster_cols$cluster)))])

cluster_cols_rect_multiple <- cluster_cols_rect[which(cluster_cols_rect$colour
                                                      %in% names(multipleChr)),]
cluster_cols_rect_multiple$colour <- factor(cluster_cols_rect_multiple$colour, 
                                            levels= as.character(
                                            cluster_cols_rect_multiple$colour))

##  rectangles for significant trait clusters 
cluster_rows <-  data.frame(ID=colnames(pbeta_map_sig_order))
cluster_rows$cluster <- sapply(cluster_rows$ID, function(x) {
    if(any(pbeta_map_traits_clusters$ID == x)) {
        pbeta_map_traits_clusters$cluster[pbeta_map_traits_clusters$ID == x]
    } else {
        "no_cluster"
    }})

cluster_rows_start <- which(!duplicated(cluster_rows$cluster))
cluster_rows_start <- cluster_rows_start[-which(
    grepl("no",unique(cluster_rows$cluster)))]
cluster_rows_end <- sapply(cluster_rows_start, function(x) 
    max(which(cluster_rows$cluster==cluster_rows$cluster[x])))
cluster_rows_rect <- data.frame(ymin=cluster_rows_start, 
                                ymax=cluster_rows_end, 
                                xmin=rep(0, length(cluster_rows_start)), 
                                xmax=rep(dim(pbeta_map_sig_order)[1],
                                         length(cluster_rows_start)), 
                                colour=unique(cluster_rows$cluster)[-which(
                                    grepl("no",unique(cluster_rows$cluster)))])



### ggplots for Figure 5.5, 5.6 (thesis) and Figure 5 (publication)

##   plotting parameters
colManhattan <- wes_palette(5, name="Darjeeling", type='continuous')[c(5,4)]
colStrength <- c(wes_palette(5, name="Darjeeling", type='continuous')[c(2)],
                 "white",
                 wes_palette(5, name="Darjeeling", type='continuous')[4])
Colcol <- wes_palette(18, name="Darjeeling", type='continuous')[-c(1:2)]
linesize <- 0.5
textsize <- 10
legendsize <- 10

# manhattan plot
pman <- manhattan(p_ldfiltered, genomewideline=-log10(fdr_single), 
                  size.y.labels=textsize, size.x.labels=textsize, xscale=TRUE, 
                  mtvsst=TRUE, cols=colManhattan, a=0.8, 
                  colorGenomewideline=colManhattan[2]) + 
    theme(panel.border = element_blank(),
          legend.text=element_text(size=textsize),
          legend.title=element_text(size=textsize),
          legend.key.size = unit(0.5,"line")) +
    labs(colour="GWAS set-up")

# chromosome location of clustered SNPs
chr_snp <-  ggplot(data.frame(SNP=snp, variable=1, chr=chr, 
                              stringsAsFactors=FALSE)) + 
    geom_tile(aes(x=snp, y=variable, fill=chr)) +
    labs(fill="Chromosome" ) +
    scale_fill_manual( values = Colcol, guide=guide_legend(nrow=2)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position ='bottom',
          axis.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.text=element_text(size=legendsize, angle = 0, hjust = 1),
          legend.title=element_text(size=legendsize),
          legend.key.size = unit(0.5,"line"),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
legendChr <- get_legend(chr_snp)
legendChr$layout$t[1] <- 1 


# dendrogram of clustered SNPs
relabeled_dendro_snps <- pbeta_map_snps_pv$hclust %>% 
as.dendrogram %>% pvclust_show_signif(pbeta_map_snps_pv, signif_type='au')
snp_dendro <- as.ggdend(relabeled_dendro_snps )

dendrogram_snp <- ggplot() + geom_blank() +
    geom_segment(data = segment(snp_dendro), 
                 aes(x = x, y = y, xend = xend, yend = yend,
                 colour= as.factor(lwd))) +
    scale_color_manual(values=c("black",'#03569b'),
                guide=FALSE) +
    labs(title="SNPs") +
    scale_x_continuous(expand=c(0,0)) +
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          plot.title =  element_text(size = textsize),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

# dendrogram of trait clustering
relabeled_dendro_traits <- pbeta_map_traits_pv$hclust %>% 
    as.dendrogram %>% pvclust_show_signif(pbeta_map_traits_pv, signif_type='au')
traits_dendro <- as.ggdend(relabeled_dendro_traits )

dendrogram_traits <- ggplot() + geom_blank() +
    geom_segment(data = segment(traits_dendro), 
                 aes(x = x, y = y, xend = xend, yend = yend,
                 colour= as.factor(lwd))) +
    scale_x_continuous(breaks = seq_along(traits_dendro$labels$label), 
                       labels = traits_dendro$labels$label, 
                       expand=c(0,0), 
                       limits=c(0,42)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values=c("black",
                                '#03569b'),
                       guide=FALSE) +
    coord_flip() + 
    theme(axis.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = textsize, angle = 0, hjust = 1))

labelTraitsDendro <- ggplot() + 
    annotate("text", label = "Traits", x = 1, y = 1, size = 4, 
             fontface="bold") +
    theme_nothing()

# label cluster 28 and 31
specificCluster <- cluster_cols_rect_multiple[
                    cluster_cols_rect_multiple$colour %in% 
                        c("cluster 28", "cluster 31"),]
specificCluster$xmean <- (specificCluster$xmax + specificCluster$xmin)/2
specificCluster$label <- c("a", "b")
labelClusters <- ggplot(data=specificCluster,aes(x=xmean, y=ymin, 
                                                 label=label)) + 
                geom_text(size=4) +
                scale_x_continuous(expand=c(0,0), limits=c(0,210)) +
                scale_y_continuous(expand=c(0,0)) +
                theme(axis.title=element_blank(),
                      axis.line = element_blank(),
                      axis.ticks=element_blank(),
                      axis.text=element_blank(),          
                      plot.margin = unit(c(-8, 0, -8, 0), "cm"))

# effect sizes clustered by SNPs and traits and 
colStrength <- c('#ca0020','#f4a582','#ffffff','#bababa','#404040')
beta_snp_hm <-  ggplot(pbeta_map_sig_order.m) + 
    geom_tile(aes(x=as.numeric(SNP), y=as.numeric(trait), 
                  fill=beta), height=1) +
    scale_fill_gradientn( colours = colorRampPalette(colStrength)(100), 
                          name = "Effect size") +
    geom_rect(data=cluster_cols_rect_multiple, aes(xmin=xmin + 0.5, 
                                                   xmax=xmax + 0.5, 
                                          ymin=ymin, 
                                          ymax=ymax), color= 'grey', 
              size=linesize, fill = NA) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0,42),
                       breaks = seq_along(traits_dendro$labels$label)) +
    theme(legend.position ='left',
          axis.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.text=element_text(size=legendsize, angle = 0, vjust = 0),
          legend.title=element_text(size=legendsize, hjust = 1),
         # legend.key.size = unit(0.5,"line"),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
legendAbs <- get_legend(beta_snp_hm )

# combined effect size plot
effect_sizes <- plot_grid(NULL, dendrogram_snp, NULL,
                          NULL , chr_snp + theme(legend.position ='none'), 
                            labelTraitsDendro, 
                          legendAbs, beta_snp_hm + theme(legend.position =
                                                             'none'), 
                          dendrogram_traits,
                          NULL , labelClusters,  NULL,
                          NULL , legendChr,  NULL, 
                          ncol=3, nrow=5,  rel_heights=c(2,0.5, 7, 0.5, 2), 
                          rel_widths=c(0.5, 4,2), 
                          align='vh')


## combine all plots
png(paste(directory, "/GWAS/manhattan_effectsizes.png", sep=""), height=12, 
    width=10, res=450, unit="in" )
plot_grid(pman, effect_sizes, nrow=2,  rel_heights=c(3,5), labels=c("", ""))
dev.off()

pman_sep <- manhattan(p_ldfiltered, genomewideline=-log10(fdr_single), 
                      size.y.labels=12, size.x.labels=12, 
                      mtvsst=TRUE, color=colManhattan, a=0.8,
                      colorGenomewide =  "gray50", chromboundaries = TRUE) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size=12),
          axis.title =element_text(size=12)) +
    labs(colour="GWAS set-up")

ggsave(plot=pman_sep, height=4, width=10,
       file=paste(directory,"/GWAS/manhattanplot.png", sep=""))

png(paste(directory, "/GWAS/effectsizes.png", sep=""), height=10, 
    width=10, res=450, unit="in" )
effect_sizes
dev.off()
