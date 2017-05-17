###############################
### Libraries and Functions ###
###############################
library("data.table")
library("GenomicRanges")
library("plyr")
library("pvclust")
source("~/GWAS/analysis/general/pvclustHM.R")

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
                   "XIII"=13, "XIV"=14, "XV"=15, "XVI"=16, "Mito"="Mito")
        }
}

############
### data ###
############

directory='~/GWAS/data/LiMMBo/feasabilityBootstrap/yeast'

fdr_multi=1.99537104453e-05
fdr_single=1.30007477173e-05

ld <- fread(paste(directory, "/inputdata/BYxRM.3kb.tags.list", sep=""),sep=" ", 
            data.table=FALSE)  
ld$ID <- gsub("\\d*_chr0?(\\d{1,2})_(\\d*)_.*", "chr\\1:\\2",ld$SNP)

pany <- fread(paste(directory, "/GWAS/pvalues_lmm_any.csv", sep=""), sep=",", 
              data.table=FALSE)
colnames(pany) <- c("ID", "P")
pany$SNP <- ld$SNP
pany$CHR <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\1", pany$ID))
pany$BP <- as.numeric(gsub("chr(\\d{1,2}):(\\d*)", "\\2", pany$ID))
pany$TYPE <- 'multitrait'

pany_ldfiltered <- do.call(rbind, lapply(seq_along(1:nrow(pany)), function(snp) 
    filterLD(snp=snp, ld=ld, p_original=pany)))
pany_ldfiltered <- pany_ldfiltered[!duplicated(pany_ldfiltered$SNP),]


bany <- fread(paste(directory, "/GWAS/betas_lmm_any.csv", sep=""), sep=",", 
              data.table=FALSE)
colnames(bany)[1] <- "SNP"
bany$ID <- ld$SNP

################
### analysis ###
################

## get yeast annotation and store as GRanges and gtf
### make granges
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

### map GWAS snps to genes
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
                paste(directory, "/GOenrichment/ReferenceGWASGenes.txt",sep=""),
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

ID2Genes <- apply(as.matrix(snpsGWASGenes.df), 1, function(x) {
    ids <- unlist(strsplit( x[9], ","))
    df <- expandRows(as.data.frame(t(x)), length(ids), count.is.col = FALSE)
    df$SNP <- ids
    return(df)
})
ID2Genes <- do.call(rbind, ID2Genes)
beta_map <- merge(ID2Genes, bany, by.x=10, by.y=1)
pbeta_map <- merge(beta_map, pany[,1:2], by=1)



# effect size analysis
## Remove all non-significant SNPs 
pbeta_map_sig <- pbeta_map[pbeta_map$P < fdr_single,]
bany_ldfiltered <- bany[which(bany[,1] %in% pany_ldfiltered[,1]),]
bany_sig <- bany_ldfiltered[pany_ldfiltered$P < fdr_single,]

## Magnitude of effect size (absolute values) 
bany_sig_abs <- abs(bany_sig[,-c(1,43)])
rownames(bany_sig_abs) <-  as.character(bany_sig[,1])

## pvclust
## cluster and determine stability of clusters with pvclust
#method.dist <- "euclidean"
method.dist <- "correlation"

### for traits
pbeta_map_traits_pv <- pvclust(pbeta_map_sig[,11:51], nboot=50000, iseed=10, 
                          method.dist=method.dist)
saveRDS(pbeta_map_traits_pv, paste(directory, "/GWAS/",  
                              strftime(Sys.time(), "%Y%m%d"), 
                              "pbeta_map_traits_pv_", 
                              method.dist, ".rds", sep=""))

### for SNPs 
pbeta_map_snps_pv <- pvclust(t(pbeta_map_sig[,11:51]), nboot=10000, iseed=10, 
                               method.dist=method.dist)
saveRDS(pbeta_map_traits_pv, paste(directory, "/GWAS/",  
                                   strftime(Sys.time(), "%Y%m%d"), 
                                   "pbeta_map_traits_pv_", 
                                   method.dist, ".rds", sep=""))


bany_snps_pv <- pvclust(t(bany_sig_abs), nboot=10000, iseed=10, 
                        method.dist=method.dist)
#tmp <- tt[match(bany_snps_pv[[1]]$labels,tt[,1]),]
#bany_snps_pv[[1]]$labels <- tmp[,2]
saveRDS(bany_snps_pv, paste(directory, "/GWAS/", strftime(Sys.time(), "%Y%m%d"), 
                            "bany_snps_pv_", method.dist, ".rds", sep=""))

pvclustPlot(bany_snps_pv, print.bp=FALSE, print.edge=TRUE, labels=FALSE, 
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(bany_snps_pv, highestEdge = 395)
bany_snps_clusters <- pvclustPick(bany_snps_pv, highestEdge = 395)

pvclustPlot(bany_traits_pv, print.bp=FALSE, print.edge=TRUE,
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(bany_traits_pv, max.only=FALSE)
pvclustPlot(pbeta_map_traits_pv, print.bp=FALSE, print.edge=TRUE,
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(pbeta_map_traits_pv, max.only=FALSE)
bany_traits_clusters <- pvclustPick(bany_traits_pv,  max.only=FALSE)

# get absolute value of effect sizes and order according to cluster order
bany_sig_abs <- t(apply(bany_sig[,-c(1,43)], 1, function(s) abs(s)))
bany_sig_abs_order <- bany_sig_abs_snp[bany_snps_pv$hclust$order,
                                       bany_traits_pv$hclust$order]
bany_sig_abs_order.m <- melt(data.frame(SNP=snp, BETA=bany_sig_abs_order))


## Ewan's script
beta = read.table("~/GWAS/data/LiMMBo/feasabilityBootstrap/yeast/beta_cleaned.txt",header = TRUE)
snp2gene = read.table("~/GWAS/data/LiMMBo/feasabilityBootstrap/yeast/snp2gene.txt",header = TRUE)

beta_gene = merge(beta,snp2gene, by.x = "snp", by.y = "snp")
beta_gene$count_01 = apply(beta_gene[,5:41],1,function(x){ sum(abs(x) > 0.01)})
                  
innermat = as.matrix(beta_gene[beta_gene$count_01 > 0,5:41])
hc = hclust(as.dist(1 - cor(t(innermat))), method = "average")
cc = hclust(as.dist(1 - cor((innermat))), method = "average")
heatmap(innermat, Colv = as.dendrogram(cc), Rowv = as.dendrogram(hc), 
        scale = "none", 
        labRow = beta_gene[beta_gene$count_01 > 0,]$sym, 
        cexRow = 0.15)
innermat_traits_pv <- pvclust(innermat, nboot=10000, iseed=10, 
                               method.dist=method.dist)
pvclustPlot(innermat_traits_pv, print.bp=FALSE, print.edge=TRUE,
            col.pv=wes_palette(4, name="Moonrise2", type='continuous')[2])
pvclustRect(innermat_traits_pv, max.only=FALSE)
## SNP with notable effect sizes across several traits
snps_across_traits <- sapply(c(0.5, 1, 1.5, 2), function(threshold) {
    apply(bany_sig_abs, 1, function(snp, thr) {
        length(which(snp >= thr))
    }, thr=threshold)
})

snps_across_traits_multiple <- apply(snps_across_traits, 2, function(thr) {
    which(thr > 3)
})

snps2genes <- lapply(snps_across_traits_multiple, function(thr) {
                   tmp <- do.call(c, lapply(names(thr), function(snp) {
                        snpsGWASGenes.df$geneID[grepl(snp, 
                                                      snpsGWASGenes.df$rsid)]
                    }))
                   unique(tmp)
})


