## parameters and input data ####
option_list <- list(
    optparse::make_option(c("-p", "--pheno"), action="store", dest="pheno",
               type="character", help="Path to fd phenotype file [default:
               %default].", default=NULL),
    optparse::make_option(c("-c", "--covs"), action="store", dest="covs",
               type="character", help="Path to covariate file
               [default: %default].", default=NULL)
)

args <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

pheno <- data.table::fread(args$pheno, data.table=FALSE, stringsAsFactors=FALSE)
rownames(pheno) <- pheno[,1]
pheno <- pheno[,-1]

covs <- data.table::fread(args$covs, data.table=FALSE, stringsAsFactors=FALSE)
rownames(covs) <- covs[,1]
covs <- covs[,-1]

phenoname <- gsub("(.*/Ysim.*)\\.csv", "\\1", args$pheno)

## regress covariates from phenotypes ####
lm_res <- sapply(pheno, function(x) {
    tmp <- lm(y ~ ., data=data.frame(y=x, covs))
    tmp$residuals
})

## save residuals of regression ####
# LIMIX/LiMMBo format
write.table(lm_res, paste(phenoname, "_reg_limmbo.csv", sep=""),
    sep=",", quote=FALSE, col.names=NA, row.names=TRUE)

# sbat format
write.table(lm_res, paste(phenoname, "_reg_sbat.txt", sep=""),
    sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

# GEMMA format
write.table(lm_res, paste(phenoname, "_reg_gemma.txt", sep=""),
    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
