## parameters and input data ####
option_list <- list(
    optparse::make_option(c("-k", "--kinship"), action="store", dest="kinship",
               type="character", help="Path to kinship file [default:
               %default].", default=NULL)
)

args <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

kinship <- data.table::fread(args$kinship, data.table=FALSE,
                             stringsAsFactors=FALSE)

name <- gsub("(.*/.*kinship).*\\.csv", "\\1", args$kinship)

## eigendecomposition ####
ed <- eigen(kinship, symmetric=TRUE)

## save eigenvalues and vectors of kinship matrix ####
# LIMIX/LiMMBo format
write.table(ed$vectors, paste(name, "_eigenvec_limmbo.csv", sep=""),
    sep=",", quote=FALSE, col.names=names(kinship), row.names=FALSE)
write.table(ed$values, paste(name, "_eigenvalue_limmbo.csv", sep=""),
    sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)

# sbat/GEMMA format
write.table(ed$vectors, paste(name, "_eigenvec.txt", sep=""),
    sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(ed$values, paste(name, "_eigenvalue.txt", sep=""),
    sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
