## Libraries ####
library("data.table")

## directory ####
directory <- "~/data/LiMMBo/rat/rawdata"

## read data ####
measures <- fread(paste(directory, "/arrayexpress/measures.txt", sep=""), 
                  data.table=FALSE, stringsAsFactors = FALSE)

pheno_ScientificData <- fread(paste(directory, 
                                 "/PublicationResults/ScientificDataS1.csv", 
                                 sep=""),
                           data.table=FALSE, stringsAsFactors = FALSE, 
                           skip=2)[, -c(6:7)]
colnames(pheno_ScientificData) <- gsub(" ", "", colnames(pheno_ScientificData))
pheno_ScientificData$gsub <- gsub("_normalized_by_batch", "",
                                  pheno_ScientificData$Measure)

pheno_NatGenet <- fread(paste(directory, 
                              "/PublicationResults/", 
                              "SupplementaryTable1_Phenotypes.csv", sep=""),
                              data.table=FALSE, stringsAsFactors = FALSE, 
                        skip=3)
colnames(pheno_NatGenet) <- gsub(" ", "", colnames(pheno_NatGenet))
pheno_NatGenet$MeasureinGSCAN <- gsub("_bc", "", pheno_NatGenet$MeasureinGSCAN)

## sorting data ####
covariate_names <- dplyr::filter(pheno_ScientificData, 
                                 Phenotypingtest == "Covariate")
phenotype_names <- dplyr::filter(pheno_ScientificData, 
                                 Phenotypingtest != "Covariate")

phenotype_names <- merge(phenotype_names, pheno_NatGenet[, -c(1:2,8)], 
                         by.x="gsub", by.y="MeasureinGSCAN", all=TRUE)
 
phenotype_names_notN <- dplyr::filter(phenotype_names, 
                                           Mappingmethod != "Mixed models",
                                           !is.na(Measure))
phenotype_names_N <- dplyr::filter(phenotype_names, 
                                  Mappingmethod == "Mixed models",
                                  !is.na(Measure))

phenotypes_normal <- measures[, colnames(measures) %in% 
                                  phenotype_names_N$Measure ]
rownames(phenotypes_normal) <- measures$SUBJECT.NAME

phenotypes_notNormal <- measures[, colnames(measures) %in% 
                                     phenotype_names_notN$Measure ]
rownames(phenotypes_notNormal) <- measures$SUBJECT.NAME

write.table(phenotypes_normal, paste(directory, "/phenotypes_normal.csv", 
                                     sep=""),
          row.names=TRUE, col.names=NA, quote=FALSE, sep=",")
write.table(phenotypes_notNormal, paste(directory, "/phenotypes_notNormal.csv", 
                                        sep=""),
            row.names=TRUE, col.names=NA, quote=FALSE, sep=",")








