###
#   File name : UpdateSampleInfo.R
#   Author    : Hyunjin Kim
#   Date      : Jun 25, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : With the MAF file, update the sample info file with POLE mutation info
#
#   Instruction
#               1. Source("UpdateSampleInfo.R")
#               2. Run the function "update_sample_info" - specify the necessary input file paths and output sample info file path
#               3. The updated sample info file will be generated as the output path
#
#   Example
#               > source("The_directory_of_UpdateSampleInfo.R/UpdateSampleInfo.R")
#               > update_sample_info(mafFilePath="./data/somatic_mutation/MAF/somatic_mutation_maf_tcga_coad_read.maf",
#                                    sampleInfoPath="./data/coadread_tcga_clinical_data_updated.txt",
#                                    outputPath="./data/coadread_tcga_clinical_data_updated2.txt")
###

update_sample_info <- function(mafFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/somatic_mutation/somatic_mutation_maf_tcga_coad_read.maf",
                               sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated.txt",
                               outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt") {
  
  ### load library
  options(java.parameters = "-Xmx8000m")
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  ### load data
  maf <- read.table(file = mafFilePath, header = TRUE, sep = "\t", quote = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
  sample_info <- read.table(file = sampleInfoPath, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### find POLE-mutated (P286R, S459F, or V411L) samples
  pole_P286R_idx <- intersect(which(maf$Hugo_Symbol == "POLE"), which(maf$HGVSp_Short == "p.P286R"))
  pole_S459F_idx <- intersect(which(maf$Hugo_Symbol == "POLE"), which(maf$HGVSp_Short == "p.S459F"))
  pole_V411L_idx <- intersect(which(maf$Hugo_Symbol == "POLE"), which(maf$HGVSp_Short == "p.V411L"))
  pole_mutant_idx <- c(pole_P286R_idx, pole_S459F_idx, pole_V411L_idx)
  pole_mutant_samples <- substr(maf$Tumor_Sample_Barcode[pole_mutant_idx], 1, 12)
  
  ### update the sample info with the pole-mutant info
  sample_info$POLE_MUTANT <- FALSE
  sample_info$POLE_MUTANT[which(sample_info$`Patient ID` %in% pole_mutant_samples)] <- TRUE
  unique_samples_in_maf <- unique(substr(maf$Tumor_Sample_Barcode, 1, 12))
  sample_info$POLE_MUTANT[which(sample_info$`Patient ID` %in% setdiff(sample_info$`Patient ID`, unique_samples_in_maf))] <- NA
  
  ### write out the updated result
  write.table(sample_info, file = outputPath, sep = "\t", row.names = FALSE)
  write.xlsx2(sample_info, file = paste0(outputPath, ".xlsx"),
              sheetName = "tcga_coadread_sample_info", row.names = FALSE)
  
}
