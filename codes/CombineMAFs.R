###
#   File name : CombineMAFs.R
#   Author    : Hyunjin Kim
#   Date      : Jun 24, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load MAF files, combine them, and preprocess to make one MAF file
#
#   Instruction
#               1. Source("CombineMAFs.R")
#               2. Run the function "combineMAFs" - specify the necessary input file paths and output MAF file path
#               3. The MAF file will be generated as the output path
#
#   Example
#               > source("The_directory_of_CombineMAFs.R/CombineMAFs.R")
#               > combineMAFs(MAF_path="./data/somatic_mutation/MAF/",
#                             outputPath="./data/somatic_mutation/MAF/somatic_mutation_maf_tcga_coad_read.maf")
###

combineMAFs <- function(MAF_path="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/somatic_mutation/",
                        outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/somatic_mutation/somatic_mutation_maf_tcga_coad_read.maf") {
  
  ### load libraries
  if(!require(maftools, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("maftools", version = "3.8")
    require(maftools, quietly = TRUE)
  }
  
  ### load MAF file paths
  f <- list.files(path = MAF_path, pattern = "*.maf$")
  
  ### load MAFs
  for(i in 1:length(f)) {
    if(i == 1) {
      maf <- read.table(file = paste0(MAF_path, f[i]), header = TRUE, sep = "\t", quote = "",
                        stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      maf <- rbind(maf, read.table(file = paste0(MAF_path, f[i]), header = TRUE, sep = "\t", quote = "",
                                   stringsAsFactors = FALSE, check.names = FALSE))
    }
  }

  ### remove duplicates in the combined MAF file
  duplicates <- duplicated(maf[,c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome",
                                  "Start_Position", "End_Position", "Strand", "Reference_Allele",
                                  "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode")])
  maf <- maf[!duplicates,]
  
  ### write out the MAF file
  write.table(maf, file = outputPath, sep = "\t", row.names = FALSE, quote = FALSE)
  
  ### these are optional comments at the top of the MAF file
  ### you could add the following to above of the file
  writeLines("#version gdc-1.0.0")
  writeLines("#filedate 20190624")
  writeLines("#annotation.spec gdc-1.0.1-public")
  writeLines(paste("#n.analyzed.samples", length(unique(maf$Tumor_Sample_Barcode))))
  writeLines(paste0("#tumor.aliquots.submitter_id ", paste(unique(maf$Tumor_Sample_Barcode), collapse = ",")))
  
}

