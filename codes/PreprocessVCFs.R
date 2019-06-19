###
#   File name : PreprocessVCFs.R
#   Author    : Hyunjin Kim
#   Date      : Jun 13, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load VCF files, combine them, and annotate to make a RDA file
#
#   Instruction
#               1. Source("PreprocessVCFs.R")
#               2. Run the function "preprocessVCFs" - specify the necessary input file paths and output RDA file path
#               3. The RDA file will be generated as the output path
#
#   Example
#               > source("The_directory_of_PreprocessVCFs.R/PreprocessVCFs.R")
#               > preprocessVCFs(vcfFilePath="./data/somatic_mutation/",
#                                sampleSheetPath="./data/somatic_mutation/gdc_sample_sheet_tcga_coad_read_mutect2_vcfs.tsv",
#                                outputPath="./data/somatic_mutation/somatic_mutation_vcfs_tcga_coad_read.rda")
###

### THE FILE SIZE IS TOO LARGE, THEREFORE, I DID NOT UPLOAD THE DATA ON THE SERVER,
### AND JUST RAN IT ON MY LOCAL COMPUTER
preprocessVCFs <- function(vcfFilePath="./data/somatic_mutation/",
                           sampleSheetPath="./data/somatic_mutation/gdc_sample_sheet_tcga_coad_read_mutect2_vcfs.tsv",
                           outputPath="./data/somatic_mutation/somatic_mutation_vcfs_tcga_coad_read.rda") {
  
  ### load library
  if(!require(vcfR)) {
    install.packages("vcfR")
    library(vcfR)
  }
  
  ### load vcf file paths
  f <- list.files(vcfFilePath, pattern = "*.vcf$")
  
  ### load vcf files
  vcf <- vector("list", length = length(f))
  names(vcf) <- f
  for(file in f) {
    temp <- read.vcfR(paste0(vcfFilePath, file), verbose = FALSE)
    vcf[[file]] <- data.frame(temp@fix, temp@gt, stringsAsFactors = FALSE, check.names = FALSE)
  }
  gc()
  
  ### load vcf sample info
  vcf_sample_info <- read.table(file = sampleSheetPath, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE, check.names = FALSE)
  
  ### refine file Name and Case ID column in the sample info
  vcf_sample_info$`File Name` <- sapply(vcf_sample_info$`File Name`, function(x) {
    substr(x, 1, nchar(x)-3)
  })
  vcf_sample_info$`Case ID` <- sapply(vcf_sample_info$`Case ID`, function(x) {
    strsplit(x, split = ", ", fixed = TRUE)[[1]][1]
  })
  
  ### remove duplicated patients
  rownames(vcf_sample_info) <- vcf_sample_info$`File Name`
  vcf_sample_info <- vcf_sample_info[names(vcf),]
  dup_idx <- duplicated(vcf_sample_info$`Case ID`)
  vcf_sample_info <- vcf_sample_info[!dup_idx,]
  vcf <- vcf[!dup_idx]
  
  ### change the names of the objects to Case ID
  rownames(vcf_sample_info) <- vcf_sample_info$`Case ID`
  names(vcf) <- rownames(vcf_sample_info)
  
  ### write out as RDA file
  save(list = c("vcf", "vcf_sample_info"), file = outputPath)
  
}
