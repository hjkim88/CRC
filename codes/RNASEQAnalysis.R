###
#   File name : RNASEQAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Apr 24, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DE analysis on the TCGA COAD/READ raw counts between different age groups
#               and between different race groups. Also perform pathway analysis on the DE genes
#
#   Instruction
#               1. Source("RNASEQAnalysis.R")
#               2. Run the function "rnaseq_parvathi" - specify raw count path, clinical info file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_RNASEQAnalysis.R/RNASEQAnalysis.R")
#               > rnaseq_parvathi(rCntPath = "C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_RAW_COUNTS.rda",
#                                 clinInfoPath_640 = "./data/coadread_tcga_clinical_data.tsv",
#                                 msiInfoPath = "./data/nationwidechildrens.org_auxiliary_coad_read.txt",
#                                 outputDir="./results/demographic/")
###

rnaseq_parvathi <- function(rCntPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_RAW_COUNTS.rda",
                            clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                            msiInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_auxiliary_coad_read.txt",
                            outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/demographic/") {
  
  
  
  
}
