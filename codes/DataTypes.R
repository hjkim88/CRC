###
#   File name : DataTypes.R
#   Author    : Hyunjin Kim
#   Date      : Apr 10, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Identify data types of all available files on TCGA portal
#               that are from the samples of our interest (the 640 samples)
#
#   Instruction
#               1. Source("DataTypes.R")
#               2. Run the function "identifyDataTypes" - specify the file path of the 640 sample info
#                  and the file path of the sample info of all the TCGA COAD-READ
#               3. The results will be printed in the console
#
#   Example
#               > source("The_directory_of_DataTypes.R/DataTypes.R")
#               > identifyDataTypes(clinInfoPath = "./data/coadread_tcga_clinical_data.tsv",
#                                   tcgaFileInfoPath = "./data/TCGA_COAD_READ_All_Files_Info.txt")
###

identifyDataTypes <- function(clinInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                              tcgaFileInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data///isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/TCGA_COAD_READ_All_Files_Info.txt") {
  
  
  
  
}


