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
                              tcgaFileInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/TCGA_COAD_READ_All_Files_Info.txt") {
  
  ### load the data
  clinInfo <- read.table(file = clinInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  tcgaFileInfo <- read.table(file = tcgaFileInfoPath, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)
  
  ### there are multiple sample IDs in some of tcgaFileInfo, so be aware of that
  
  ### get indicies in tcgaFileInfo that overlap with the sample of interest
  interest_idx <- NULL
  for(i in 1:nrow(tcgaFileInfo)) {
    for(j in 1:nrow(clinInfo)) {
      if(grepl(clinInfo$`Sample ID`[j], tcgaFileInfo$`Sample ID`[i])) {
        interest_idx <- c(interest_idx, i)
        break
      }
    }
    if(i %% 100 == 0) {
      writeLines(paste(i, "/", nrow(tcgaFileInfo)))
    }
  }
  
  ### only retain the info of the 640 samples
  info_interest <- tcgaFileInfo[interest_idx,]
  
  ### get data types of all the available files
  total_data_types <- unique(info_interest$`Data Type`)
  
  ### print out the data types
  for(i in 1:length(total_data_types)) {
    writeLines(paste("# of samples for", total_data_types[i], "=", length(which(info_interest$`Data Type` == total_data_types[i]))))
  }
  
  ### Gene Expression Quantification - there are unnormalized & normalized, so we should distinguish them
  writeLines(paste("# of samples for RNA-Seq counts =",
                   length(intersect(which(info_interest$`Data Type` == "Gene Expression Quantification"),
                                    grep("counts", info_interest$`File Name`)))))
  
  ### Aligned Reads - WXS, RNA-Seq, and miRNA-Seq
  writeLines(paste("# of samples for RNA-Seq reads =",
                   length(intersect(which(info_interest$`Data Type` == "Aligned Reads"),
                                    grep("rehead", info_interest$`File Name`)))))
  writeLines(paste("# of samples for miRNA-Seq reads =",
                   length(intersect(which(info_interest$`Data Type` == "Aligned Reads"),
                                    grep("mirna", info_interest$`File Name`)))))
  writeLines(paste("# of samples for WXS reads =",
                   length(which(info_interest$`Data Type` == "Aligned Reads"))-
                     length(intersect(which(info_interest$`Data Type` == "Aligned Reads"),
                                      grep("rehead", info_interest$`File Name`)))-
                     length(intersect(which(info_interest$`Data Type` == "Aligned Reads"),
                                      grep("mirna", info_interest$`File Name`)))))
}
