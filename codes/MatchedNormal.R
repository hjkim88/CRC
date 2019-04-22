###
#   File name : MatchedNormal.R
#   Author    : Hyunjin Kim
#   Date      : Apr 11, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Identify how many normal Bam files exist in TCGA portal
#               that are matched to cancer Bam files
#
#   Instruction
#               1. Source("MatchedNormal.R")
#               2. Run the function "identifyNormals" - specify the file path of the 640 sample info
#                  and the file path of the TCGA WXS file info
#               3. The results will be printed in the console
#
#   Example
#               > source("The_directory_of_MatchedNormal.R/MatchedNormal.R")
#               > identifyNormals(clinInfoPath = "./data/coadread_tcga_clinical_data.tsv",
#                                 wxsInfoPath = "./data/TCGA_COAD_READ_WXS_Files_Info.tsv",
#                                 outputFilePath = "./data/MSIsensor/WXS_matched_tumor_normal_sample_info.txt")
###

identifyNormals <- function(clinInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                            wxsInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/TCGA_COAD_READ_WXS_Files_Info.tsv",
                            outputFilePath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/WXS_matched_tumor_normal_sample_info.txt") {
  
  ### load the data
  clinInfo <- read.table(file = clinInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  wxsInfo <- read.table(file = wxsInfoPath, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### get wxs info of samples of our interest (the 640 samples)
  wxsInfo_interest <- wxsInfo[which(wxsInfo$`Case ID` %in% clinInfo$`Patient ID`),]
  
  ### get info separately for cancer and normal
  wxsInfo_interest_cancer <- wxsInfo_interest[union(which(wxsInfo_interest$`Sample Type` == "Primary Tumor"),
                                                    which(wxsInfo_interest$`Sample Type` == "Recurrent Tumor")),]
  wxsInfo_interest_normal <- wxsInfo_interest[union(which(wxsInfo_interest$`Sample Type` == "Blood Derived Normal"),
                                                    which(wxsInfo_interest$`Sample Type` == "Solid Tissue Normal")),]
  
  ### print out the result
  matched_case <- intersect(wxsInfo_interest_cancer$`Case ID`, wxsInfo_interest_normal$`Case ID`)
  writeLines(paste("# of samples that have Bam files of the matched cancer-normal pair:",
                   length(matched_case)))
  
  ### make a table for the matched cancer and normal
  matched_info <- NULL
  for(i in 1:length(matched_case)) {
    ### get the indicies of the specific case
    a <- which(wxsInfo_interest_cancer$`Case ID` == matched_case[i])
    b <- which(wxsInfo_interest_normal$`Case ID` == matched_case[i])
    
    ### organize the info for the case
    temp <- c(wxsInfo_interest_cancer$`Case ID`[a[1]],
              wxsInfo_interest_cancer$`Project ID`[a[1]],
              paste(wxsInfo_interest_cancer$`Sample ID`[a], collapse = ","),
              paste(wxsInfo_interest_normal$`Sample ID`[b], collapse = ","),
              paste(wxsInfo_interest_cancer$`Sample Type`[a], collapse = ","),
              paste(wxsInfo_interest_normal$`Sample Type`[b], collapse = ","),
              paste(wxsInfo_interest_cancer$`File ID`[a], collapse = ","),
              paste(wxsInfo_interest_normal$`File ID`[b], collapse = ","),
              paste(wxsInfo_interest_cancer$`File Name`[a], collapse = ","),
              paste(wxsInfo_interest_normal$`File Name`[b], collapse = ","))
    matched_info <- rbind(matched_info, temp)
  }
  colnames(matched_info) <- c("Case_ID", "Project_ID", "Tumor_Sample_ID", "Normal_Sample_ID",
                              "Tumor_Sample_Type", "Normal_Sample_Type",
                              "Tumor_File_ID", "Normal_File_ID",
                              "Tumor_File_Name", "Normal_File_Name")
  
  ### write out the result
  write.table(matched_info, file = outputFilePath, sep = "\t", row.names = FALSE)
  
}
