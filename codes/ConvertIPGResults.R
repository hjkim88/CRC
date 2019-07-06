###
#   File name : ConvertIPGResults.R
#   Author    : Hyunjin Kim
#   Date      : July 6, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Convert CSV files (iPathwayGuide Results) to Excel files
#
#   Instruction
#               1. Source("ConvertIPGResults.R")
#               2. Run the function "convertIPG" - specify the input directory (csv)
#               3. The converted files will be generated under the input directory
#
#   Example
#               > source("The_directory_of_ConvertIPGResults.R/ConvertIPGResults.R")
#               > convertIPG(inputDir="./results/rnaseq/preprocessed_raw_counts3/")
###

convertIPG <- function(inputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Parvathi_Myer/results/rnaseq/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  ### get a list of csv files
  csvFiles <- list.files(inputDir)
  csvFiles <- csvFiles[which(endsWith(csvFiles, ".csv"))]
  
  ### iteratively convert the csv files
  if(length(csvFiles) > 0) {
    ### make new columns for the files
    mirnaCols <- c("miRNA_Name", "Down-regulated_Associated_DE_Genes", "All_Associated_DE_Genes", "Down-regulated_Associated_Genes", "All_Associated_Genes", "FDR")
    goTableCols <- c("GO_ID", "GO_Name", "Associated_DE_Genes", "All_Associated_Genes", "FDR")
    pathwayCols <- c("Pathway_Name", "FDR")
    diseaseCols <- c("Disease_Name", "Associated_DE_Genes", "All_Associated_Genes", "FDR")
    
    for(i in 1:length(csvFiles)) {
      ### load a csv file
      csvFile <- read.csv(file = paste0(inputDir, csvFiles[i]))
      
      ### convert the columns
      if(endsWith(csvFiles[i], "mirnaTable.csv")) {
        colnames(csvFile) <- mirnaCols
        type <- strsplit(csvFiles[i], "_", fixed = TRUE)[[1]][4]
        type <- substr(type, 1, nchar(type)-4)
      } else if(endsWith(csvFiles[i], "goTable.csv")) {
        colnames(csvFile) <- goTableCols
        type <- strsplit(csvFiles[i], "_", fixed = TRUE)[[1]][4]
        type <- substr(type, 1, nchar(type)-4)
      } else if(endsWith(csvFiles[i], "pathwaysTable_fdr.csv")) {
        csvFile <- csvFile[,c(1,2)]
        colnames(csvFile) <- pathwayCols
        type <- strsplit(csvFiles[i], "_", fixed = TRUE)[[1]][4]
      } else if(endsWith(csvFiles[i], "diseaseTable.csv")) {
        colnames(csvFile) <- diseaseCols
        type <- strsplit(csvFiles[i], "_", fixed = TRUE)[[1]][4]
        type <- substr(type, 1, nchar(type)-4)
      } else {
        stop(paste0("The CSV file has unknown name: ", inputDir, csvFiles[i]))
      }
      
      ### save the result in Excel format
      write.xlsx2(csvFile, file = paste0(inputDir, "iPathwayGuide_", basename(inputDir), "_",
                                         type, ".xlsx"),
                  sheetName = basename(inputDir), row.names = FALSE)
      
      ### erase the csv files
      file.remove(paste0(inputDir, csvFiles[i]))
    }
  }
  
}


### iteratively perform the function "convertIPG" with many directories
directories <- list.dirs("./results/rnaseq/preprocessed_raw_counts3/")
for(i in 2:length(directories)) {
  convertIPG(inputDir = paste0(directories[i], "/"))
}
