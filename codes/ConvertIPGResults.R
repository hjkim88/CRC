###
#   File name : ConvertIPGResults.R
#   Author    : Hyunjin Kim
#   Date      : Apr 27, 2019
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
#               > convertIPG(inputDir="./results/rnaseq/preprocessed_raw_counts/iPathwayGuide/")
###

convertIPG <- function(inputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Parvathi_Myer/results/rnaseq/iPathwayGuide/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  ### get a list of csv files
  csvFiles <- list.files(inputDir)
  csvFiles <- csvFiles[which(endsWith(csvFiles, ".csv"))]
  
  ### make new columns for the files
  mirnaCols <- c("miRNA_Name", "Down-regulated_Associated_DE_Genes", "All_Associated_DE_Genes", "Down-regulated_Associated_Genes", "All_Associated_Genes", "FDR")
  goTableCols <- c("GO_ID", "GO_Name", "Associated_DE_Genes", "All_Associated_Genes", "FDR")
  pathwayCols <- c("Pathway_Name", "FDR")
  diseaseCols <- c("Disease_Name", "Associated_DE_Genes", "All_Associated_Genes", "FDR")
  
  ### iteratively convert the csv files
  if(length(csvFiles) > 0) {
    for(i in 1:length(csvFiles)) {
      ### load a csv file
      csvFile <- read.csv(file = paste0(inputDir, csvFiles[i]))
      
      ### convert the columns
      if(endsWith(csvFiles[i], "mirnaTable.csv")) {
        colnames(csvFile) <- mirnaCols
      } else if(endsWith(csvFiles[i], "goTable.csv")) {
        colnames(csvFile) <- goTableCols
      } else if(endsWith(csvFiles[i], "pathwaysTable.csv")) {
        csvFile <- csvFile[,c(1,2)]
        colnames(csvFile) <- pathwayCols
      } else if(endsWith(csvFiles[i], "diseaseTable.csv")) {
        colnames(csvFile) <- diseaseCols
      } else {
        stop(paste0("The CSV file has unknown name: ", inputDir, csvFiles[i]))
      }
      
      ### save the result in Excel format
      write.xlsx2(csvFile, file = paste0(inputDir, substr(csvFiles[i], 1, nchar(csvFiles[i])-3), "xlsx"),
                  sheetName = substr(csvFiles[i], 1, nchar(csvFiles[i])-4), row.names = FALSE)
      
      ### erase the csv files
      file.remove(paste0(inputDir, csvFiles[i]))
    }
  }
  
}


# ### iteratively perform the function "convertIPG" with many directories
# directories <- list.dirs("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Zhenpeng/2018/DEA/drug_treatments/results/Pathway/")
# for(i in 2:length(directories)) {
#   convertIPG(inputDir = paste0(directories[i], "/"))
# }
