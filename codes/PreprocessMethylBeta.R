###
#   File name : PreprocessMethylBeta.R
#   Author    : Hyunjin Kim
#   Date      : May 29, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load beta files, combine them, remove duplicated samples, remove NA rows,
#               remove probes that lie on sex chromosomes, remove cross-reactive probes,
#               and perform quality control analysis
#
#   Workflow
#               https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
#
#   Instruction
#               1. Source("PreprocessMethylBeta.R")
#               2. Run the function "preprocessMethylBeta" - specify the necessary input file paths and output directory
#               3. The combined and normalized methylation data will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PreprocessMethylBeta.R/PreprocessMethylBeta.R")
#               > preprocessMethylBeta(betaFileDir="./data/methylation/",
#                                      sampleInfoPath="./data/methylation/gdc_sample_sheet_coad_read_methylation.tsv",
#                                      xReactiveProbePath="E:/Methylation/cross_reactive_probes/xReactiveProbes_450k.txt",
#                                      outputDir="./data/methylation/preprocessed/")
###

### THE RAW METHYLATION BETA VALUES ARE TOO BIG (60GB), SO I DID NOT UPLOAD THE DATA ON THE SERVER AND RAN THIS CODE ON MY LOCAL MACHINE
preprocessMethylBeta <- function(betaFileDir="./data/methylation/",
                                 sampleInfoPath="./data/methylation/gdc_sample_sheet_coad_read_methylation.tsv",
                                 xReactiveProbePath="E:/Methylation/cross_reactive_probes/xReactiveProbes_450k.txt",
                                 outputDir="./data/methylation/preprocessed/") {
  
  ### load libraries
  if(!require(ChAMP)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ChAMP")
    library(ChAMP)
  }
  
  ### save current working directory
  ### ChAMP sometimes set working directory on its own
  currentWD <- getwd()
  
  ### load beta file paths
  f <- list.files(betaFileDir, pattern = "jhu-usc.edu_*")
  
  ### load beta values and combine them
  methyl_beta_tcga_coad_read <- read.table(file = paste0(betaFileDir, f[1]), header = TRUE,
                                           sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  CpGInfo <- methyl_beta_tcga_coad_read[,-2]
  rownames(CpGInfo) <- CpGInfo$`Composite Element REF`
  methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[,"Beta_value"]
  for(fileName in f[-1]) {
    methyl_beta_tcga_coad_read <- data.frame(methyl_beta_tcga_coad_read,
                                             read.table(file = paste0(betaFileDir, fileName), header = TRUE,
                                                        sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)[,"Beta_value"],
                                             stringsAsFactors = FALSE, check.names = FALSE)
  }
  colnames(methyl_beta_tcga_coad_read) <- sapply(f, function(x) strsplit(x, ".", fixed = TRUE)[[1]][6], USE.NAMES = FALSE)
  rownames(methyl_beta_tcga_coad_read) <- CpGInfo[,"Composite Element REF"]
  
  ### load sample info
  methyl_beta_sample_info <- read.table(file = sampleInfoPath, header = TRUE, sep = "\t",
                                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(methyl_beta_sample_info) <- sapply(methyl_beta_sample_info$`File Name`, function(x) strsplit(x, ".", fixed = TRUE)[[1]][6], USE.NAMES = FALSE)
  
  ### order(colnames(beta_value)) = order(rownames(sample_info))
  methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[,rownames(methyl_beta_sample_info)]
  
  ### only keep primary and recurrent tumors
  keep_idx <- union(which(methyl_beta_sample_info$`Sample Type` == "Primary Tumor"),
                    which(methyl_beta_sample_info$`Sample Type` == "Recurrent Tumor"))
  methyl_beta_sample_info <- methyl_beta_sample_info[keep_idx,]
  methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[,keep_idx]
  
  ### remove duplicated samples (1:15 of TCGA barcode, 16 means vials)
  methyl_beta_sample_info$`Sample ID` <- substr(methyl_beta_sample_info$`Sample ID`, 1, 15)
  remove_idx <- which(duplicated(methyl_beta_sample_info$`Sample ID`))
  if(length(remove_idx) > 0) {
    methyl_beta_sample_info <- methyl_beta_sample_info[-remove_idx,]
    methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[,-remove_idx]
  }
  
  ### change column names as sample ID
  colnames(methyl_beta_tcga_coad_read) <- methyl_beta_sample_info$`Sample ID`
  
  ### remove all [NA across all samples] rows
  remove_idx <- which(apply(apply(methyl_beta_tcga_coad_read, 1, is.na), 2, sum) == ncol(methyl_beta_tcga_coad_read))
  if(length(remove_idx) > 0) {
    methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[-remove_idx,]
    CpGInfo <- CpGInfo[-remove_idx,]
  }
  
  ### remove cross-reactive probes
  xReactiveProbes <- as.character(read.table(xReactiveProbePath, header = FALSE, sep = "\t")[,1])
  keep_idx <- which(!(rownames(methyl_beta_tcga_coad_read) %in% xReactiveProbes))
  methyl_beta_tcga_coad_read <- methyl_beta_tcga_coad_read[keep_idx,]
  CpGInfo <- CpGInfo[keep_idx,]
  
  ### make pd (sample sheet) matrix as an input for ChAMP
  pd <- data.frame(Sample_Name=colnames(methyl_beta_tcga_coad_read),
                   Sample_Plate=substr(rownames(methyl_beta_sample_info),22,25),
                   Sample_Group="CRC",
                   Pool_ID=NA,
                   Project=methyl_beta_sample_info$`Project ID`,
                   Sample_Well=NA,
                   Slide=NA,
                   Array=NA,
                   stringsAsFactors = FALSE, check.names = FALSE)
  rownames(pd) <- pd$Sample_Name
  
  ### filter the beta matrix with ChAMP
  filteredB <- champ.filter(beta = methyl_beta_tcga_coad_read, pd = pd,
                    fixOutlier = FALSE, arraytype = "450K")
  CpGInfo <- CpGInfo[rownames(filteredB$beta),]
  
  ### impute beta values only for the NAs (BMIQ does not take NAs)
  filteredB <- champ.impute(beta = as.matrix(filteredB$beta), pd = filteredB$pd)
  CpGInfo <- CpGInfo[rownames(filteredB$beta),]
  
  ### the champ.impute() changed the original working directory
  setwd(currentWD)
  
  ### generate QC plots with the filtered raw data
  png(paste0(outputDir, "raw_beta_qc_plots.png"), width = 2000, height = 1000, res = 120)
  par(mfrow=c(1,2))
  colors = rainbow(length(unique(filteredB$pd$Project)))
  names(colors) = unique(as.character(filteredB$pd$Project))
  plotMDS(filteredB$beta, top = 1000, pch = 19, col = colors[as.character(filteredB$pd$Project)],
          xlab = "Dimension1", ylab = "Dimension2", main = "MDS_Raw_Beta")
  legend("topright", legend = unique(as.character(filteredB$pd$Project)),
         col = colors[unique(as.character(filteredB$pd$Project))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  plot(density(as.numeric(filteredB$beta[,1])), main = "Density_Raw_Beta", ylim = c(0, 6),
       col = colors[as.character(filteredB$pd$Project[1])])
  for(i in 2:ncol(filteredB$beta)) {
    lines(density(as.numeric(filteredB$beta[,i])), col = colors[as.character(filteredB$pd$Project[i])])
  }
  legend("topright", legend = unique(as.character(filteredB$pd$Project)),
         col = colors[unique(as.character(filteredB$pd$Project))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  dev.off()
  
  ### BMIQ normalization
  normB <- champ.norm(beta = as.matrix(filteredB$beta), rgSet = NULL, mset = NULL,
                      resultsDir = outputDir, method = "BMIQ", plotBMIQ = FALSE,
                      arraytype = "450K", cores = 4)
  
  ### generate QC plots with the filtered raw data
  png(paste0(outputDir, "BMIQ_beta_qc_plots.png"), width = 2000, height = 1000, res = 120)
  par(mfrow=c(1,2))
  colors = rainbow(length(unique(filteredB$pd$Project)))
  names(colors) = unique(as.character(filteredB$pd$Project))
  plotMDS(normB, top = 1000, pch = 19, col = colors[as.character(filteredB$pd$Project)],
          xlab = "Dimension1", ylab = "Dimension2", main = "MDS_BMIQ_Beta")
  legend("topright", legend = unique(as.character(filteredB$pd$Project)),
         col = colors[unique(as.character(filteredB$pd$Project))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  plot(density(as.numeric(normB[,1])), main = "Density_BMIQ_Beta", ylim = c(0, 6),
       col = colors[as.character(filteredB$pd$Project[1])])
  for(i in 2:ncol(normB)) {
    lines(density(as.numeric(normB[,i])), col = colors[as.character(filteredB$pd$Project[i])])
  }
  legend("topright", legend = unique(as.character(filteredB$pd$Project)),
         col = colors[unique(as.character(filteredB$pd$Project))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  dev.off()
  
  ### save the analysis-ready data in RDA
  filteredB$beta <- normB
  normB <- filteredB
  save(list = c("normB", "CpGInfo"), file = paste0(outputDir, "norm_beta_tcga_coad_read.rda"))
  
  ### save the analysis-ready data in txt
  write.table(data.frame(CpG_site=rownames(normB$beta), normB$beta,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "norm_beta_tcga_coad_read.txt"), sep = "\t", row.names = FALSE)
  
}
