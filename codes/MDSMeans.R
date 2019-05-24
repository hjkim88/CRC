###
#   File name : MDSMeans.R
#   Author    : Hyunjin Kim
#   Date      : May 14, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : We found there are many DE genes between MSS_AA and MSS_CC from DE analysis.
#               But when we looked at MDS plots, there was no big difference between MSS_AA and MSS_CC.
#               Looking by eye can tell us there is no big difference but it actually can be.
#               So, we want to show mean and median value of each cluster in the MDS plots for
#               more exact comparison between two clusters.
#
#   Instruction
#               1. Source("MDSMeans.R")
#               2. Run the function "mds_mean" - specify raw count path, clinical info file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MDSMeans.R/MDSMeans.R")
#               > mds_mean(rCntPath = "C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_RAW_COUNTS.rda",
#                          clinInfoPath_640 = "./data/coadread_tcga_clinical_data.tsv",
#                          msiInfoPath = "./data/nationwidechildrens.org_auxiliary_coad_read.txt",
#                          predictedRaceInfoPath="./results/genotyping/Race_Prediction_Result.txt",
#                          outputDir="./results/rnaseq/")
###

mds_mean <- function(rCntPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_RAW_COUNTS.rda",
                     clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                     msiInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_auxiliary_coad_read.txt",
                     predictedRaceInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/genotyping/Race_Prediction_Result.txt",
                     outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/rnaseq/") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  ### load the data
  load(rCntPath)
  clinicalInfo_640 <- read.table(file = clinInfoPath_640, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
  if(length(which(duplicated(clinicalInfo_640$`Patient ID`))) > 0) {
    clinicalInfo_640 <- clinicalInfo_640[-which(duplicated(clinicalInfo_640$`Patient ID`)),]
  }
  rownames(clinicalInfo_640) <- clinicalInfo_640$`Patient ID`
  msiInfo <- read.table(file = msiInfoPath, sep = "\t",
                        header = TRUE, row.names = 1,
                        stringsAsFactors = FALSE, check.names = FALSE)
  predicted_race_info <- read.table(file = predictedRaceInfoPath, header = TRUE, sep = "\t",
                                    row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  predicted_race_info <- data.frame(Sample_ID=substr(rownames(predicted_race_info), 1, 15), predicted_race_info,
                                    stringsAsFactors = FALSE, check.names = FALSE)
  predicted_race_info <- predicted_race_info[!duplicated(predicted_race_info$Sample_ID),]
  
  ### retain only primary tumor & new primary tumor
  tcga_sample_info <- tcga_sample_info[union(union(which(tcga_sample_info[,"Sample Type"] == "Primary Tumor"),
                                                   which(tcga_sample_info[,"Sample Type"] == "Additional - New Primary")),
                                             which(tcga_sample_info[,"Sample Type"] == "Primary Blood Derived Cancer - Peripheral Blood")),]
  
  ### remove FFPE samples
  tcga_sample_info <- tcga_sample_info[which(tcga_sample_info[,"is_derived_from_ffpe"] == "NO"),]
  
  ### raw counts with the filtered samples
  htseq_raw_counts <- htseq_raw_counts[,rownames(tcga_sample_info)]
  
  ### add msi status
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  common_cases <- intersect(rownames(clinicalInfo_640), rownames(msiInfo))
  clinicalInfo_640[common_cases, "MSI"] <- msiInfo[common_cases, "mononucleotide_and_dinucleotide_marker_panel_analysis_status"]
  
  ### add predicted race info (from the Eigenstrat result)
  clinicalInfo_640 <- merge(clinicalInfo_640, predicted_race_info[,c("Sample_ID","Prediction_Filtered")],
                            by.x = "Sample ID", by.y = "Sample_ID", all.x = TRUE)
  rownames(clinicalInfo_640) <- clinicalInfo_640$`Patient ID`
  
  ### add group info (MSI/MSS - Young/Old)
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI_AGE_Status=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                            which(clinicalInfo_640$`Diagnosis Age` < 50))] <- "MSI-H_Young"
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                            which(clinicalInfo_640$`Diagnosis Age` >= 50))] <- "MSI-H_Old"
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                            which(clinicalInfo_640$`Diagnosis Age` < 50))] <- "MSS_Young"
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                            which(clinicalInfo_640$`Diagnosis Age` >= 50))] <- "MSS_Old"
  
  ### add group info (MSI/MSS - African American/Caucasian)
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI_RACE_Status=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_640$MSI_RACE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                             which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSI-H_AA"
  clinicalInfo_640$MSI_RACE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                             which(clinicalInfo_640$`Race Category` == "WHITE"))] <- "MSI-H_CC"
  clinicalInfo_640$MSI_RACE_Status[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                             which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSS_AA"
  clinicalInfo_640$MSI_RACE_Status[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                             which(clinicalInfo_640$`Race Category` == "WHITE"))] <- "MSS_CC"
  
  ### add group info (MSI/MSS - African/Caucasian) - based on the predicted info
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI_RACE_Status2=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                             which(clinicalInfo_640$Prediction_Filtered == "African"))] <- "MSI-H_AA"
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                             which(clinicalInfo_640$Prediction_Filtered == "Caucasian"))] <- "MSI-H_CC"
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                             which(clinicalInfo_640$Prediction_Filtered == "African"))] <- "MSS_AA"
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                             which(clinicalInfo_640$Prediction_Filtered == "Caucasian"))] <- "MSS_CC"
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  
  ### WITH SELF-REPORTED AGE INFO
  
  ### MSI-H: Young vs Old
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Young"),which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Old"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "msi-h_young_vs_old"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 235")
  
  ### calculate means for each group cluster
  young_mean <- c(mean(fit$points[which(grp == "Young"),1]), mean(fit$points[which(grp == "Young"),2]))
  old_mean <- c(mean(fit$points[which(grp == "Old"),1]), mean(fit$points[which(grp == "Old"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("Young_Mean", "Old_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(young_mean[1], young_mean[2], col = colors["Young"], pch = 15, cex = 2)
  points(old_mean[1], old_mean[2], col = colors["Old"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(young_mean, old_mean))), digits = 5)))
  
  dev.off()
  
  
  ### MSS: Young vs Old
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_AGE_Status == "MSS_Young"),which(clinicalInfo_640$MSI_AGE_Status == "MSS_Old"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "mss_young_vs_old"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 652")
  
  ### calculate means for each group cluster
  young_mean <- c(mean(fit$points[which(grp == "Young"),1]), mean(fit$points[which(grp == "Young"),2]))
  old_mean <- c(mean(fit$points[which(grp == "Old"),1]), mean(fit$points[which(grp == "Old"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("Young_Mean", "Old_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(young_mean[1], young_mean[2], col = colors["Young"], pch = 15, cex = 2)
  points(old_mean[1], old_mean[2], col = colors["Old"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(young_mean, old_mean))), digits = 5)))
  
  dev.off()
  
  
  
  ### WITH SELF-REPORTED RACE INFO
  
  ### MSI-H: AA vs CC
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status == "MSI-H_AA"),which(clinicalInfo_640$MSI_RACE_Status == "MSI-H_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "msi-h_AA_vs_CC"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 63")
  
  ### calculate means for each group cluster
  aa_mean <- c(mean(fit$points[which(grp == "AA"),1]), mean(fit$points[which(grp == "AA"),2]))
  cc_mean <- c(mean(fit$points[which(grp == "CC"),1]), mean(fit$points[which(grp == "CC"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  dev.off()
  
  
  ### MSS: AA vs CC
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status == "MSS_AA"),which(clinicalInfo_640$MSI_RACE_Status == "MSS_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "mss_AA_vs_CC"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 3769")
  
  ### calculate means for each group cluster
  aa_mean <- c(mean(fit$points[which(grp == "AA"),1]), mean(fit$points[which(grp == "AA"),2]))
  cc_mean <- c(mean(fit$points[which(grp == "CC"),1]), mean(fit$points[which(grp == "CC"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  dev.off()
  
  
  
  ### WITH PREDICTED RACE INFO
  
  ### MSI-H: AA vs CC
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_AA"),which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status2"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "msi-h_AA_vs_CC_predicted"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 77")
  
  ### calculate means for each group cluster
  aa_mean <- c(mean(fit$points[which(grp == "AA"),1]), mean(fit$points[which(grp == "AA"),2]))
  cc_mean <- c(mean(fit$points[which(grp == "CC"),1]), mean(fit$points[which(grp == "CC"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  dev.off()
  
  
  ### MSS: AA vs CC
  
  ### extract raw counts of samples of our interest and set group labels
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_AA"),which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status2"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  
  ### plot file name
  fileName <- "mss_AA_vs_CC_predicted"
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### set colors for groups
  colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(grp)))]
  names(colors) = unique(as.character(grp))
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 6384")
  
  ### calculate means for each group cluster
  aa_mean <- c(mean(fit$points[which(grp == "AA"),1]), mean(fit$points[which(grp == "AA"),2]))
  cc_mean <- c(mean(fit$points[which(grp == "CC"),1]), mean(fit$points[which(grp == "CC"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  dev.off()
  
  
  ### remove samples in the right side of the PCA plot
  rIdx <- which(fit$points[,1] > 25)
  plot(fit$points[-rIdx,1], fit$points[-rIdx,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp[-rIdx])], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  plot(fit$points[-rIdx,1], fit$points[-rIdx,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  aa_mean <- c(mean(fit$points[-rIdx,1][which(grp[-rIdx] == "AA")]), mean(fit$points[-rIdx,2][which(grp[-rIdx] == "AA")]))
  cc_mean <- c(mean(fit$points[-rIdx,1][which(grp[-rIdx] == "CC")]), mean(fit$points[-rIdx,2][which(grp[-rIdx] == "CC")]))
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  rIdx <- union(which(fit$points[,1] > 0), which(fit$points[,2] < 10))
  plot(fit$points[-rIdx,1], fit$points[-rIdx,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp[-rIdx])], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  plot(fit$points[-rIdx,1], fit$points[-rIdx,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  aa_mean <- c(mean(fit$points[-rIdx,1][which(grp[-rIdx] == "AA")]), mean(fit$points[-rIdx,2][which(grp[-rIdx] == "AA")]))
  cc_mean <- c(mean(fit$points[-rIdx,1][which(grp[-rIdx] == "CC")]), mean(fit$points[-rIdx,2][which(grp[-rIdx] == "CC")]))
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  ### remove the samples from the raw counts
  rCnt <- rCnt[,-rIdx]
  grp <- grp[-rIdx]
  
  ### normalize the raw counts with VST
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  ### get PCA eigen values
  d <- dist(t(normCnt), method = "euclidean")
  fit <- cmdscale(d,eig=TRUE, k=2)
  
  ### plot file name
  fileName <- "mss_AA_vs_CC_predicted"
  
  
  #####################################################
  ### A function to perform DE analysis with DESeq2 ###
  #####################################################
  #' @title deseqWithComparisons
  #' @param rCnt raw count matrix
  #' @param grp a character vector of class info of the samples
  #' @param exp_class a string of the experiment group's name
  #' @param ctrl_class a string of the control group's name
  #' @param bat_eff a character vector of batch effect info of the samples
  #' @param thresh numeric. Filters out from the results genes with adjusted
  #' 		p-value larger than this value
  #' @return data.frame
  #' @export
  #' @author Hyunjin Kim
  ####################################################
  deseqWithComparisons <- function(rCnt, grp, exp_class, ctrl_class, bat_eff=NULL, thresh = 1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DE analysis
    sampleType <- as.character(grp)
    
    if(is.null(bat_eff)) {
      Coldata <- data.frame(sampleType)
    } else {
      batch_eff <- as.character(bat_eff)
      Coldata <- data.frame(sampleType, batch_eff)
    }
    
    rownames(Coldata) <- colnames(rCnt)
    Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
    
    ### data preparation for DE analysis
    if(is.null(bat_eff)) {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
    } else {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
    }
    
    deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### run DE analysis
    dea <- DESeq(deSeqData)
    deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    deresults <- deresults[order(deresults$padj, na.last = TRUE), ,drop = FALSE]
    deresults <- deresults[deresults$padj <= thresh, ,drop = FALSE]
    
    return(data.frame(deresults))
  }
  
  ### A function to print volcano plot of DE analysis with DESeq2 result
  volPlotWithDeseq <- function(deresult, outputFilePath, pvalue=0.05) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    deresult$padj[which(is.na(deresult$padj))] <- 1
    volcanoData <- as.data.frame(cbind(deresult$log2FoldChange, -log10(deresult$padj), as.character(as.factor(deresult$padj < pvalue))))
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    s <- strsplit(outputFilePath, "/")
    f1 <- s[[1]][length(s[[1]])]
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) + xlab("log2_Fold_Change") + ylab("-log10(FDR)") + ggtitle(paste("Significant (adj.p < ", pvalue, " ) DE genes -", sum(volcanoData$Significance == "TRUE"))) + theme_classic(base_size = 16) + geom_point(alpha=0.4)
    ggsave(filename = outputFilePath)
  }
  
  ### DE analysis without the samples in the right side of the PCA plot
  deresult <- deseqWithComparisons(rCnt = rCnt, grp = grp,
                                   exp_class = "AA", ctrl_class = "CC",
                                   bat_eff = NULL, thresh = 1)
  
  ### gene mapping list
  ensembl2eg <- unlist(as.list(org.Hs.egENSEMBL2EG))
  eg2gs <- unlist(as.list(org.Hs.egSYMBOL))
  
  ### annotate the genes with gene symbols
  entrez_id <- ensembl2eg[sapply(rownames(deresult), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])]
  deresult <- data.frame(Ensembl_ID=rownames(deresult),
                         Gene_Symbol=eg2gs[entrez_id], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result table, draw a volcano plot
  write.xlsx2(deresult, file = paste0(outputDir, "deresult_", fileName, "_left_only.xlsx"), sheetName = fileName, row.names = FALSE)
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "volplot_", fileName, "._left_only.png"), pvalue = 0.05)
  
  ### make a MDS plot
  png(paste0(outputDir, "mdsplot_", fileName, "_means_left_only.png"), width = 2000, height = 1000, res = 130)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### draw the original PCA
  plot(fit$points[,1], fit$points[,2], main=fileName, xlab="Dimension1", ylab="Dimension2",
       col = colors[as.character(grp)], pch = 19)
  legend("topright", legend = unique(as.character(grp)),
         col = colors[unique(as.character(grp))], pch = 19,
         title = "Sample Groups", cex = 0.7)
  
  ### the number of DE genes
  mtext("The number of DE genes: 381")
  
  ### calculate means for each group cluster
  aa_mean <- c(mean(fit$points[which(grp == "AA"),1]), mean(fit$points[which(grp == "AA"),2]))
  cc_mean <- c(mean(fit$points[which(grp == "CC"),1]), mean(fit$points[which(grp == "CC"),2]))
  
  ### print the original PCA with "white" - keep the plot scale 
  plot(fit$points[,1], fit$points[,2], main=paste0(fileName, "_means"),
       xlab="Dimension1", ylab="Dimension2",
       col = "white", pch = 19)
  legend("topright", legend = c("CC_Mean", "AA_Mean"),
         col = c("black", "red"), pch = 15,
         title = "Sample Groups", cex = 0.7)
  
  ### print the cluster means
  points(aa_mean[1], aa_mean[2], col = colors["AA"], pch = 15, cex = 2)
  points(cc_mean[1], cc_mean[2], col = colors["CC"], pch = 15, cex = 2)
  
  ### print the distance between the two means
  mtext(paste("Euclidean distance between the two means:",
              signif(as.numeric(dist(rbind(aa_mean, cc_mean))), digits = 5)))
  
  dev.off()
  
  ###
  
  
  
  
  
}
