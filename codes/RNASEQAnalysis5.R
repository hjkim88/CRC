###
#   File name : RNASEQAnalysis5.R
#   Author    : Hyunjin Kim
#   Date      : Jul 2, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DE analysis on the TCGA COAD/READ raw counts between
#               different race groups. Also perform pathway analysis on the DE genes
#
#   * THIS CODE USES PREPROCESSED RAW COUNT DATA WHICH IS DIFFERENT THAN THE ONES USED IN
#     RNASEQAnalysis.R AND IN RNASEQAnalysis2.R.
#
#   * NOW POLE-MUTATED SAMPLES ARE REMOVED FROM THE ANALYSIS AND
#     LOCATION-ANALYSIS IS ALSO ADDED. THIS IS THE DIFFERENCE FROM RNASEQAnalysis3.R.
#
#   * MSI-L SAMPLES ARE NOW REGARDED AS MSS.
#
#   Instruction
#               1. Source("RNASEQAnalysis5.R")
#               2. Run the function "rnaseq_parvathi5" - specify raw count path, clinical info file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_RNASEQAnalysis5.R/RNASEQAnalysis5.R")
#               > rnaseq_parvathi5(rCntPath = "./data/rnaseq/raw_count_tcga_coad_read.rda",
#                                  clinInfoPath_640 = "./data/coadread_tcga_clinical_data_updated2.txt",
#                                  padj_thres = 0.05,
#                                  outputDir="./results/rnaseq/preprocessed_raw_counts3/")
###

rnaseq_parvathi5 <- function(rCntPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_RAW_COUNTS.rda",
                             clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt",
                             padj_thres = 0.05,
                             outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/rnaseq/") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(limma, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("limma", version = "3.8")
    require(limma, quietly = TRUE)
  }
  
  ### load the data
  load(rCntPath)
  clinicalInfo_640 <- read.table(file = clinInfoPath_640, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinicalInfo_640) <- clinicalInfo_640$`Sample ID`
  
  ### only retain info of the samples that have raw counts
  clinicalInfo_640 <- clinicalInfo_640[colnames(raw_count_tcga_coad_read),]
  
  ### add tissue location info to the clinical info
  clinicalInfo_640$TUMOR_LOCATION <- NA
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Cecum")] <- "Proximal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Ascending Colon")] <- "Proximal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Hepatic Flexure")] <- "Proximal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Transverse Colon")] <- "Proximal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Splenic Flexure")] <- "Distal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Descending Colon")] <- "Distal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Sigmoid Colon")] <- "Distal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Rectum")] <- "Distal"
  clinicalInfo_640$TUMOR_LOCATION[which(clinicalInfo_640$`Patient Primary Tumor Site` == "Rectosigmoid Junction")] <- "Distal"
  
  ### add new MSI info to the sample info since MSI-L should be treated as MSS
  clinicalInfo_640$NEW_MSI <- clinicalInfo_640$MSI
  clinicalInfo_640$NEW_MSI[which(clinicalInfo_640$NEW_MSI == "MSI-L")] <- "MSS"
  
  ### change other MSI-related info
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                            which(clinicalInfo_640$`Diagnosis Age` < 50))] <- "MSS_Young"
  clinicalInfo_640$MSI_AGE_Status[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                            which(clinicalInfo_640$`Diagnosis Age` >= 50))] <- "MSS_Old"
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                              which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSS_AA"
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                              which(clinicalInfo_640$`Race Category` == "WHITE"))] <- "MSS_CC"
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                              which(clinicalInfo_640$Prediction_Filtered == "African"))] <- "MSS_AA"
  clinicalInfo_640$MSI_RACE_Status2[intersect(which(clinicalInfo_640$MSI == "MSI-L"),
                                              which(clinicalInfo_640$Prediction_Filtered == "Caucasian"))] <- "MSS_CC"
  
  ### remove POLE-muated samples
  clinicalInfo_640 <- clinicalInfo_640[-which(clinicalInfo_640$POLE_MUTANT == TRUE),]
  raw_count_tcga_coad_read <- raw_count_tcga_coad_read[,rownames(clinicalInfo_640)]
  
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
    
    ### add baseMean for each group
    nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
    exp_rowMeans <- apply(nCnt[,which(Coldata$sampleType == exp_class), drop=FALSE], 1, mean)
    ctrl_rowMeans <- apply(nCnt[,which(Coldata$sampleType == ctrl_class), drop=FALSE], 1, mean)
    deresults <- data.frame(baseMean=deresults[,1],
                            V1=exp_rowMeans[rownames(deresults)],
                            V2=ctrl_rowMeans[rownames(deresults)],
                            deresults[,2:6],
                            stringsAsFactors = FALSE, check.names = FALSE)
    colnames(deresults)[2:3] <- c(paste0("baseMean_", exp_class), paste0("baseMean_", ctrl_class))
    
    return(deresults)
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
  
  ### a function to perform RNA-Seq analysis (MDS plot, DEA, volcano plot) with given input
  rnaseq_analysis <- function(target_idx, grp, exp_class, ctrl_class, fileName) {
    
    ### create a directory
    dir.create(paste0(outputDir, fileName))
    
    ### run DE analysis
    deresult <- deseqWithComparisons(rCnt = raw_count_tcga_coad_read[,target_idx], grp = grp,
                                     exp_class = exp_class, ctrl_class = ctrl_class,
                                     bat_eff = NULL, thresh = 1)
    
    ### write out the DE result table, draw a volcano plot, and perform pathway analysis
    write.xlsx2(data.frame(Gene_Symbol=rownames(deresult), deresult, stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir, fileName, "/deresult_", fileName, ".xlsx"),
                sheetName = fileName, row.names = FALSE)
    volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, fileName, "/volplot_", fileName, ".png"), pvalue = padj_thres)
    
    ### QC - MDS plots
    normCnt <- normalizeRNASEQwithVST(raw_count_tcga_coad_read[,target_idx])
    
    ### select the top genes for MDS using limma
    png(paste0(outputDir, fileName, "/mdsplot_", fileName, ".png"), width = 2000, height = 1000, res = 130)
    par(mfrow=c(1,2))
    colors = rainbow(length(unique(as.character(grp))))
    names(colors) = unique(as.character(grp))
    plotMDS(normCnt, top = 500, pch = 19, col = colors[as.character(grp)],
            xlab = "Dimension1", ylab = "Dimension2", main = paste0(fileName, "_500"))
    legend("topright", legend = unique(as.character(grp)),
           col = colors[unique(as.character(grp))], pch = 19,
           title = "Sample Groups", cex = 0.7)
    plotMDS(normCnt, top = 2000, pch = 19, col = colors[as.character(grp)],
            xlab = "Dimension1", ylab = "Dimension2", main = paste0(fileName, "_2000"))
    legend("topright", legend = unique(as.character(grp)),
           col = colors[unique(as.character(grp))], pch = 19,
           title = "Sample Groups", cex = 0.7)
    dev.off()
    
  }
  
  
  ### MSI-H: Young vs Old
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Young"),
                      which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Old"))
  grp <- clinicalInfo_640$MSI_AGE_Status[target_idx]
  grp <- sapply(grp, function(x) {
    if(grepl("-", x)) {
      return(paste(strsplit(x, "-", fixed = TRUE)[[1]], collapse = "_"))
    } else {
      return(x)
    }
  }, USE.NAMES = FALSE)
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSI_H_Young", "MSI_H_Old", "MSI-H_Young_vs_Old")
  
  
  ### MSS: Young vs Old
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_AGE_Status == "MSS_Young"),
                      which(clinicalInfo_640$MSI_AGE_Status == "MSS_Old"))
  grp <- clinicalInfo_640$MSI_AGE_Status[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSS_Young", "MSS_Old", "MSS_Young_vs_Old")
  
  
  ### MSI-H: AA vs CC
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_RACE_Status1 == "MSI-H_AA"),
                      which(clinicalInfo_640$MSI_RACE_Status1 == "MSI-H_CC"))
  grp <- clinicalInfo_640$MSI_RACE_Status1[target_idx]
  grp <- sapply(grp, function(x) {
    if(grepl("-", x)) {
      return(paste(strsplit(x, "-", fixed = TRUE)[[1]], collapse = "_"))
    } else {
      return(x)
    }
  }, USE.NAMES = FALSE)
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSI_H_AA", "MSI_H_CC", "MSI-H_AA_vs_CC")
  
  
  ### MSS: AA vs CC
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_RACE_Status1 == "MSS_AA"),
                      which(clinicalInfo_640$MSI_RACE_Status1 == "MSS_CC"))
  grp <- clinicalInfo_640$MSI_RACE_Status1[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSS_AA", "MSS_CC", "MSS_AA_vs_CC")
  
  
  ### MSI-H: AA vs CC based on the predicted values
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_AA"),
                      which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_CC"))
  grp <- clinicalInfo_640$MSI_RACE_Status2[target_idx]
  grp <- sapply(grp, function(x) {
    if(grepl("-", x)) {
      return(paste(strsplit(x, "-", fixed = TRUE)[[1]], collapse = "_"))
    } else {
      return(x)
    }
  }, USE.NAMES = FALSE)
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSI_H_AA", "MSI_H_CC", "MSI-H_AA_vs_CC_predicted")
  
  
  ### MSS: AA vs CC based on the predicted values
  ### set group info
  target_idx <- union(which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_AA"),
                      which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_CC"))
  grp <- clinicalInfo_640$MSI_RACE_Status2[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "MSS_AA", "MSS_CC", "MSS_AA_vs_CC_predicted")
  
  
  ### MSI-H & Young - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Young"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_Young_Proximal_vs_Distal")
  
  
  ### MSI-H & Old - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Old"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_Old_Proximal_vs_Distal")
  
  
  ### MSS & Young - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_AGE_Status == "MSS_Young"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_Young_Proximal_vs_Distal")
  
  
  ### MSS & Old - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_AGE_Status == "MSS_Old"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_Old_Proximal_vs_Distal")
  
  
  ### MSI-H & AA - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status1 == "MSI-H_AA"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_AA_Proximal_vs_Distal")
  
  
  ### MSI-H & CC - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status1 == "MSI-H_CC"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_CC_Proximal_vs_Distal")
  
  
  ### MSS & AA - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status1 == "MSS_AA"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_AA_Proximal_vs_Distal")
  
  
  ### MSS & CC - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status1 == "MSS_CC"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_CC_Proximal_vs_Distal")
  
  
  ### MSI-H & AA_predicted - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_AA"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_AA_predicted_Proximal_vs_Distal")
  
  
  ### MSI-H & CC_predicted - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status2 == "MSI-H_CC"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSI-H_CC_predicted_Proximal_vs_Distal")
  
  
  ### MSS & AA_predicted - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_AA"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_AA_predicted_Proximal_vs_Distal")
  
  
  ### MSS & CC_predicted - Proximal vs Distal
  ### set group info
  target_idx <- intersect(union(which(clinicalInfo_640$TUMOR_LOCATION == "Proximal"),
                                which(clinicalInfo_640$TUMOR_LOCATION == "Distal")),
                          which(clinicalInfo_640$MSI_RACE_Status2 == "MSS_CC"))
  grp <- clinicalInfo_640$TUMOR_LOCATION[target_idx]
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Proximal", "Distal", "MSS_CC_predicted_Proximal_vs_Distal")
  
  
  ### Global: Young vs Old
  ### set group info
  target_idx <- union(which(clinicalInfo_640$`Diagnosis Age` < 50),
                      which(clinicalInfo_640$`Diagnosis Age` >= 50))
  grp <- clinicalInfo_640$`Diagnosis Age`[target_idx]
  tIdx1 <- which(grp < 50)
  tIdx2 <- which(grp >= 50)
  grp[tIdx1] <- "Young"
  grp[tIdx2] <- "Old"
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "Young", "Old", "Young_vs_Old")
  
  
  ### Global: AA vs CC (self-reported)
  ### set group info
  target_idx <- union(which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"),
                      which(clinicalInfo_640$`Race Category` == "WHITE"))
  grp <- clinicalInfo_640$`Race Category`[target_idx]
  grp[which(grp == "BLACK OR AFRICAN AMERICAN")] <- "AA"
  grp[which(grp == "WHITE")] <- "CC"
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "AA", "CC", "AA_vs_CC")
  
  
  ### Global: AA vs CC (predicted)
  ### set group info
  target_idx <- union(which(clinicalInfo_640$Prediction_Filtered == "African"),
                      which(clinicalInfo_640$Prediction_Filtered == "Caucasian"))
  grp <- clinicalInfo_640$Prediction_Filtered[target_idx]
  grp[which(grp == "African")] <- "AA"
  grp[which(grp == "Caucasian")] <- "CC"
  ### run the rnaseq function
  rnaseq_analysis(target_idx, grp,  "AA", "CC", "AA_vs_CC_predicted")
  
}
