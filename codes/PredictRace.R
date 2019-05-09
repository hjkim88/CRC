###
#   File name : PredictRace.R
#   Author    : Hyunjin Kim
#   Date      : Apr 30, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Some samples in TCGA has no information of race. We use genotype data to
#               predict race of the samples. The 1000 genome project genotype data is
#               used as reference and EIGENSTRAT/PCoC are employed to get cluster info.
#
#   Instruction
#               1. Source("PredictRace.R")
#               2. Run the function "predict_race" - specify reference-test 012 matrix file and output result directory
#               3. The results will be generated in the output result directory
#
#   Example
#               > source("The_directory_of_PredictRace.R/PredictRace.R")
#               > predict_race(combinedDataPath="./data/genotyping/1000g_tcga_combined_genotype.RDA",
#                              outputResultDir="./results/genotyping/")
###

predict_race <- function(combinedDataPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/genotyping/1000g_tcga_combined_genotype.RDA",
                         outputResultDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/genotyping/") {
  
  ### load library
  if(!require(AssocTests, quietly = TRUE)) {
    install.packages("AssocTests")
    require(AssocTests, quietly = TRUE)
  }
  
  ### load the genotype data
  load(combinedDataPath)
  
  ### add Race column to the sample_info
  ### based on http://www.internationalgenome.org/category/population/
  sample_info <- data.frame(sample_info, Race1="Test_Samples", stringsAsFactors = FALSE, check.names = FALSE)
  sample_info$Race1[which(sample_info$Ethnicity %in% c("CHB", "JPT", "CHS", "CDX", "KHV"))] <- "East_Asian"
  sample_info$Race1[which(sample_info$Ethnicity %in% c("CEU", "TSI", "FIN", "GBR", "IBS"))] <- "Caucasian"
  sample_info$Race1[which(sample_info$Ethnicity %in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"))] <- "African"
  sample_info$Race1[which(sample_info$Ethnicity %in% c("MXL", "PUR", "CLM", "PEL"))] <- "Latino"
  sample_info$Race1[which(sample_info$Ethnicity %in% c("GIH", "PJL", "BEB", "STU", "ITU"))] <- "South_Asian"
  
  sample_info <- data.frame(sample_info, Race2="NOT REPORTED", stringsAsFactors = FALSE, check.names = FALSE)
  sample_info$Race2[which(sample_info$Ethnicity %in% c("CHB", "JPT", "CHS", "CDX", "KHV"))] <- "East_Asian"
  sample_info$Race2[which(sample_info$Ethnicity %in% c("CEU", "TSI", "FIN", "GBR", "IBS"))] <- "Caucasian"
  sample_info$Race2[which(sample_info$Ethnicity %in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"))] <- "African"
  sample_info$Race2[which(sample_info$Ethnicity %in% c("MXL", "PUR", "CLM", "PEL"))] <- "Latino"
  sample_info$Race2[which(sample_info$Ethnicity %in% c("GIH", "PJL", "BEB", "STU", "ITU"))] <- "South_Asian"
  sample_info$Race2[which(sample_info$Ethnicity == "WHITE")] <- "WHITE"
  sample_info$Race2[which(sample_info$Ethnicity == "BLACK OR AFRICAN AMERICAN")] <- "AFRICAN AMERICAN"
  sample_info$Race2[which(sample_info$Ethnicity == "ASIAN")] <- "ASIAN"
  sample_info$Race2[which(sample_info$Ethnicity == "AMERICAN INDIAN OR ALASKA NATIVE")] <- "AMERICAN INDIAN OR ALASKA NATIVE"
  
  # ### select rows that have top variances from a given matrix/data frame
  # selectTopV <- function(x, selectNum) {
  #   v <- apply(x, 1, var)
  #   x <- x[order(-v),]
  #   x <- x[1:selectNum,]
  #   
  #   return (x)
  # }
  
  ### A function to perform t-test for DEA
  ### mat = input matrix, row = observations, column = samples
  ### idx1 = column indices of condition1 (e.g., tumor)
  ### idx2 = column indices of condition2 (e.g., control)
  simple_t_test <- function(mat, idx1, idx2) {
    ### remove useless genes
    mat <- mat[rowSums(mat) > 1,]
    
    ### compuate raw mean
    mean <- apply(mat, 1, mean)
    
    ### normalization
    mat <- log2(mat+1)
    
    ### compute fold changes and t-statistic values
    lfc <- 0
    t <- 0
    for(i in 1:nrow(mat)) {
      lfc[i] <- mean(as.numeric(mat[i,idx1])) - mean(as.numeric(mat[i,idx2]))
      t[i] <- lfc[i] / sqrt((sd(mat[i,idx1])^2)/length(idx1) + (sd(mat[i,idx2])^2)/length(idx2))
    }
    
    ### compute p-value
    p <- 2*pt(-abs(t), df=length(union(idx1, idx2))-2)
    
    ### correct the p-value with Bonferroni correction
    bon <- 1-((1-p)^nrow(mat))
    
    ### A function to correct p-values with Benjamini-Hochberg approach
    correct_bh <- function(pvs) {
      
      temp <- cbind(p=pvs, No=seq(1:length(pvs)))
      
      if(length(which(is.na(temp[,"p"]))) > 0) {
        temp[which(is.na(temp[,"p"])),"p"] <- 1
      }
      
      temp <- temp[order(temp[,"p"]),]
      temp <- cbind(temp, Rank=seq(1:length(pvs)), BH=1)
      
      
      temp[length(pvs), "BH"] <- temp[length(pvs), "p"]
      for(i in (length(pvs)-1):1) {
        temp[i,"BH"] <- min(temp[i+1, "BH"], temp[i,"p"]*length(pvs)/temp[i,"Rank"])
      }
      
      temp <- temp[order(temp[,"No"]),]
      
      return(as.numeric(temp[,"BH"]))
    }  
    
    ### correct the p-value with Benjamini-Hochberg correction (FDR)
    bh <- correct_bh(p)
    
    ### combine mean, fold changes, t values, p-values, bonferroni, and benjamini-hochberg
    result <- data.frame(cbind(Gene=rownames(mat), mean=mean, log2FC=lfc, t=t, pVal=p, bon_pVal=bon, bh_pVal=bh))
    
    ### sort the result in order of bh p-value
    result <- result[order(result$bh_pVal, result$pVal),]
    
    ### if there are factor columns, change them to character columns
    idx <- sapply(result, is.factor)
    result[idx] <- lapply(result[idx], function(x) as.character(x))
    
    ### change numeric columns as numeric
    result[2:ncol(result)] <- lapply(result[2:ncol(result)], function(x) as.numeric(x))
    
    return(result)
  }
  
  ### remove rows that have very different values between TCGA and reference
  ### this is kind of a batch effect and we want to remove it
  ### perform t-test between the two and remove the rows with p-value < threshold
  ### idx1 = column indices of condition1
  ### idx2 = column indices of condition2
  filterConfoundingRows <- function(x, threshold, idx1, idx2) {
    ### run t-test
    result <- simple_t_test(x, idx1, idx2)
    
    ### filter the input matrix with the threshold
    result <- result[which(result$bon_pVal >= threshold),]
    x <- x[rownames(result),]
    
    return(x)
  }
  
  # ### downsize the genotype matrix with 2000 rows that have the top variance
  # geno <- selectTopV(combined_geno, 500)
  
  ### filter the genotype matrix to correct batch effect
  geno <- filterConfoundingRows(combined_geno, 0.01, 1:2504, 2505:3814)
  
  ### write out the 012 file - the eigenstrat function needs it
  ### we will delete the file later
  write.table(geno, file = "eigenstratG.eg.txt", quote = FALSE,
              sep = "", row.names = FALSE, col.names = FALSE)
  
  ### run eigenstrat on the data
  eigen_result <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = NULL, outFile.txt = NULL)$topK.eigenvectors
  pcoc_result <- pcoc(genoFile = "eigenstratG.eg.txt", outFile.txt = "pcoc.result.txt",
                      n.MonteCarlo = 100, num.splits = (length(unique(sample_info$Race1))-1))
  
  ### remove the temporary 012 file
  file.remove("eigenstratG.eg.txt")
  
  
  ### first plot labeling with Test_Samples
  png(paste0(outputResultDir, "PCA_1000g_TCGA_COAD_READ_Test_Samples.png"),
      width = 2200, height = 1200, res = 130)
  
  ### set the colors for labeling for each sample group
  # colors <- rainbow(length(unique(sample_info$Race1)), alpha = 0.5)
  colors <- c("red", "yellow", "green", "skyblue", "blue", "purple")
  names(colors) <- unique(sample_info$Race1)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### make a plot with Race1 - PC1 & PC2
  plot(eigen_result[,1], eigen_result[,2], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC2",
       xlab="Component1", ylab="Component2",
       col = colors[sample_info$Race1], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race1)[1:5], "_1000GP"), paste0(unique(sample_info$Race1)[6], "_TCGA")),
         col = colors[unique(sample_info$Race1)], pch = c(rep(19, 5), 15),
         title = "Sample Groups", cex = 0.8)
  
  ### make a plot with Race1 - PC1 & PC3
  plot(eigen_result[,1], eigen_result[,3], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC3",
       xlab="Component1", ylab="Component3",
       col = colors[sample_info$Race1], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race1)[1:5], "_1000GP"), paste0(unique(sample_info$Race1)[6], "_TCGA")),
         col = colors[unique(sample_info$Race1)], pch = c(rep(19, 5), 15),
         title = "Sample Groups", cex = 0.8)
  
  dev.off()
  
  
  ### second plot labeling with self-reported info
  png(paste0(outputResultDir, "PCA_1000g_TCGA_COAD_READ_Self_Reported.png"),
      width = 2200, height = 1200, res = 130)
  
  ### set the colors for labeling for each sample group
  # colors <- rainbow(length(unique(sample_info$Race2)), alpha = 0.5)
  colors <- c("red", "yellow", "green", "skyblue", "blue", "orange", "black", "pink", "purple", "gray")
  names(colors) <- unique(sample_info$Race2)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### make a plot with Race2 - PC1 & PC2
  plot(eigen_result[,1], eigen_result[,2], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC2",
       xlab="Component1", ylab="Component2",
       col = colors[sample_info$Race2], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race2)[1:5], "_1000GP"), paste0(unique(sample_info$Race2)[6:10], "_TCGA")),
         col = colors[unique(sample_info$Race2)], pch = c(rep(19, 5), rep(15, 5)),
         title = "Sample Groups", cex = 0.7)
  
  ### make a plot with Race2 - PC1 & PC3
  plot(eigen_result[,1], eigen_result[,3], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC3",
       xlab="Component1", ylab="Component3",
       col = colors[sample_info$Race2], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race2)[1:5], "_1000GP"), paste0(unique(sample_info$Race2)[6:10], "_TCGA")),
         col = colors[unique(sample_info$Race2)], pch = c(rep(19, 5), rep(15, 5)),
         title = "Sample Groups", cex = 0.7)
  
  dev.off()
  
  
  ### third plot labeling with Test_Samples (only with AA and CC)
  png(paste0(outputResultDir, "PCA_1000g_TCGA_COAD_READ_Test_Samples2.png"),
      width = 2200, height = 1200, res = 130)
  
  ### set the colors for labeling for each sample group
  colors <- c("red", "blue", "purple")
  names(colors) <- c("Caucasian", "African", "Test_Samples")
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### make a plot with Race1 - PC1 & PC2
  plot(eigen_result[,1], eigen_result[,2], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC2",
       xlab="Component1", ylab="Component2",
       col = colors[sample_info$Race1], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c("Caucasian_1000GP", "African_1000GP", "Test_Samples_TCGA"),
         col = c("red", "blue", "purple"), pch = c(19, 19, 15),
         title = "Sample Groups", cex = 0.8)
  
  ### make a plot with Race1 - PC1 & PC3
  plot(eigen_result[,1], eigen_result[,3], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC3",
       xlab="Component1", ylab="Component3",
       col = colors[sample_info$Race1], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c("Caucasian_1000GP", "African_1000GP", "Test_Samples_TCGA"),
         col = c("red", "blue", "purple"), pch = c(19, 19, 15),
         title = "Sample Groups", cex = 0.8)
  
  dev.off()
  
  
  ### fourth plot labeling with self-reported info (only with AA and CC)
  png(paste0(outputResultDir, "PCA_1000g_TCGA_COAD_READ_Self_Reported2.png"),
      width = 2200, height = 1200, res = 130)
  
  ### set the colors for labeling for each sample group
  colors <- c("red", "blue", "orange", "black", "purple")
  names(colors) <- c("Caucasian", "African", "WHITE", "AFRICAN AMERICAN", "NOT REPORTED")
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### make a plot with Race2 - PC1 & PC2
  plot(eigen_result[,1], eigen_result[,2], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC2",
       xlab="Component1", ylab="Component2",
       col = colors[sample_info$Race2], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race2)[c(1, 5)], "_1000GP"), paste0(unique(sample_info$Race2)[c(6, 7, 9)], "_TCGA")),
         col = c("red", "blue", "orange", "black", "purple"), pch = c(19, 19, 15, 15, 15),
         title = "Sample Groups", cex = 0.7)
  
  ### make a plot with Race2 - PC1 & PC3
  plot(eigen_result[,1], eigen_result[,3], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC3",
       xlab="Component1", ylab="Component3",
       col = colors[sample_info$Race2], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = c(paste0(unique(sample_info$Race2)[c(1, 5)], "_1000GP"), paste0(unique(sample_info$Race2)[c(6, 7, 9)], "_TCGA")),
         col = c("red", "blue", "orange", "black", "purple"), pch = c(19, 19, 15, 15, 15),
         title = "Sample Groups", cex = 0.7)
  
  dev.off()
  
  
  ### race prediction based on the cluster of 1000GP
  
  ### calculate cluster means (or medians)
  medians <- vector("list", (length(unique(sample_info$Race1)) - 1))
  names(medians) <- unique(sample_info$Race1)[1:length(medians)]
  for(i in 1:length(medians)) {
    ### mean of 3 dimensions - 3 points per row - the number of components PC1, PC2, and PC3
    medians[[i]] <- apply(eigen_result[which(sample_info$Race1 == names(medians)[i]),1:3], 2, mean)
  }
  
  ### add prediction result columns to the sample info
  sample_info <- data.frame(sample_info, Dist_CC=NA, Dist_EA=NA, Dist_LT=NA, Dist_SA=NA, Dist_AF=NA, Prediction=NA,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### Euclidean distance function
  euc.dist <- function(x1, x2) sqrt(sum((x1-x2)^2))
  
  ### predict the race of the test samples (TCGA)
  testIdx <- which(sample_info$Race1 == "Test_Samples")
  for(idx in testIdx) {
    ### calculate Euclidean distance between samples and medians
    ### only using PC1, PC2, and PC3
    sample_info$Dist_CC[idx] <- euc.dist(eigen_result[idx,1:3], medians[[1]])
    sample_info$Dist_EA[idx] <- euc.dist(eigen_result[idx,1:3], medians[[2]])
    sample_info$Dist_LT[idx] <- euc.dist(eigen_result[idx,1:3], medians[[3]])
    sample_info$Dist_SA[idx] <- euc.dist(eigen_result[idx,1:3], medians[[4]])
    sample_info$Dist_AF[idx] <- euc.dist(eigen_result[idx,1:3], medians[[5]])
    
    ### predict race of the samples
    sample_info$Prediction[idx] <- names(medians)[which(sample_info[idx,] == min(sample_info[idx,8:12]))-7]
  }
  
  ### accuracy
  compareIdx <- union(union(which(sample_info$Race2 == "WHITE"),
                            which(sample_info$Race2 == "AFRICAN AMERICAN")),
                      which(sample_info$Race2 == "ASIAN"))
  
  acc <- sapply(compareIdx, function(x) {
    if(sample_info$Race2[x] == "WHITE" && sample_info$Prediction[x] == "Caucasian") {
      return(TRUE)
    } else if(sample_info$Race2[x] == "AFRICAN AMERICAN" && sample_info$Prediction[x] == "African") {
      return(TRUE)
    } else if(sample_info$Race2[x] == "ASIAN" && sample_info$Prediction[x] == "East_Asian") {
      return(TRUE)
    } else if(sample_info$Race2[x] == "ASIAN" && sample_info$Prediction[x] == "South_Asian") {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  writeLines(paste("Accuracy (Caucasian, African, Asian) = ", (length(which(acc))/length(acc))*100, "%"))
  writeLines(paste(length(which(acc)), "cases are correct out of", length(acc), "cases"))
  
  ### write out the result table
  write.table(data.frame(sample_info[testIdx,c(3,4)], Race=sample_info$Race2[testIdx],
                         sample_info[testIdx,8:13], stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputResultDir, "Race_Prediction_Result.txt"),
              sep = "\t", row.names = FALSE)
  
  ### draw a plot with the predicted results
  png(paste0(outputResultDir, "PCA_1000g_TCGA_COAD_READ_Prediction.png"),
      width = 2200, height = 1200, res = 130)
  
  ### labels
  lbls <- c(paste0(sample_info$Race2[1:2504], "_1000GP"), paste0(sample_info$Prediction[2505:3814], "_TCGA"))
  
  ### set the colors for labeling for each sample group
  colors <- c("red", "yellow", "green", "skyblue", "blue", "orange", "gray", "pink", "purple", "black")
  names(colors) <- unique(lbls)
  
  ### two plots in one file
  par(mfrow=c(1,2))
  
  ### make a plot with Race2 - PC1 & PC2
  plot(eigen_result[,1], eigen_result[,2], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC2",
       xlab="Component1", ylab="Component2",
       col = colors[lbls], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = unique(lbls),
         col = colors[unique(lbls)], pch = c(rep(19, 5), rep(15, 5)),
         title = "Sample Groups", cex = 0.7)
  
  ### make a plot with Race2 - PC1 & PC3
  plot(eigen_result[,1], eigen_result[,3], main="1000 Genomes Project + TCGA COAD/READ - PC1 & PC3",
       xlab="Component1", ylab="Component3",
       col = colors[lbls], pch = c(rep(19, 2504), rep(15, 1310)))
  legend("topright", legend = unique(lbls),
         col = colors[unique(lbls)], pch = c(rep(19, 5), rep(15, 5)),
         title = "Sample Groups", cex = 0.7)
  
  dev.off()
  
}
