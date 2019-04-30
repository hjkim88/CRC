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
#                                 padj_thres = 0.05,
#                                 outputDir="./results/rnaseq/")
###

rnaseq_parvathi <- function(rCntPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_RAW_COUNTS.rda",
                            clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                            msiInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_auxiliary_coad_read.txt",
                            padj_thres = 0.05,
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
  if(length(which(duplicated(clinicalInfo_640$`Patient ID`))) > 0) {
    clinicalInfo_640 <- clinicalInfo_640[-which(duplicated(clinicalInfo_640$`Patient ID`)),]
  }
  rownames(clinicalInfo_640) <- clinicalInfo_640$`Patient ID`
  msiInfo <- read.table(file = msiInfoPath, sep = "\t",
                        header = TRUE, row.names = 1,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
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
  
  #******************************************************************************************
  # Generate MDS plot for matrix columns.
  #
  # * mat:		matrix whose columns are to be clustered.
  # * plot_names:	if TRUE, name the points on the MDS plot.
  # * alt_names:	String vector. By default, if plot_names == TRUE, the plotted points are 
  #		named using colnames(mat). However, if alt_names in not NULL, then this is used
  #		instead. In that case, length(alt_names) should be equal to ncol(mat) and the point
  #		mat[, i] will be named after alt_names[i].
  # * groups:		if not NULL, it must be a string vector such that names(groups) == colnames(mat)
  #		and groups[i] is the name of the group where the i-th column of "mat" belongs. This
  #		info will be used to color members of the same group using the same color. If the
  #		value of this argument is NULL, then no group-based coloring is performed.
  # * dist_fun:	distance function to use for computing distances between the column vectors in
  #		mat. This should take as input a matrix object and return an object of class "dist". If
  #		the value of this argument is NULL, then the value of the argument dist_options below
  # 		is used for determining how to compute the distances.
  # * dist_options:	if dist_fun is NULL the the standard "dist" function in R is used for 
  #		computing distances. In that case, the value of dist_options specifies which distance 
  #		option to use when calling "dist". The default option is "euclidan".
  # * save:		if TRUE, save plot to file. Otherwise just plot on screen.
  # * f_name:		if save == TRUE, this is the full pathname of the file were the plot will be saved.
  # * width, height, res:	values of graphical parameters to use when generating the plot.
  # * xlab, ylab, main:	titles for x-axis, y-axis, and entire plot, respectively
  # * pch:		plot parameter, specifying the plot point type. The default pch = 1 drwaws open circles.
  #		Another useful setting is pch = 19, this draws full circles.
  mdsPlot <- function(mat, plot_names = FALSE, alt_names = NULL, groups = NULL, dist_fun = NULL,
                      dist_options = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"), 
                      save = FALSE, f_name = "./mds_plot.png", width = 1000, height = 1000, res = 130,
                      xlab = "", ylab="", main="", pch = 1){
    
    ### load library
    if(!require(ArgumentCheck, quietly = TRUE)) {
      install.packages("ArgumentCheck")
      require(ArgumentCheck, quietly = TRUE)
    }
    if(!require(randomcoloR, quietly = TRUE)) {
      install.packages("randomcoloR")
      require(randomcoloR, quietly = TRUE)
    }
    
    ### argument checking
    check <- ArgumentCheck::newArgCheck()
    if(is.null(mat)) {
      ArgumentCheck::addError(
        msg = "[mat]: matrix whose columns are to be clustered should exists",
        argcheck = check
      )
    }
    if((!is.null(alt_names)) && (!length(alt_names) == ncol(mat))) {
      ArgumentCheck::addError(
        msg = "[alt_names]: length(alt_names) should be equal to ncol(mat)",
        argcheck = check
      )
    }
    if((!is.null(groups)) && (!length(groups) == ncol(mat))) {
      ArgumentCheck::addError(
        msg = "[groups]: length(groups) should be equal to ncol(mat)",
        argcheck = check
      )
    }
    if((!is.null(dist_fun)) && (!class(dist_fun(mat)) == "dist")) {
      ArgumentCheck::addError(
        msg = "[dist_fun]: it should be a function which generates \"dist\" object",
        argcheck = check
      )
    }
    if((is.null(dist_fun)) && (is.null(dist_options))) {
      ArgumentCheck::addError(
        msg = "[dist_options]: if the \"[dist_fun]\" was not provided, \"[dist_options]\" should be provided and should be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\", \"pearson\", \"spearman\", or \"kendall\"",
        argcheck = check
      )
    } else if((is.null(dist_fun)) && (!dist_options %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"))) {
      ArgumentCheck::addError(
        msg = "[dist_options]: if the \"[dist_fun]\" was not provided, \"[dist_options]\" should be provided and should be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\", \"pearson\", \"spearman\", or \"kendall\"",
        argcheck = check
      )
    }
    if((save == TRUE) && is.null(f_name)) {
      ArgumentCheck::addError(
        msg = "[f_name]: if \"[save]\" == TRUE, then \"[f_name]\" should be provided",
        argcheck = check
      )
    }
    ArgumentCheck::finishArgCheck(check)
    
    
    ### save the mds plot as a png format
    if(save){
      png(f_name, width = width, height = height, res = res)
    }
    
    ### make a distance matrix
    if(!is.null(dist_fun)) {
      d <- dist_fun(mat)
    } else if(dist_options[1] == "euclidean") {
      d <- dist(t(mat), method = "euclidean")
    } else if(dist_options[1] == "maximum") {
      d <- dist(t(mat), method = "maximum")
    } else if(dist_options[1] == "manhattan") {
      d <- dist(t(mat), method = "manhattan")
    } else if(dist_options[1] == "canberra") {
      d <- dist(t(mat), method = "canberra")
    } else if(dist_options[1] == "binary") {
      d <- dist(t(mat), method = "binary")
    } else if(dist_options[1] == "minkowski") {
      d <- dist(t(mat), method = "minkowski")
    } else if(dist_options[1] == "pearson") {
      d <- as.dist(1-cor(mat, method = "pearson"))
    } else if(dist_options[1] == "spearman") {
      d <- as.dist(1-cor(mat, method = "spearman"))
    } else if(dist_options[1] == "kendall") {
      d <- as.dist(1-cor(mat, method = "kendall"))
    }
    
    ### get MDS points
    fit <- cmdscale(d,eig=TRUE, k=2)
    Dimension1 <- fit$points[,1]
    Dimension2 <- fit$points[,2]
    
    ### group coloring
    if(is.null(groups)) {
      ### make a MDS plot
      plot(Dimension1, Dimension2, main=main, xlab=xlab, ylab=ylab)
      
      ### print point names
      if(plot_names == TRUE) {
        if(is.null(alt_names)) {
          text(Dimension1, Dimension2, labels = labels(d), cex=.7, pos=3)
        } else {
          text(Dimension1, Dimension2, labels = alt_names, cex=.7, pos=3)
        }
      }  
    } else {
      ###
      if(length(unique(as.character(groups))) < 6)
        colors = c("black", "red", "blue", "magenta", "green")[1:length(unique(as.character(groups)))]
      else {
        require(randomcoloR)
        colors = distinctColorPalette(length(unique(as.character(groups))))
        # colors = rainbow(length(unique(as.character(groups))))
      }
      names(colors) = unique(as.character(groups))
      
      ### make a MDS plot
      plot(Dimension1, Dimension2, main=main, xlab=xlab, ylab=ylab,
           col = colors[as.character(groups)], pch = pch)
      legend("topright", legend = unique(as.character(groups)),
             col = colors[unique(as.character(groups))], pch = 15,
             title = "Sample Groups", cex = 0.7)
      
      ### print point names
      if(plot_names == TRUE) {
        if(is.null(alt_names)) {
          text(Dimension1, Dimension2, labels = labels(d), cex=0.7, pos=3, col = colors[as.character(groups)])
        } else {
          text(Dimension1, Dimension2, labels = alt_names, cex=0.7, pos=3, col = colors[as.character(groups)])
        }
      }
    }
    
    ### print out the plot
    if(save)
      dev.off()
    
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
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("clusterProfiler")
      library(clusterProfiler)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### colect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title))
              
              png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[1]])
              dev.off()
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title))
              
              png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[2]])
              dev.off()
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  ### gene mapping list
  ensembl2eg <- unlist(as.list(org.Hs.egENSEMBL2EG))
  eg2gs <- unlist(as.list(org.Hs.egSYMBOL))
  
  
  ### MSI-H: Young vs Old
  
  ### extract raw counts of samples of our interest and run DE analysis
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Young"),which(clinicalInfo_640$MSI_AGE_Status == "MSI-H_Old"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  deresult <- deseqWithComparisons(rCnt = rCnt, grp = grp,
                                   exp_class = "Young", ctrl_class = "Old",
                                   bat_eff = NULL, thresh = 1)
  
  ### annotate the genes with gene symbols
  entrez_id <- ensembl2eg[sapply(rownames(deresult), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])]
  deresult <- data.frame(Ensembl_ID=rownames(deresult),
                         Gene_Symbol=eg2gs[entrez_id], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result table, draw a volcano plot, and perform pathway analysis
  fileName <- "msi-h_young_vs_old"
  write.xlsx2(deresult, file = paste0(outputDir, "deresult_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "volplot_", fileName, ".png"), pvalue = padj_thres)
  pathresult <- pathwayAnalysis_CP(geneList = entrez_id[which(deresult$padj < padj_thres)], org = "human", database = "GO",
                                   displayNum = 50, title = paste0("pathway_", fileName),
                                   pv_threshold = padj_thres, dir = outputDir)
  write.xlsx2(pathresult, file = paste0(outputDir, "pathway_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  
  ### QC - MDS plots
  normCnt <- normalizeRNASEQwithVST(rCnt)
  ### original MDS
  mdsPlot(normCnt, groups = grp, save = TRUE, pch = 19,
          main = fileName, xlab = "Dimension1", ylab = "Dimension2",
          f_name = paste0(outputDir, "mdsplot_", fileName, ".png"))
  
  ### select the top genes for MDS using limma
  ### age
  png(paste0(outputDir, "mdsplot2_", fileName, "_age.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### race
  png(paste0(outputDir, "mdsplot2_", fileName, "_race.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### gender
  png(paste0(outputDir, "mdsplot2_", fileName, "_gender.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"Sex"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### stage
  png(paste0(outputDir, "mdsplot2_", fileName, "_stage.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"American Joint Committee on Cancer Tumor Stage Code"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  
  
  ### MSS: Young vs Old
  
  ### extract raw counts of samples of our interest and run DE analysis
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_AGE_Status == "MSS_Young"),which(clinicalInfo_640$MSI_AGE_Status == "MSS_Old"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  deresult <- deseqWithComparisons(rCnt = rCnt, grp = grp,
                                   exp_class = "Young", ctrl_class = "Old",
                                   bat_eff = NULL, thresh = 1)
  
  ### annotate the genes with gene symbols
  entrez_id <- ensembl2eg[sapply(rownames(deresult), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])]
  deresult <- data.frame(Ensembl_ID=rownames(deresult),
                         Gene_Symbol=eg2gs[entrez_id], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result table, draw a volcano plot, and perform pathway analysis
  fileName <- "mss_young_vs_old"
  write.xlsx2(deresult, file = paste0(outputDir, "deresult_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "volplot_", fileName, ".png"), pvalue = padj_thres)
  pathresult <- pathwayAnalysis_CP(geneList = entrez_id[which(deresult$padj < padj_thres)], org = "human", database = "GO",
                                   displayNum = 50, title = paste0("pathway_", fileName),
                                   pv_threshold = padj_thres, dir = outputDir)
  write.xlsx2(pathresult, file = paste0(outputDir, "pathway_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  
  ### QC - MDS plots
  normCnt <- normalizeRNASEQwithVST(rCnt)
  ### original MDS
  mdsPlot(normCnt, groups = grp, save = TRUE, pch = 19,
          main = fileName, xlab = "Dimension1", ylab = "Dimension2",
          f_name = paste0(outputDir, "mdsplot_", fileName, ".png"))
  
  ### select the top genes for MDS using limma
  ### age
  png(paste0(outputDir, "mdsplot2_", fileName, "_age.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### race
  png(paste0(outputDir, "mdsplot2_", fileName, "_race.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### gender
  png(paste0(outputDir, "mdsplot2_", fileName, "_gender.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"Sex"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### stage
  png(paste0(outputDir, "mdsplot2_", fileName, "_stage.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"American Joint Committee on Cancer Tumor Stage Code"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  
  
  ### MSI-H: AA vs CC
  
  ### extract raw counts of samples of our interest and run DE analysis
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status == "MSI-H_AA"),which(clinicalInfo_640$MSI_RACE_Status == "MSI-H_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  deresult <- deseqWithComparisons(rCnt = rCnt, grp = grp,
                                   exp_class = "AA", ctrl_class = "CC",
                                   bat_eff = NULL, thresh = 1)
  
  ### annotate the genes with gene symbols
  entrez_id <- ensembl2eg[sapply(rownames(deresult), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])]
  deresult <- data.frame(Ensembl_ID=rownames(deresult),
                         Gene_Symbol=eg2gs[entrez_id], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result table, draw a volcano plot, and perform pathway analysis
  fileName <- "msi-h_AA_vs_CC"
  write.xlsx2(deresult, file = paste0(outputDir, "deresult_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "volplot_", fileName, ".png"), pvalue = padj_thres)
  pathresult <- pathwayAnalysis_CP(geneList = entrez_id[which(deresult$padj < padj_thres)], org = "human", database = "GO",
                                   displayNum = 50, title = paste0("pathway_", fileName),
                                   pv_threshold = padj_thres, dir = outputDir)
  write.xlsx2(pathresult, file = paste0(outputDir, "pathway_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  
  ### QC - MDS plots
  normCnt <- normalizeRNASEQwithVST(rCnt)
  ### original MDS
  mdsPlot(normCnt, groups = grp, save = TRUE, pch = 19,
          main = fileName, xlab = "Dimension1", ylab = "Dimension2",
          f_name = paste0(outputDir, "mdsplot_", fileName, ".png"))
  
  ### select the top genes for MDS using limma
  ### age
  png(paste0(outputDir, "mdsplot2_", fileName, "_age.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### race
  png(paste0(outputDir, "mdsplot2_", fileName, "_race.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### gender
  png(paste0(outputDir, "mdsplot2_", fileName, "_gender.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"Sex"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### stage
  png(paste0(outputDir, "mdsplot2_", fileName, "_stage.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"American Joint Committee on Cancer Tumor Stage Code"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  
  
  ### MSS: AA vs CC
  
  ### extract raw counts of samples of our interest and run DE analysis
  rCnt <- htseq_raw_counts[,which(tcga_sample_info$`Case ID` %in% clinicalInfo_640$`Patient ID`[union(which(clinicalInfo_640$MSI_RACE_Status == "MSS_AA"),which(clinicalInfo_640$MSI_RACE_Status == "MSS_CC"))])]
  grp <- clinicalInfo_640[tcga_sample_info[colnames(rCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  deresult <- deseqWithComparisons(rCnt = rCnt, grp = grp,
                                   exp_class = "AA", ctrl_class = "CC",
                                   bat_eff = NULL, thresh = 1)
  
  ### annotate the genes with gene symbols
  entrez_id <- ensembl2eg[sapply(rownames(deresult), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])]
  deresult <- data.frame(Ensembl_ID=rownames(deresult),
                         Gene_Symbol=eg2gs[entrez_id], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result table, draw a volcano plot, and perform pathway analysis
  fileName <- "mss_AA_vs_CC"
  write.xlsx2(deresult, file = paste0(outputDir, "deresult_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "volplot_", fileName, ".png"), pvalue = padj_thres)
  pathresult <- pathwayAnalysis_CP(geneList = entrez_id[which(deresult$padj < padj_thres)], org = "human", database = "GO",
                                   displayNum = 50, title = paste0("pathway_", fileName),
                                   pv_threshold = padj_thres, dir = outputDir)
  write.xlsx2(pathresult, file = paste0(outputDir, "pathway_", fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
  
  ### QC - MDS plots
  normCnt <- normalizeRNASEQwithVST(rCnt)
  ### original MDS
  mdsPlot(normCnt, groups = grp, save = TRUE, pch = 19,
          main = fileName, xlab = "Dimension1", ylab = "Dimension2",
          f_name = paste0(outputDir, "mdsplot_", fileName, ".png"))
  
  ### select the top genes for MDS using limma
  ### age
  png(paste0(outputDir, "mdsplot2_", fileName, "_age.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_AGE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### race
  png(paste0(outputDir, "mdsplot2_", fileName, "_race.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"MSI_RACE_Status"]
  grp <- sapply(grp, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### gender
  png(paste0(outputDir, "mdsplot2_", fileName, "_gender.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"Sex"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
  ### stage
  png(paste0(outputDir, "mdsplot2_", fileName, "_stage.png"), width = 2000, height = 1000, res = 130)
  par(mfrow=c(1,2))
  grp <- clinicalInfo_640[tcga_sample_info[colnames(normCnt),"Case ID"],"American Joint Committee on Cancer Tumor Stage Code"]
  if(length(which(is.na(grp))) > 0) {
    grp[which(is.na(grp))] <- "NA"
  }
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
