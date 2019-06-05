###
#   File name : DifferentialMethylation.R
#   Author    : Hyunjin Kim
#   Date      : May 31, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find differentially methylated probes and regions,
#               perform pathway analysis on differentially methylated genes.
#
#   Workflow
#               https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
#
#   Instruction
#               1. Source("DifferentialMethylation.R")
#               2. Run the function "dma" - specify the necessary input file paths and output directory
#               3. The differentially methylated results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DifferentialMethylation.R/DifferentialMethylation.R")
#               > dma(preprocessedBetaPath="./data/methylation/preprocessed/norm_beta_tcga_coad_read.rda",
#                     clinInfoPath_640 = "./data/coadread_tcga_clinical_data.tsv",
#                     msiInfoPath = "./data/nationwidechildrens.org_auxiliary_coad_read.txt",
#                     predictedRaceInfoPath="./results/genotyping/Race_Prediction_Result.txt",
#                     pvalThreshold = 0.05,
#                     cpg_cutoff = 10,
#                     dmrPrintNum = 3,
#                     outputDir="./results/methylation/")
###

dma <- function(preprocessedBetaPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/norm_beta_tcga_coad_read.rda",
                clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                msiInfoPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_auxiliary_coad_read.txt",
                predictedRaceInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/genotyping/Race_Prediction_Result.txt",
                pvalThreshold = 0.05,
                cpg_cutoff = 10,
                dmrPrintNum = 3,
                dmrPrintSampleNum = 20,
                outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/methylation/") {
  
  ### load libraries
  if(!require(ChAMP)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ChAMP")
    library(ChAMP)
  }
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  ### load the data
  load(preprocessedBetaPath)
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
  
  ### add group info (MSI/MSS - African American/Caucasian) - based on self-reported info
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI_RACE_Status1=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                              which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSI-H_AA"
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSI-H"),
                                              which(clinicalInfo_640$`Race Category` == "WHITE"))] <- "MSI-H_CC"
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                              which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSS_AA"
  clinicalInfo_640$MSI_RACE_Status1[intersect(which(clinicalInfo_640$MSI == "MSS"),
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
              
              png(paste0(dir, "kegg_", title, ".png"), width = 2000, height = 1000)
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
              
              png(paste0(dir, "go_", title, ".png"), width = 2000, height = 1000)
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
  gs2eg <- unlist(as.list(org.Hs.egSYMBOL2EG))
  
  ### set sample groups for DMA - age
  grp <- clinicalInfo_640[substr(colnames(normB$beta), 1, 12),"MSI_AGE_Status"]
  grp[which(grp == "MSI-H_Young")] <- "MSIHYOUNG"
  grp[which(grp == "MSI-H_Old")] <- "MSIHOLD"
  grp[which(grp == "MSS_Young")] <- "MSSYOUNG"
  grp[which(grp == "MSS_Old")] <- "MSSOLD"
  grp[which(is.na(grp))] <- "NOTHING"
  
  ### Young vs Old
  ### MSI-H
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSIHYOUNG", "MSIHOLD"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSI-H_Young_vs_Old.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSI-H_Young_vs_Old",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSI-H_Young_vs_Old.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### MSS
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSSYOUNG", "MSSOLD"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSS_Young_vs_Old.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSS_Young_vs_Old",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSS_Young_vs_Old.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### DMRCate and plots
  ### make a design matrix
  design <- model.matrix(~0+grp)
  colnames(design) <- levels(as.factor(grp))
  ### make a contrast matrix
  contrastMat <- makeContrasts(MSIHYOUNG-MSIHOLD,
                               MSSYOUNG-MSSOLD,
                               levels = design)
  ### MSI-H
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSIHYOUNG - MSIHOLD")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSI-H_Young_vs_Old.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSIHYOUNG"), which(grp == "MSIHOLD"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSIHYOUNG", "MSIHOLD")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSI-H_Young_vs_Old.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  ### MSS
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSSYOUNG - MSSOLD")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSS_Young_vs_Old.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSSYOUNG"), which(grp == "MSSOLD"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSSYOUNG", "MSSOLD")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSS_Young_vs_Old.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  
  
  ### set sample groups for DMA - self-reported race
  grp <- clinicalInfo_640[substr(colnames(normB$beta), 1, 12),"MSI_RACE_Status1"]
  grp[which(grp == "MSI-H_AA")] <- "MSIHAA"
  grp[which(grp == "MSI-H_CC")] <- "MSIHCC"
  grp[which(grp == "MSS_AA")] <- "MSSAA"
  grp[which(grp == "MSS_CC")] <- "MSSCC"
  grp[which(is.na(grp))] <- "NOTHING"
  
  ### AA vs CC - self-reported
  ### MSI-H
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSIHAA", "MSIHCC"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSI-H_AA_vs_CC.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSI-H_AA_vs_CC",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSI-H_AA_vs_CC.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### MSS
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSSAA", "MSSCC"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSS_AA_vs_CC.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSS_AA_vs_CC",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSS_AA_vs_CC.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### DMRCate and plots
  ### make a design matrix
  design <- model.matrix(~0+grp)
  colnames(design) <- levels(as.factor(grp))
  ### make a contrast matrix
  contrastMat <- makeContrasts(MSIHAA-MSIHCC,
                               MSSAA-MSSCC,
                               levels = design)
  ### MSI-H
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSIHAA - MSIHCC")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSI-H_AA_vs_CC.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSIHAA"), which(grp == "MSIHCC"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSIHAA", "MSIHCC")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSI-H_AA_vs_CC.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  ### MSS
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSSAA - MSSCC")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSS_AA_vs_CC.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSSAA"), which(grp == "MSSCC"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSSAA", "MSSCC")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSS_AA_vs_CC.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  
  
  ### set sample groups for DMA - predicted race
  grp <- clinicalInfo_640[substr(colnames(normB$beta), 1, 12),"MSI_RACE_Status2"]
  grp[which(grp == "MSI-H_AA")] <- "MSIHAA"
  grp[which(grp == "MSI-H_CC")] <- "MSIHCC"
  grp[which(grp == "MSS_AA")] <- "MSSAA"
  grp[which(grp == "MSS_CC")] <- "MSSCC"
  grp[which(is.na(grp))] <- "NOTHING"
  
  ### AA vs CC - predicted
  ### MSI-H
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSIHAA", "MSIHCC"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSI-H_AA_vs_CC_predicted.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSI-H_AA_vs_CC_predicted",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSI-H_AA_vs_CC_predicted.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### MSS
  ### differentially methylated positions
  dmps <- champ.DMP(beta = normB$beta, pheno = grp, compare.group = c("MSSAA", "MSSCC"),
                    adjPVal = pvalThreshold, adjust.method = "BH", arraytype = "450K")[[1]]
  write.xlsx2(data.frame(CpG_Site=rownames(dmps), dmps),
              file = paste0(outputDir, "DMPs_MSS_AA_vs_CC_predicted.xlsx"),
              sheetName = "DMPs", row.names = FALSE)
  ### pathway analysis with the DMPs
  dm_genes <- as.character(dmps$gene)
  dm_genes <- dm_genes[which(dm_genes != "")]
  dm_genes <- gs2eg[dm_genes]
  dm_genes <- dm_genes[which(!is.na(dm_genes))]
  pathways <- pathwayAnalysis_CP(geneList = dm_genes, org = "human", database = "GO", imgPrint = TRUE,
                                 title = "Top_50_DMG-associated_Pathways_MSS_AA_vs_CC_predicted",
                                 displayNum = 50, dir = outputDir)
  write.xlsx2(pathways, file = paste0(outputDir, "go_DMG-associated_Pathways_MSS_AA_vs_CC_predicted.xlsx"),
              sheetName = "DMG-associated_Pathways", row.names = FALSE)
  ### DMRCate and plots
  ### make a design matrix
  design <- model.matrix(~0+grp)
  colnames(design) <- levels(as.factor(grp))
  ### make a contrast matrix
  contrastMat <- makeContrasts(MSIHAA-MSIHCC,
                               MSSAA-MSSCC,
                               levels = design)
  ### MSI-H
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSIHAA - MSIHCC")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSI-H_AA_vs_CC_predicted.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSIHAA"), which(grp == "MSIHCC"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSIHAA", "MSIHCC")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSI-H_AA_vs_CC_predicted.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  ### MSS
  ### annotation for DMR(Differentially Methylated Region)s
  myAnno <- cpg.annotate(object = normB$beta, datatype = "array", what = "Beta",
                         arraytype = "450K", fdr = pvalThreshold,
                         annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "MSSAA - MSSCC")
  ### get DMRs and save
  if(length(which(myAnno$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs <- dmrcate(myAnno, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges <- extractRanges(DMRs, genome = "hg19")
    ### filter with Stouffer combined p-values and save the DMR info
    if(length(results.ranges) > 0) {
      results.ranges <- results.ranges[which(results.ranges$Stouffer < pvalThreshold),]
      results.ranges <- results.ranges[which(results.ranges$no.cpgs > cpg_cutoff)]
      write.xlsx2(data.frame(DMR=paste0(rep("DMR", length(results.ranges)), 1:length(results.ranges)), results.ranges),
                  file = paste0(outputDir, "DMRs_MSS_AA_vs_CC_predicted.xlsx"),
                  sheetName = "DMRs", row.names = FALSE)
      idx <- union(which(grp == "MSSAA"), which(grp == "MSSCC"))
      pheno <- c("skyblue", "pink")
      names(pheno) <- c("MSSAA", "MSSCC")
      cols <- pheno[grp]
      min_num <- min(dmrPrintSampleNum, min(length(which(cols[idx] == pheno[1])), length(which(cols[idx] == pheno[2]))))
      set.seed(1234)
      for(i in 1:min(dmrPrintNum, length(results.ranges))) {
        png(paste0(outputDir, "DMR", i, "_MSS_AA_vs_CC_predicted.png"), width = 1800, height = 1000)
        DMR.plot(ranges=results.ranges, dmr=i, CpGs=normB$beta[,idx], phen.col=cols[idx], what = "Beta",
                 arraytype = "450K", pch=19, toscale=TRUE, plotmedians=TRUE, genome="hg19",
                 samps = union(sample(which(cols[idx] == pheno[1]), min_num),
                               sample(which(cols[idx] == pheno[2]), min_num)))
        dev.off()
      }
    }
  }
  
}
