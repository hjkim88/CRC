###
#   File name : SomaticMutationAnalysis2.R
#   Author    : Hyunjin Kim
#   Date      : Jun 25, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Using maftools on MAF file, perform somatic mutation analysis
#
#   * THIS IS DIFFERENCE FROM SomaticMutationAnalysis.R SINCE THAT ONE IS BASED ON VCF FILE,
#     AND THIS ONE IS BASED ON MAF FILE (+MAFTOOLS)
#   * Also, POLE-mutated samples were removed and MSI-L samples were treated as MSS
#     Comparisons based on location (Proximal/Distal) are also added
#
#   Instruction
#               1. Source("SomaticMutationAnalysis2.R")
#               2. Run the function "sm_analysis2" - specify the necessary input file paths and output directory
#               3. The analysis results will be generated in the output directory
#
#   Example
#               > source("The_directory_of_SomaticMutationAnalysis2.R/SomaticMutationAnalysis2.R")
#               > sm_analysis2(mafFilePath="./data/somatic_mutation/MAF/somatic_mutation_maf_tcga_coad_read.maf",
#                              sampleInfoPath="./data/coadread_tcga_clinical_data_updated2.txt",
#                              driver_gene_fdr_cutoff=0.05,
#                              outputDir="./results/somatic_mutation2/")
###

sm_analysis2 <- function(mafFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/somatic_mutation/somatic_mutation_maf_tcga_coad_read.maf",
                         sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt",
                         driver_gene_fdr_cutoff=0.05,
                         outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/somatic_mutation2/") {
  
  ### load libraries
  options(java.parameters = "-Xmx8000m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(maftools, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("maftools", version = "3.8")
    require(maftools, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  
  ### load clinical info
  clinInfo <- read.table(sampleInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinInfo) <- clinInfo$`Sample ID`
  
  ### remove POLE-mutated samples
  clinInfo <- clinInfo[-which(clinInfo$POLE_MUTANT == TRUE),]
  
  ### add tissue location info to the clinical info
  clinInfo$TUMOR_LOCATION <- NA
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Cecum")] <- "Proximal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Ascending Colon")] <- "Proximal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Hepatic Flexure")] <- "Proximal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Transverse Colon")] <- "Proximal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Splenic Flexure")] <- "Distal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Descending Colon")] <- "Distal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Sigmoid Colon")] <- "Distal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Rectum")] <- "Distal"
  clinInfo$TUMOR_LOCATION[which(clinInfo$`Patient Primary Tumor Site` == "Rectosigmoid Junction")] <- "Distal"
  
  ### add new MSI info to the sample info since MSI-L should be treated as MSS
  clinInfo$NEW_MSI <- clinInfo$MSI
  clinInfo$NEW_MSI[which(clinInfo$NEW_MSI == "MSI-L")] <- "MSS"
  
  ### change other MSI-related info
  clinInfo$MSI_AGE_Status[intersect(which(clinInfo$MSI == "MSI-L"),
                                    which(clinInfo$`Diagnosis Age` < 50))] <- "MSS_Young"
  clinInfo$MSI_AGE_Status[intersect(which(clinInfo$MSI == "MSI-L"),
                                    which(clinInfo$`Diagnosis Age` >= 50))] <- "MSS_Old"
  clinInfo$MSI_RACE_Status1[intersect(which(clinInfo$MSI == "MSI-L"),
                                      which(clinInfo$`Race Category` == "BLACK OR AFRICAN AMERICAN"))] <- "MSS_AA"
  clinInfo$MSI_RACE_Status1[intersect(which(clinInfo$MSI == "MSI-L"),
                                      which(clinInfo$`Race Category` == "WHITE"))] <- "MSS_CC"
  clinInfo$MSI_RACE_Status2[intersect(which(clinInfo$MSI == "MSI-L"),
                                      which(clinInfo$Prediction_Filtered == "African"))] <- "MSS_AA"
  clinInfo$MSI_RACE_Status2[intersect(which(clinInfo$MSI == "MSI-L"),
                                      which(clinInfo$Prediction_Filtered == "Caucasian"))] <- "MSS_CC"
  
  ### msi status Indeterminate -> NA
  clinInfo$NEW_MSI[which(clinInfo$NEW_MSI == "Indeterminate")] <- NA
  
  ### add age categorical info (Old/Young)
  clinInfo$Age <- NA
  clinInfo$Age[which(clinInfo$`Diagnosis Age` >= 50)] <- "Old"
  clinInfo$Age[which(clinInfo$`Diagnosis Age` < 50)] <- "Young"
  
  ### add race categorical info
  clinInfo$Self_Reported_Race <- NA
  clinInfo$Self_Reported_Race[which(clinInfo$`Race Category` == "BLACK OR AFRICAN AMERICAN")] <- "AA"
  clinInfo$Self_Reported_Race[which(clinInfo$`Race Category` == "WHITE")] <- "CC"
  clinInfo$Predicted_Race <- NA
  clinInfo$Predicted_Race[which(clinInfo$Prediction_Filtered == "African")] <- "AA"
  clinInfo$Predicted_Race[which(clinInfo$Prediction_Filtered == "Caucasian")] <- "CC"
  
  ### load MAF file with read.table
  maf_table <- read.table(file = mafFilePath, header = TRUE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### make new maf_table only with the samples in the clinInfo (POLE-mutated samples removed)
  maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo)),]
  
  ### write out the new maf table
  write.table(maf_table, file = paste0(outputDir, "somatic_mutation_maf_tcga_coad_read_filtered.maf"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  ### a function to generate group-specific maf file
  specific_maf_file <- function(group) {
    ### split the group into pieces
    sp <- strsplit(group, "_", fixed = TRUE)[[1]]
    
    ### get maf for the specific group
    if((group == "Young") || (group == "Old")) {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$Age == group),])),]
    } else if((group == "AA") || (group == "CC")) {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$Self_Reported_Race == group),])),]
    } else if((group == "AA_Predicted") || (group == "CC_Predicted")) {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$Predicted_Race == sp[1]),])),]
    } else if((sp[length(sp)] == "Young") || (sp[length(sp)] == "Old")) {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$MSI_AGE_Status == group),])),]
    } else if((sp[length(sp)] == "AA") || (sp[length(sp)] == "CC")) {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$MSI_RACE_Status1 == group),])),]
    } else if(sp[length(sp)] == "Predicted") {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[which(clinInfo$MSI_RACE_Status2 == paste0(sp[1], "_", sp[2])),])),]
    } else if(sp[3] == "Predicted") {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[intersect(which(clinInfo$MSI_RACE_Status2 == paste0(sp[1], "_", sp[2])),
                                                                                                                      which(clinInfo$TUMOR_LOCATION == sp[4])),])),]
    } else if ((sp[2] == "AA") || (sp[2] == "CC")){
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[intersect(which(clinInfo$MSI_RACE_Status1 == paste0(sp[1], "_", sp[2])),
                                                                                                                      which(clinInfo$TUMOR_LOCATION == sp[3])),])),]
    } else {
      new_maf_table <- maf_table[which(substr(maf_table$Tumor_Sample_Barcode, 1, 15) %in% rownames(clinInfo[intersect(which(clinInfo$MSI_AGE_Status == paste0(sp[1], "_", sp[2])),
                                                                                                                      which(clinInfo$TUMOR_LOCATION == sp[3])),])),]
    }
    
    ### write out the specific maf table
    write.table(new_maf_table, file = paste0(outputDir, "somatic_mutation_maf_tcga_coad_read_filtered_", group, ".maf"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  ### run the specific_maf_file function for each group
  specific_maf_file("MSI-H_Young")
  specific_maf_file("MSI-H_Old")
  specific_maf_file("MSS_Young")
  specific_maf_file("MSS_Old")
  specific_maf_file("MSI-H_AA")
  specific_maf_file("MSI-H_CC")
  specific_maf_file("MSS_AA")
  specific_maf_file("MSS_CC")
  specific_maf_file("MSI-H_AA_Predicted")
  specific_maf_file("MSI-H_CC_Predicted")
  specific_maf_file("MSS_AA_Predicted")
  specific_maf_file("MSS_CC_Predicted")
  specific_maf_file("MSI-H_Young_Proximal")
  specific_maf_file("MSI-H_Old_Proximal")
  specific_maf_file("MSS_Young_Proximal")
  specific_maf_file("MSS_Old_Proximal")
  specific_maf_file("MSI-H_AA_Proximal")
  specific_maf_file("MSI-H_CC_Proximal")
  specific_maf_file("MSS_AA_Proximal")
  specific_maf_file("MSS_CC_Proximal")
  specific_maf_file("MSI-H_AA_Predicted_Proximal")
  specific_maf_file("MSI-H_CC_Predicted_Proximal")
  specific_maf_file("MSS_AA_Predicted_Proximal")
  specific_maf_file("MSS_CC_Predicted_Proximal")
  specific_maf_file("MSI-H_Young_Distal")
  specific_maf_file("MSI-H_Old_Distal")
  specific_maf_file("MSS_Young_Distal")
  specific_maf_file("MSS_Old_Distal")
  specific_maf_file("MSI-H_AA_Distal")
  specific_maf_file("MSI-H_CC_Distal")
  specific_maf_file("MSS_AA_Distal")
  specific_maf_file("MSS_CC_Distal")
  specific_maf_file("MSI-H_AA_Predicted_Distal")
  specific_maf_file("MSI-H_CC_Predicted_Distal")
  specific_maf_file("MSS_AA_Predicted_Distal")
  specific_maf_file("MSS_CC_Predicted_Distal")
  specific_maf_file("Young")
  specific_maf_file("Old")
  specific_maf_file("AA")
  specific_maf_file("CC")
  specific_maf_file("AA_Predicted")
  specific_maf_file("CC_Predicted")
  
  
  # Create variables to be used by the methods that map back and forth betweer symbols 
  # and entrez ids
  if(!exists("utils.env")){
    utils.env = new.env()
    utils.env$gene_symbols = keys(org.Hs.eg.db,  keytype="SYMBOL")
    utils.env$gene_aliases = keys(org.Hs.eg.db,  keytype="ALIAS")
  }
  
  # *****************************************************************************
  #
  # Map Gene Symbols to Gene IDs
  #
  # geneSymbol:	A single string or a vector of strings representing gene symbol(s).
  # 
  # Returns a vector of the same size as "geneSymbol" where the i-th entry is the 
  # gene ID (as an integer) corresponsing to the i-th gene symbol in "geneSymbol". 
  # The return vector entries are named using the gene symbols. For gene symbols mapped
  # to more than one entrez ids, only the first id is returned. If geneSymbol[i] is
  # a gene alias, the entrez Id of its corresponding official gene symbol is returned.
  #
  # ATTENTION: 
  # geneSymbol is first stripped of symbols which are not mapped to 
  # at least one entrez ID. If no gene symbols remain after this stripping, a NA 
  # value is returned.
  # *****************************************************************************
  geneSymbolToEntrezId <- function(geneSymbol) {
    mapped_symbol = geneSymbol[geneSymbol %in% utils.env$gene_symbols]
    mapped_alias = geneSymbol[geneSymbol %in% utils.env$gene_aliases]
    if (length(mapped_symbol) > 0 && length(mapped_alias) > 0)
      mapped_alias = setdiff(mapped_alias, mapped_symbol)
    ms = ma = matrix(nrow=0, ncol=2)
    if (length(mapped_symbol) > 0)
      ms = as.matrix(select(org.Hs.eg.db, keys=mapped_symbol, keytype="SYMBOL", columns=c("ENTREZID")))
    if (length(mapped_alias) > 0)
      ma = as.matrix(select(org.Hs.eg.db, keys=mapped_alias, keytype="ALIAS", columns=c("ENTREZID")))
    tmp = rbind(ms, ma)
    if (nrow(tmp) == 0)
      return(NA)		
    
    # Remove duplicate mappings, if any
    t = which(duplicated(tmp[,1]))
    if (length(t) > 0)
      tmp = tmp[-t, ]
    res = as.integer(tmp[,2])
    names(res) = tmp[, 1]
    return(res)	
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
  
  ### load total MAF file
  maf <- read.maf(maf = paste0(outputDir, "somatic_mutation_maf_tcga_coad_read_filtered.maf"))
  
  ### create a result directory
  resultDir <- paste0(outputDir, "Total/")
  dir.create(resultDir)
  
  ### create a summary plot
  png(paste0(resultDir, "Summary_Plot_Total.png"), width = 2000, height = 1000, res = 200)
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  ### get driver genes
  driver_genes <- oncodrive(maf = maf, AACol = "HGVSp_Short", minMut = 5)
  
  ### write out the driver gene result
  write.xlsx2(driver_genes, file = paste0(resultDir, "Driver_Genes_Table_Total.xlsx"),
              sheetName = "Driver_Genes", row.names = FALSE)
  
  ### if there are significant genes, draw a plot and perform pathway analysis
  if(length(which(driver_genes$fdr < driver_gene_fdr_cutoff)) > 0) {
    ### driver gene plot
    png(paste0(resultDir, "Driver_Genes_Plot_Total.png"), width = 1200, height = 1000, res = 130)
    plotOncodrive(res = driver_genes, fdrCutOff = driver_gene_fdr_cutoff, useFraction = TRUE)
    dev.off()
    
    ### pathway analysis
    pathways <- pathwayAnalysis_CP(geneList = geneSymbolToEntrezId(as.character(driver_genes$Hugo_Symbol)),
                                   org = "human", database = "GO", imgPrint = TRUE,
                                   title = "Driver_Genes_Pathways_Total",
                                   displayNum = 50, dir = resultDir)
    if(nrow(pathways) > 0) {
      write.xlsx2(pathways, file = paste0(resultDir, "Driver_Genes_Pathways_Total.xlsx"),
                  sheetName = "Driver_Genes_Pathways", row.names = FALSE)
    }
  }
  
  ### Total TMB
  tmb_df <- maf@variant.classification.summary
  tmb_df <- tmb_df[,-1]
  tmb_df <- tmb_df / 40
  rownames(tmb_df) <- maf@variant.classification.summary$Tumor_Sample_Barcode
  tmb_df$MSI_Info <- clinInfo[substr(rownames(tmb_df), 1, 15), "NEW_MSI"]
  tmb_df$Age_Info <- clinInfo[substr(rownames(tmb_df), 1, 15), "MSI_AGE_Status"]
  tmb_df$Race_Info <- clinInfo[substr(rownames(tmb_df), 1, 15), "MSI_RACE_Status2"]
  tmb_df$Survival <- clinInfo[substr(rownames(tmb_df), 1, 15), "Overall Survival (Months)"]
  
  ### write out the TMB table
  write.xlsx2(data.frame(Sample_Name=rownames(tmb_df), tmb_df, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(resultDir, "Tumor_Mutation_Burden_Total.xlsx"),
              sheetName = "Tumor_Mutation_Burden", row.names = FALSE)
  
  ### Overall TMB summary variable
  overall_tmb <- NULL
  overall_tmb_mis <- NULL
  
  ### a function to compare TMB for given comparison
  ### for each group also perform summary_plot, driver_gene, and pathway analyses
  mutation_analysis_for_comparison <- function(group1, group2) {
    
    ### load MAF files for each group
    group1_maf <- read.maf(maf = paste0(outputDir, "somatic_mutation_maf_tcga_coad_read_filtered_", group1, ".maf"))
    group2_maf <- read.maf(maf = paste0(outputDir, "somatic_mutation_maf_tcga_coad_read_filtered_", group2, ".maf"))
    
    ### summary plot for each group
    png(paste0(outputDir, "Summary/Summary_Plot_", group1, ".png"), width = 2000, height = 1000, res = 200)
    plotmafSummary(maf = group1_maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    png(paste0(outputDir, "Summary/Summary_Plot_", group2, ".png"), width = 2000, height = 1000, res = 200)
    plotmafSummary(maf = group2_maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    ### run driver gene and pathway analyses for each group
    ### group1 - run the below if there are at least 3 samples in the group
    if(nrow(group1_maf@variants.per.sample) > 2) {
      dir.create(paste0(outputDir, "Driver_Genes/", group1))
      driver_genes <- oncodrive(maf = group1_maf, AACol = "HGVSp_Short", minMut = 5)
      
      ### write out the driver gene result
      write.xlsx2(driver_genes, file = paste0(outputDir, "Driver_Genes/", group1, "/Driver_Genes_Table_", group1, ".xlsx"),
                  sheetName = "Driver_Genes", row.names = FALSE)
      
      if(length(which(driver_genes$fdr < driver_gene_fdr_cutoff)) > 0) {
        ### driver gene plot
        png(paste0(outputDir, "Driver_Genes/", group1, "/Driver_Genes_Plot_", group1, ".png"), width = 1200, height = 1000, res = 130)
        plotOncodrive(res = driver_genes, fdrCutOff = driver_gene_fdr_cutoff, useFraction = TRUE)
        dev.off()
        
        ### pathway analysis
        pathways <- pathwayAnalysis_CP(geneList = geneSymbolToEntrezId(as.character(driver_genes$Hugo_Symbol)),
                                       org = "human", database = "GO", imgPrint = TRUE,
                                       title = paste0("Driver_Genes_Pathways_", group1),
                                       displayNum = 50, dir = paste0(outputDir, "Driver_Genes/", group1, "/"))
        if(nrow(pathways) > 0) {
          write.xlsx2(pathways, file = paste0(outputDir, "Driver_Genes/", group1, "/Driver_Genes_Pathways_", group1, ".xlsx"),
                      sheetName = "Driver_Genes_Pathways", row.names = FALSE)
        }
      }
    }
    ### group2 - run the below if there are at least 3 samples in the group
    if(nrow(group2_maf@variants.per.sample) > 2) {
      dir.create(paste0(outputDir, "Driver_Genes/", group2))
      driver_genes <- oncodrive(maf = group2_maf, AACol = "HGVSp_Short", minMut = 5)
      
      ### write out the driver gene result
      write.xlsx2(driver_genes, file = paste0(outputDir, "Driver_Genes/", group2, "/Driver_Genes_Table_", group2, ".xlsx"),
                  sheetName = "Driver_Genes", row.names = FALSE)
      
      if(length(which(driver_genes$fdr < driver_gene_fdr_cutoff)) > 0) {
        ### driver gene plot
        png(paste0(outputDir, "Driver_Genes/", group2, "/Driver_Genes_Plot_", group2, ".png"), width = 1200, height = 1000, res = 130)
        plotOncodrive(res = driver_genes, fdrCutOff = driver_gene_fdr_cutoff, useFraction = TRUE)
        dev.off()
        
        ### pathway analysis
        pathways <- pathwayAnalysis_CP(geneList = geneSymbolToEntrezId(as.character(driver_genes$Hugo_Symbol)),
                                       org = "human", database = "GO", imgPrint = TRUE,
                                       title = paste0("Driver_Genes_Pathways_", group2),
                                       displayNum = 50, dir = paste0(outputDir, "Driver_Genes/", group2, "/"))
        if(nrow(pathways) > 0) {
          write.xlsx2(pathways, file = paste0(outputDir, "Driver_Genes/", group2, "/Driver_Genes_Pathways_", group2, ".xlsx"),
                      sheetName = "Driver_Genes_Pathways", row.names = FALSE)
        }
      }
    }
    
    ### TMB
    dir.create(paste0(outputDir, "TMB/", group1, "_vs_", group2))
    ### combine mutation count of each group
    df1 <- group1_maf@variant.classification.summary
    df1 <- df1[,-1]
    df1 <- df1 / 40
    df1$Sample_Group <- group1
    rownames(df1) <- group1_maf@variant.classification.summary$Tumor_Sample_Barcode
    df2 <- group2_maf@variant.classification.summary
    df2 <- df2[,-1]
    df2 <- df2 / 40
    df2$Sample_Group <- group2
    rownames(df2) <- group2_maf@variant.classification.summary$Tumor_Sample_Barcode
    mdf <- merge(df1, df2, all = TRUE, sort = FALSE)
    mdf[is.na(mdf)] <- 0
    rownames(mdf) <- c(rownames(df1), rownames(df2))
    if(grepl("Old", group1) || grepl("Young", group2)) {
      mdf$Additional_Info <- clinInfo[substr(rownames(mdf), 1, 15), "MSI_RACE_Status2"]
    } else {
      mdf$Additional_Info <- clinInfo[substr(rownames(mdf), 1, 15), "MSI_AGE_Status"]
    }
    mdf$Survival <- clinInfo[substr(rownames(mdf), 1, 15), "Overall Survival (Months)"]
    
    ### overall tmb
    overall_tmb <- c(get("overall_tmb"), mean(df1$total), mean(df2$total))
    overall_tmb_mis <- c(get("overall_tmb_mis"), mean(df1$Missense_Mutation), mean(df2$Missense_Mutation))
    names(overall_tmb)[(length(overall_tmb)-1):length(overall_tmb)] <- c(group1, group2)
    names(overall_tmb_mis)[(length(overall_tmb_mis)-1):length(overall_tmb_mis)] <- c(group1, group2)
    assign("overall_tmb", overall_tmb, envir = globalenv())
    assign("overall_tmb_mis", overall_tmb_mis, envir = globalenv())
    
    ### draw a beeswarm plot for TMB
    ### file name
    fName <- paste0("Beeswarm_Plot_TMB_", group1, "_vs_", group2)
    ### draw a plot with the data frame
    ggplot(mdf, aes(x=Sample_Group, y=total)) +
      ggtitle(fName) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes(color=Additional_Info), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Tumor Mutation Burden",
           subtitle = paste(group1, "Mean:", signif(mean(df1$total), 4), ",", group2, "Mean:", signif(mean(df2$total), 4))) +
      theme(legend.position="right", plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
    ggsave(filename = paste0(outputDir, "TMB/", group1, "_vs_", group2, "/", fName, ".png"), width = 15, height = 12)
    
    ### draw a beeswarm plot for TMB (missense only)
    ### file name
    fName <- paste0("Beeswarm_Plot_TMB_Missense_", group1, "_vs_", group2)
    ### draw a plot with the data frame
    ggplot(mdf, aes(x=Sample_Group, y=Missense_Mutation)) +
      ggtitle(fName) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes(color=Additional_Info), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Tumor Mutation Burden (Missense Only)",
           subtitle = paste(group1, "Mean:", signif(mean(df1$Missense_Mutation), 4), ",", group2, "Mean:", signif(mean(df2$Missense_Mutation), 4))) +
      theme(legend.position="right", plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
    ggsave(filename = paste0(outputDir, "TMB/", group1, "_vs_", group2, "/", fName, ".png"), width = 15, height = 12)
    
    ### write out the TMB table
    write.xlsx2(data.frame(Sample_Name=rownames(mdf), mdf, stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir, "TMB/", group1, "_vs_", group2, "/", "TMB_Table_", group1, "_vs_", group2, ".xlsx"),
                sheetName = "Tumor_Mutation_Burden", row.names = FALSE)
    
    ### correlation between survival and TMB
    if((length(which(mdf$Sample_Group == group1)) > 2) && (length(which(mdf$Sample_Group == group2)) > 2)) {
      fName <- paste0("Correlation_Plot_TMB_Survival_", group1, "_vs_", group2)
      colors <- c("skyblue", "pink")
      names(colors) <- c(group1, group2)
      png(paste0(outputDir, "TMB/", group1, "_vs_", group2, "/", fName, ".png"), width = 1500, height = 1200, res = 130)
      par(oma = c(0,0,3,0))
      plot(mdf$total, mdf$Survival, pch = 19,
           col = colors[mdf$Sample_Group],
           main = paste0(names(colors)[1], " P.Cor = ", round(cor(mdf$total[which(mdf$Sample_Group == names(colors)[1])],
                                                                  mdf$Survival[which(mdf$Sample_Group == names(colors)[1])], use = "pairwise.complete.obs"), 5),
                         ", p-value = ", signif(cor.test(mdf$total[which(mdf$Sample_Group == names(colors)[1])],
                                                         mdf$Survival[which(mdf$Sample_Group == names(colors)[1])])$p.value, 5),
                         "\n", names(colors)[2], " P.Cor = ", round(cor(mdf$total[which(mdf$Sample_Group == names(colors)[2])],
                                                                        mdf$Survival[which(mdf$Sample_Group == names(colors)[2])], use = "pairwise.complete.obs"), 5),
                         ", p-value = ", signif(cor.test(mdf$total[which(mdf$Sample_Group == names(colors)[2])],
                                                         mdf$Survival[which(mdf$Sample_Group == names(colors)[2])])$p.value, 5)),
           xlab = "Tumor Mutation Burden",
           ylab = "Overall Survival (Months)")
      abline(lm(mdf$Survival[which(mdf$Sample_Group == names(colors)[1])]~mdf$total[which(mdf$Sample_Group == names(colors)[1])]), col=colors[1], lwd=2)
      abline(lm(mdf$Survival[which(mdf$Sample_Group == names(colors)[2])]~mdf$total[which(mdf$Sample_Group == names(colors)[2])]), col=colors[2], lwd=2)
      legend("topright", legend = names(colors),
             col = colors, pch = 19,
             title = "Sample Groups", cex = 0.8)
      mtext(fName, outer = TRUE, cex = 2)
      dev.off()
    }
    
  }
  
  ### run mutation_analysis_for_comparison for each comparison
  dir.create(paste0(outputDir, "Summary"))
  dir.create(paste0(outputDir, "TMB"))
  dir.create(paste0(outputDir, "Driver_Genes"))
  mutation_analysis_for_comparison("MSI-H_Young", "MSI-H_Old")
  mutation_analysis_for_comparison("MSS_Young", "MSS_Old")
  mutation_analysis_for_comparison("MSI-H_AA", "MSI-H_CC")
  mutation_analysis_for_comparison("MSS_AA", "MSS_CC")
  mutation_analysis_for_comparison("MSI-H_AA_Predicted", "MSI-H_CC_Predicted")
  mutation_analysis_for_comparison("MSS_AA_Predicted", "MSS_CC_Predicted")
  mutation_analysis_for_comparison("MSI-H_Young_Proximal", "MSI-H_Young_Distal")
  mutation_analysis_for_comparison("MSS_Young_Proximal", "MSS_Young_Distal")
  mutation_analysis_for_comparison("MSI-H_Old_Proximal", "MSI-H_Old_Distal")
  mutation_analysis_for_comparison("MSS_Old_Proximal", "MSS_Old_Distal")
  mutation_analysis_for_comparison("MSI-H_AA_Proximal", "MSI-H_AA_Distal")
  mutation_analysis_for_comparison("MSS_AA_Proximal", "MSS_AA_Distal")
  mutation_analysis_for_comparison("MSI-H_CC_Proximal", "MSI-H_CC_Distal")
  mutation_analysis_for_comparison("MSS_CC_Proximal", "MSS_CC_Distal")
  mutation_analysis_for_comparison("MSI-H_AA_Predicted_Proximal", "MSI-H_AA_Predicted_Distal")
  mutation_analysis_for_comparison("MSS_AA_Predicted_Proximal", "MSS_AA_Predicted_Distal")
  mutation_analysis_for_comparison("MSI-H_CC_Predicted_Proximal", "MSI-H_CC_Predicted_Distal")
  mutation_analysis_for_comparison("MSS_CC_Predicted_Proximal", "MSS_CC_Predicted_Distal")
  mutation_analysis_for_comparison("Young", "Old")
  mutation_analysis_for_comparison("AA", "CC")
  mutation_analysis_for_comparison("AA_Predicted", "CC_Predicted")
  
  ### write out the overall TMB result table
  write.xlsx2(data.frame(Sample_Group=names(overall_tmb), TMB=overall_tmb, TMB_Missense=overall_tmb_mis,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(resultDir, "Tumor_Mutation_Burden_Total_Per_Group.xlsx"),
              sheetName = "Tumor_Mutation_Burden", row.names = FALSE)
  
}
