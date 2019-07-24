###
#   File name : CombinedAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Jun 5, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Analyze using both Methylation and expression levels of genes
#               It is an analysis using both methylation data and RNA-Seq
#
#   Instruction
#               1. Source("CombinedAnalysis.R")
#               2. Run the function "combAnalysis" - specify the input files (Differential) of methylation and RNA-Seq and output directory
#               3. Various results will be generated in the output path
#
#   Example
#               > source("The_directory_of_CombinedAnalysis.R/CombinedAnalysis.R")
#               > combAnalysis(methylPath="./results/methylation/differential_methylation/",
#                              methylRDAPath="./data/methylation/preprocessed/norm_beta_tcga_coad_read.rda",
#                              rnaseqPath="./results/rnaseq/preprocessed_raw_counts3/",
#                              rnaseqRDAPath="./data/rnaseq/raw_count_tcga_coad_read.rda",
#                              geneRIFPath="../Periodontitis/data/generifs_basic.txt",
#                              clinInfoPath="./data/coadread_tcga_clinical_data_updated2.txt",
#                              fdrThreshold=0.05,
#                              logFCThreshold=0,
#                              outputDir="./results/combined/")
###

combAnalysis <- function(methylPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/methylation/differential_methylation/",
                         methylRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/norm_beta_tcga_coad_read.rda",
                         rnaseqPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/rnaseq/",
                         rnaseqRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/raw_count_tcga_coad_read.rda",
                         geneRIFPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/generifs_basic.txt",
                         clinInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt",
                         fdrThreshold=0.05,
                         logFCThreshold=0,
                         outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/combined/") {
  
  ### load libraries
  options(java.parameters = "-Xmx8000m")
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(ggrepel)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ggrepel")
    library(ggrepel)
  }
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  
  ### get comparison list from the methylation result path
  comparisons <- list.dirs(methylPath, full.names = FALSE, recursive = FALSE)
  
  ### load geneRIF
  geneRIF <- read.table(file = geneRIFPath, header = FALSE, sep = "\t", check.names = FALSE,
                        quote = "", stringsAsFactors = FALSE)
  
  ### load clinical info
  clinInfo <- read.table(clinInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinInfo) <- clinInfo$`Sample ID`
  
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
  
  ### remove POLE-muated samples
  clinicalInfo_640 <- clinicalInfo_640[-which(clinicalInfo_640$POLE_MUTANT == TRUE),]
  
  ### load the data
  load(methylRDAPath)
  load(rnaseqRDAPath)
  
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
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
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
  
  ### get normalized data
  gexpLev <- normalizeRNASEQwithVST(raw_count_tcga_coad_read, filtering = FALSE)
  methylLev <- normB$beta
  
  ### should keep same sample order
  common_samples <- intersect(colnames(gexpLev), colnames(methylLev))
  gexpLev <- gexpLev[,common_samples]
  methylLev <- methylLev[,common_samples]
  
  ### remove unneccesarry objects
  rm(raw_count_tcga_coad_read)
  rm(normB)
  gc()
  
  ### for each comparison, compare methylation result and rnaseq result
  for(comp in comparisons) {
    
    ### print progress
    writeLines(paste0(comp))
    
    ### create a directory for the comparison
    dir.create(path = paste0(outputDir, comp))
    
    ### load DMP result
    dmps <- read.xlsx2(file = paste0(methylPath, comp, "/DMP/DMPs_", comp, ".xlsx"),
                       sheetIndex = 1, row.names = 1, colClasses = NA,
                       stringsAsFactors = FALSE, check.names = FALSE)
    
    ### load DMR result
    dmrs <- read.xlsx2(file = paste0(methylPath, comp, "/DMR/DMRs_", comp, ".xlsx"),
                       sheetIndex = 1, row.names = 1, colClasses = NA,
                       stringsAsFactors = FALSE, check.names = FALSE)
    
    ### load DE result
    degs <- read.xlsx2(file = paste0(rnaseqPath, comp, "/deresult_", comp, ".xlsx"),
                       sheetIndex = 1, row.names = 1, colClasses = NA,
                       stringsAsFactors = FALSE, check.names = FALSE)
    
    ### using results with the cutoff
    ### filter the data with the given cutoff
    dmps <- dmps[intersect(which(dmps$adj.P.Val < fdrThreshold), which(abs(dmps$logFC) > logFCThreshold)),]
    dmrs <- dmrs[which(dmrs$Stouffer < fdrThreshold),]
    degs <- degs[intersect(which(degs$padj < fdrThreshold), which(abs(degs$log2FoldChange) > logFCThreshold)),]
    
    ### annotate gene info to the DMRs
    dmrs$genes <- NA
    for(i in 1:nrow(dmrs)) {
      genes <- unique(sapply(strsplit(dmrs$overlapping.promoters[i], split = ", ", fixed = TRUE)[[1]],
                             function(x) strsplit(x, split = "-", fixed = TRUE)[[1]][1], USE.NAMES = FALSE))
      dmrs$genes[i] <- paste(genes, collapse = ";")
    }
    
    ### create new columns and put NAs
    degs$CpGs <- NA
    degs$DMRs <- NA
    
    ### get DE gene associated CpG sites & DMRs, put the info on the data
    for(i in 1:nrow(degs)) {
      ### find indicies of differentially methylated CpGs associated with the DE genes
      idx <- which(dmps$gene == rownames(degs)[i])
      
      ### add the CpG info to the gene
      if(length(idx) > 0) {
        degs$CpGs[i] <- paste(rownames(dmps)[idx], collapse = ";")
      }
      
      ### find indicies of DMRs associated with the DE genes
      if(nrow(dmrs) > 0) {
        idx <- NULL
        for(j in 1:nrow(dmrs)) {
          genes <- strsplit(dmrs$genes[j], split = ";", fixed = TRUE)[[1]]
          if(length(genes) > 0) {
            for(k in 1:length(genes)) {
              if(rownames(degs)[i] == genes[k]) {
                idx <- c(idx, j)
                break
              }
            }
          }
        }
      }
      
      ### add the DMR info to the gene
      if(length(idx) > 0) {
        degs$DMRs[i] <- paste(rownames(dmrs)[idx], collapse = ";")
      }
    }
    
    ### run only if there are (more than 2) DE genes - associated CpGs
    assoc_idx <- which(!is.na(degs$CpGs))
    if(length(assoc_idx) > 2) {
    
      ### create new column and put NAs
      degs$CpG <- NA
      degs$MethylFDR <- NA
      
      ### get one CpG per gene that has the lowest FDR
      for(i in assoc_idx) {
        temp <- strsplit(degs$CpGs[i], ";")[[1]]
        min_value <- Inf
        for(j in 1:length(temp)) {
          if(min_value > dmps[temp[j],"adj.P.Val"]) {
            min_value <- dmps[temp[j],"adj.P.Val"]
            degs$CpG[i] <- temp[j]
          }
        }
        degs$MethylFDR[i] <- min_value
      }
      
      ### make a correlation data frame for a plot
      cor_data <- data.frame(Expression=degs$padj[assoc_idx], Methylation=degs$MethylFDR[assoc_idx],
                             Gene_Name=rownames(degs)[assoc_idx],
                             stringsAsFactors = FALSE, check.names = FALSE)
      
      ### mark the significant genes
      higher_thresh <- 1
      idx <- 1:nrow(cor_data)
      while(length(idx) > 20) {
        higher_thresh <- higher_thresh / 10
        cor_data$Label <- ""
        idx <- intersect(which(cor_data$Expression < higher_thresh), which(cor_data$Methylation < higher_thresh))
        cor_data$Label[idx] <- cor_data$Gene_Name[idx]
      }
      if(higher_thresh == 1) {
        higher_thresh <- fdrThreshold
        cor_data$Label <- cor_data$Gene_Name
      }
      
      ### log10 transformation
      cor_data[,c(1,2)] <- -log10(cor_data[,c(1,2)])
      
      fName <- paste0("Plot_-log10(FDR)_DMP_", comp, ".png")
      ggplot(data = cor_data, aes(x=Expression, y=Methylation)) +
        geom_point(color = "black", size = 1) +
        geom_label_repel(aes(Expression, Methylation, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
        labs(title=substr(fName, 1, nchar(fName)-4),
             subtitle=sprintf("P.Cor = %s, p-value = %s",
                              round(cor(cor_data[,1], cor_data[,2], use = "pairwise.complete.obs"), 5),
                              signif(cor.test(cor_data[,1], cor_data[,2])$p.value, 5))) +
        xlab("Expression") +
        ylab("Methylation") +
        # geom_smooth(method = lm, color="gray", se=FALSE) +
        geom_vline(aes(xintercept = -log10(higher_thresh), linetype = paste0("-log10(", higher_thresh, ")")), color = "red") +
        geom_hline(aes(yintercept = -log10(higher_thresh), linetype = paste0("-log10(", higher_thresh, ")")), color = "red") +
        scale_linetype_manual(name = "FDR Threshold", values = c(2),
                              guide = guide_legend(override.aes = list(color = c("red")))) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, comp, "/", fName), width = 10, height = 10)
      
      ### make a table
      result_table <- matrix(NA, 1, 20)
      colnames(result_table) <- c(paste0("DEG_", c("Gene_Name", colnames(degs)[c(1:6, 9)])), paste0("DMP_", colnames(dmps)[c(1:6, 10:13, 15:16)]))
      for(idx in assoc_idx) {
        cpgs <- strsplit(degs$CpGs[idx], ";", fixed = TRUE)[[1]]
        for(cpg in cpgs) {
          result_table <- rbind(result_table, c(rownames(degs)[idx], degs[idx,1:6], cpg, dmps[cpg, c(1:6, 10:13, 15:16)]))
        }
      }
      result_table <- result_table[-1,]
      
      ### write out the table result
      write.xlsx2(result_table, file = paste0(outputDir, comp, "/", "DEG_DMP_combined_", comp, ".xlsx"),
                  sheetName = comp, row.names = FALSE)
      
      ### GeneRIF analysis
      
      ### get Entrez IDs
      entIDs <- geneSymbolToEntrezId(unique(cor_data$Label[which(cor_data$Label != "")]))
      
      ### extract geneRIF for the specific genes
      add_info <- geneRIF[which(geneRIF$V2 %in% entIDs),c(2,5)]
      colnames(add_info) <- c("Gene_Name", "Gene_RIF")
      
      ### change Entrez IDs to gene symbols
      name_map <- names(entIDs)
      names(name_map) <- entIDs
      add_info[,1] <- name_map[as.character(add_info[,1])]
      
      ### remove duplicates
      rIdx <- duplicated.data.frame(add_info)
      add_info <- add_info[!rIdx,]
      
      ### write out the result
      write.xlsx2(add_info, file = paste0(outputDir, comp, "/", "GeneRIF_DMP_", comp, ".xlsx"), row.names = FALSE)
      
      ### pathway analysis
      
      ### get pathways
      pathways <- pathwayAnalysis_CP(geneList = geneSymbolToEntrezId(unique(cor_data$Gene_Name)), org = "human", database = "GO", imgPrint = TRUE,
                                     title = paste0("Top_50_DEMG_Pathways_DMP_", comp),
                                     displayNum = 50, dir = paste0(outputDir, comp, "/"))
      write.xlsx2(pathways, file = paste0(outputDir, comp, "/", "go_DEMG-associated_Pathways_DMP_", comp, ".xlsx"),
                  sheetName = "DEMG-associated_Pathways", row.names = FALSE)
      
      ### level correlation
      
      ### collect genes that are DE & DM
      important_genes <- unique(cor_data$Label[which(cor_data$Label != "")])
       
      ### organize samples for the plot
      if(grepl("Proximal", comp)) {
        ### get target class
        compare <- strsplit(comp, "_", fixed = TRUE)[[1]]
        group <- paste0(compare[1], "_", compare[2])
        
        ### get sample names for each class
        if(compare[3] == "predicted") {
          compare <- c(paste0(group, "_", compare[4]), paste0(group, "_", compare[6]))
          samps <- vector("list", length = 3)
          names(samps) <- c(compare, "Total")
          
          samps[[compare[1]]] <- intersect(common_samples,
                                           clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_RACE_Status2 == group),
                                                                          which(clinInfo$TUMOR_LOCATION == "Proximal"))])
          samps[[compare[2]]] <- intersect(common_samples,
                                           clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_RACE_Status2 == group),
                                                                          which(clinInfo$TUMOR_LOCATION == "Distal"))])
        } else {
          compare <- c(paste0(group, "_", compare[3]), paste0(group, "_", compare[5]))
          samps <- vector("list", length = 3)
          names(samps) <- c(compare, "Total")
          
          if(grepl("Young", group) || grepl("Old", group)) {
            samps[[compare[1]]] <- intersect(common_samples,
                                             clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_AGE_Status == group),
                                                                            which(clinInfo$TUMOR_LOCATION == "Proximal"))])
            samps[[compare[2]]] <- intersect(common_samples,
                                             clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_AGE_Status == group),
                                                                            which(clinInfo$TUMOR_LOCATION == "Distal"))])
          } else {
            samps[[compare[1]]] <- intersect(common_samples,
                                             clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_RACE_Status1 == group),
                                                                            which(clinInfo$TUMOR_LOCATION == "Proximal"))])
            samps[[compare[2]]] <- intersect(common_samples,
                                             clinInfo$`Sample ID`[intersect(which(clinInfo$MSI_RACE_Status1 == group),
                                                                            which(clinInfo$TUMOR_LOCATION == "Distal"))])
          }
        }
        samps[["Total"]] <- c(samps[[compare[1]]], samps[[compare[2]]])
      } else {
        ### get class info of the comparison
        compare <- strsplit(comp, "_", fixed = TRUE)[[1]]
        compare <- c(paste0(compare[1], "_", compare[2]), paste0(compare[1], "_", compare[4]))
        
        ### get sample names for each class
        samps <- vector("list", length = 3)
        names(samps) <- c(compare, "Total")
        for(one in compare) {
          if(endsWith(comp, "predicted")) {
            samps[[one]] <- intersect(common_samples, clinInfo$`Sample ID`[which(clinInfo$MSI_RACE_Status2 == one)])
          } else if(endsWith(comp, "CC")) {
            samps[[one]] <- intersect(common_samples, clinInfo$`Sample ID`[which(clinInfo$MSI_RACE_Status1 == one)])
          } else {
            samps[[one]] <- intersect(common_samples, clinInfo$`Sample ID`[which(clinInfo$MSI_AGE_Status == one)])
          }
        }
        samps[["Total"]] <- c(samps[[compare[1]]], samps[[compare[2]]])
      }
      
      ### print level correlation plots
      if((length(samps[[1]]) > 2) && (length(samps[[2]]) > 2)) {
        png(paste0(outputDir, comp, "/", "DEMG_Level_Correlations_DMP_", comp, ".png"), width = 2000, height = 1200, res = 120)
        
        ### decide grid based on the number of the plots
        if(length(important_genes) == 2) {
          par(mfrow=c(1,2), oma = c(0,0,3,0))
        } else if(length(important_genes) <= 4) {
          par(mfrow=c(2,2), oma = c(0,0,3,0))
        } else if(length(important_genes) <= 6) {
          par(mfrow=c(2,3), oma = c(0,0,3,0))
        } else if(length(important_genes) <= 9) {
          par(mfrow=c(3,3), oma = c(0,0,3,0))
        } else if(length(important_genes) <= 12) {
          par(mfrow=c(3,4), oma = c(0,0,3,0))
        } else {
          par(mfrow=c(4,5), oma = c(0,0,3,0))
        }
        
        ### color for the sample based on their classes
        ### based on - samps[["Total"]] <- c(samps[[compare[1]]], samps[[compare[2]]])
        colors <- c(rep("blue", length(samps[[compare[1]]])), rep("green", length(samps[[compare[2]]])))
        
        ### for each gene draw a plot
        for(gene in important_genes) {
          cpg <- degs[gene,"CpG"]
          plot(as.numeric(gexpLev[gene,samps[["Total"]]]), methylLev[cpg,samps[["Total"]]], pch = 19,
               col = colors,
               main = paste0(compare[1], " P.Cor = ", round(cor(as.numeric(gexpLev[gene,samps[[compare[1]]]]), methylLev[cpg,samps[[compare[1]]]], use = "pairwise.complete.obs"), 5),
                             ", p-value = ", signif(cor.test(as.numeric(gexpLev[gene,samps[[compare[1]]]]), methylLev[cpg,samps[[compare[1]]]])$p.value, 5),
                             "\n", compare[2], " P.Cor = ", round(cor(as.numeric(gexpLev[gene,samps[[compare[2]]]]), methylLev[cpg,samps[[compare[2]]]], use = "pairwise.complete.obs"), 5),
                             ", p-value = ", signif(cor.test(as.numeric(gexpLev[gene,samps[[compare[2]]]]), methylLev[cpg,samps[[compare[2]]]])$p.value, 5)),
               xlab = paste("Normalized Expression of", gene),
               ylab = paste("Normalized Beta value of", cpg))
          abline(lm(methylLev[cpg,samps[[compare[1]]]]~as.numeric(gexpLev[gene,samps[[compare[1]]]])), col="blue", lwd=2)
          abline(lm(methylLev[cpg,samps[[compare[2]]]]~as.numeric(gexpLev[gene,samps[[compare[2]]]])), col="green", lwd=2)
          legend("topright", legend = compare,
                 col = c("blue", "green"), pch = 19,
                 title = "Sample Groups", cex = 0.8)
        }
        mtext(paste("Level Correlations", comp), outer = TRUE, cex = 2)
        dev.off()
      }
      
    }
    
    ### run only if there are (more than 1) DE genes - associated DMRs
    assoc_idx <- which(!is.na(degs$DMRs))
    if(length(assoc_idx) > 2) {
      
      ### create new column and put NAs
      degs$DMR <- NA
      degs$MethylStouffer <- NA
      
      ### calculate min FDR of DMRs
      for(i in assoc_idx) {
        temp <- strsplit(degs$DMRs[i], ";")[[1]]
        min_value <- Inf
        for(j in 1:length(temp)) {
          if(min_value > dmrs[temp[j],"Stouffer"]) {
            min_value <- dmrs[temp[j],"Stouffer"]
            degs$DMR[i] <- temp[j]
          }
        }
        degs$MethylStouffer[i] <- min_value
      }
      
      ### make a correlation data frame for a plot
      cor_data <- data.frame(Expression=degs$padj[assoc_idx], Methylation=degs$MethylStouffer[assoc_idx],
                             Gene_Name=rownames(degs)[assoc_idx],
                             stringsAsFactors = FALSE, check.names = FALSE)
      
      ### mark the significant genes
      higher_thresh <- 1
      idx <- 1:nrow(cor_data)
      while(length(idx) > 20) {
        higher_thresh <- higher_thresh / 10
        cor_data$Label <- ""
        idx <- intersect(which(cor_data$Expression < higher_thresh), which(cor_data$Methylation < higher_thresh))
        cor_data$Label[idx] <- cor_data$Gene_Name[idx]
      }
      if(higher_thresh == 1) {
        higher_thresh <- fdrThreshold
        cor_data$Label <- cor_data$Gene_Name
      }
      
      ### log10 transformation
      cor_data[,c(1,2)] <- -log10(cor_data[,c(1,2)])
      
      fName <- paste0("Plot_-log10(FDR)_DMR_", comp, ".png")
      ggplot(data = cor_data, aes(x=Expression, y=Methylation)) +
        geom_point(color = "black", size = 1) +
        geom_label_repel(aes(Expression, Methylation, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
        labs(title=substr(fName, 1, nchar(fName)-4),
             subtitle=sprintf("P.Cor = %s, p-value = %s",
                              round(cor(cor_data[,1], cor_data[,2], use = "pairwise.complete.obs"), 5),
                              signif(cor.test(cor_data[,1], cor_data[,2])$p.value, 5))) +
        xlab("Expression") +
        ylab("Methylation") +
        # geom_smooth(method = lm, color="gray", se=FALSE) +
        geom_vline(aes(xintercept = -log10(higher_thresh), linetype = paste0("-log10(", higher_thresh, ")")), color = "red") +
        geom_hline(aes(yintercept = -log10(higher_thresh), linetype = paste0("-log10(", higher_thresh, ")")), color = "red") +
        scale_linetype_manual(name = "FDR Threshold", values = c(2),
                              guide = guide_legend(override.aes = list(color = c("red")))) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, comp, "/", fName), width = 10, height = 10)
      
      ### make a table
      result_table <- matrix(NA, 1, 15)
      colnames(result_table) <- c(paste0("DEG_", c("Gene_Name", colnames(degs)[c(1:6, 11)])), paste0("DMR_", colnames(dmrs)[c(1:6, 8)]))
      for(idx in assoc_idx) {
        dmrss <- strsplit(degs$DMRs[idx], ";", fixed = TRUE)[[1]]
        for(dmr in dmrss) {
          result_table <- rbind(result_table, c(rownames(degs)[idx], degs[idx,1:6], dmr, dmrs[dmr, c(1:6, 8)]))
        }
      }
      result_table <- result_table[-1,]
      
      ### write out the table result
      write.xlsx2(result_table, file = paste0(outputDir, comp, "/", "DEG_DMR_combined_", comp, ".xlsx"),
                  sheetName = comp, row.names = FALSE)
      
      ### GeneRIF analysis
      
      ### get Entrez IDs
      entIDs <- geneSymbolToEntrezId(unique(cor_data$Label[which(cor_data$Label != "")]))
      
      ### extract geneRIF for the specific genes
      add_info <- geneRIF[which(geneRIF$V2 %in% entIDs),c(2,5)]
      colnames(add_info) <- c("Gene_Name", "Gene_RIF")
      
      ### change Entrez IDs to gene symbols
      name_map <- names(entIDs)
      names(name_map) <- entIDs
      add_info[,1] <- name_map[as.character(add_info[,1])]
      
      ### remove duplicates
      rIdx <- duplicated.data.frame(add_info)
      add_info <- add_info[!rIdx,]
      
      ### write out the result
      write.xlsx2(add_info, file = paste0(outputDir, comp, "/", "GeneRIF_DMR_", comp, ".xlsx"), row.names = FALSE)
      
      ### pathway analysis
      
      ### get pathways
      pathways <- pathwayAnalysis_CP(geneList = geneSymbolToEntrezId(unique(cor_data$Gene_Name)), org = "human", database = "GO", imgPrint = TRUE,
                                     title = paste0("Top_50_DEMG_Pathways_DMR_", comp),
                                     displayNum = 50, dir = paste0(outputDir, comp, "/"))
      write.xlsx2(pathways, file = paste0(outputDir, comp, "/", "go_DEMG-associated_Pathways_DMR_", comp, ".xlsx"),
                  sheetName = "DEMG-associated_Pathways", row.names = FALSE)
      
    }
    
  }
  
}
