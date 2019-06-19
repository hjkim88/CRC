###
#   File name : SomaticMutationAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Jun 18, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load VCF RDA file and compare the number of somatic mutations in each comparison
#
#   Instruction
#               1. Source("SomaticMutationAnalysis.R")
#               2. Run the function "sm_analysis" - specify the necessary input file paths and output directory
#               3. The analysis results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_SomaticMutationAnalysis.R/SomaticMutationAnalysis.R")
#               > sm_analysis(vcfRDAfilePath="./data/somatic_mutation/somatic_mutation_vcfs_tcga_coad_read.rda",
#                             clinInfoPath="./data/coadread_tcga_clinical_data_updated.txt",
#                             outputDir="./results/somatic_mutation/")
###

sm_analysis <- function(vcfRDAfilePath="./data/somatic_mutation/somatic_mutation_vcfs_tcga_coad_read.rda",
                        clinInfoPath="./data/coadread_tcga_clinical_data_updated.txt",
                        outputDir="./results/somatic_mutation/") {
  
  ### load libraries
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  ### load data
  load(vcfRDAfilePath)
  
  ### load clinical info
  clinInfo <- read.table(clinInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinInfo) <- clinInfo$`Patient ID`
  
  ### get the number of somatic mutation of the samples
  sm_num <- sapply(vcf, nrow, USE.NAMES = TRUE)
  
  ### create a data frame for the beeswarm plot
  plot_df <- data.frame(The_Number_Of_Somatic_Mutations=sm_num,
                        Age=clinInfo[names(sm_num),"MSI_AGE_Status"],
                        Self_Reported_Race=clinInfo[names(sm_num),"MSI_RACE_Status1"],
                        Predicted_Race=clinInfo[names(sm_num),"MSI_RACE_Status2"],
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### Age
  ### file name
  fName <- "Beeswarm_Plot_Somatic_Mutation_Age"
  ### draw a plot with the data frame
  ggplot(plot_df, aes(x=Age, y=The_Number_Of_Somatic_Mutations)) +
    ggtitle(fName) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Self_Reported_Race), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "The Number Of Somatic Mutations") +
    theme(legend.position="right", plot.title=element_text(hjust = 0.5))
  ggsave(filename = paste0(outputDir, fName, ".png"), width = 20, height = 12)
  
  ### Self-reported race
  ### file name
  fName <- "Beeswarm_Plot_Somatic_Mutation_Race"
  ### draw a plot with the data frame
  ggplot(plot_df, aes(x=Self_Reported_Race, y=The_Number_Of_Somatic_Mutations)) +
    ggtitle(fName) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Age), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "The Number Of Somatic Mutations") +
    theme(legend.position="right", plot.title=element_text(hjust = 0.5))
  ggsave(filename = paste0(outputDir, fName, ".png"), width = 20, height = 12)
  
  ### Predicted race
  ### file name
  fName <- "Beeswarm_Plot_Somatic_Mutation_Predicted_Race"
  ### draw a plot with the data frame
  ggplot(plot_df, aes(x=Predicted_Race, y=The_Number_Of_Somatic_Mutations)) +
    ggtitle(fName) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Age), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "The Number Of Somatic Mutations") +
    theme(legend.position="right", plot.title=element_text(hjust = 0.5))
  ggsave(filename = paste0(outputDir, fName, ".png"), width = 20, height = 12)
  
  
  ### survival analysis
  
  
}
