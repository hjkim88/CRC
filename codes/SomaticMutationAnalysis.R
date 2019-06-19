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

sm_analysis <- function(vcfRDAfilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/somatic_mutation/somatic_mutation_vcfs_tcga_coad_read.rda",
                        clinInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated.txt",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/somatic_mutation/") {
  
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
  if(!require(survminer, quietly = TRUE)) {
    install.packages("survminer")
    require(survminer, quietly = TRUE)
  }
  if(!require(survival, quietly = TRUE)) {
    install.packages("survival")
    require(survival, quietly = TRUE)
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
                        MSI=clinInfo[names(sm_num),"MSI"],
                        Age=clinInfo[names(sm_num),"MSI_AGE_Status"],
                        Self_Reported_Race=clinInfo[names(sm_num),"MSI_RACE_Status1"],
                        Predicted_Race=clinInfo[names(sm_num),"MSI_RACE_Status2"],
                        Survival=clinInfo[names(sm_num),"Overall Survival (Months)"],
                        Status=clinInfo[names(sm_num),"Patient Vital Status"],
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
  
  
  ### if we divide the samples based on the number of somatic mutation,
  ### does it affect the survival?
  
  ### it seems like MSI-H samples have larger number of somatic mutations than MSS samples
  ### there's no big difference between young vs old in general
  ### no difference in MSI-H_AA vs MSI-H_CC
  ### MSS_CC has more somatic mutations than MSS_AA
  
  ### order based on the number of somatic mutation
  plot_df <- plot_df[order(plot_df$The_Number_Of_Somatic_Mutations),]
  
  ### get the first and the third quantile and make a new data frame
  q <- quantile(plot_df$The_Number_Of_Somatic_Mutations)
  plot_df$Group <- c(rep("Small", length(which(plot_df$The_Number_Of_Somatic_Mutations < q[2]))),
                     rep("Small-Intermediate", length(intersect(which(plot_df$The_Number_Of_Somatic_Mutations >= q[2]),
                                                                        which(plot_df$The_Number_Of_Somatic_Mutations < q[3])))),
                     rep("Intermediate-Large", length(intersect(which(plot_df$The_Number_Of_Somatic_Mutations >= q[3]),
                                                                        which(plot_df$The_Number_Of_Somatic_Mutations < q[4])))),
                     rep("Large", length(which(plot_df$The_Number_Of_Somatic_Mutations >= q[4]))))
  
  ### survival analysis
  plot_df$Survival <- as.numeric(plot_df$Survival)
  plot_df$Status[plot_df$Status == "Alive"] <- 0
  plot_df$Status[plot_df$Status == "Dead"] <- 1
  plot_df$Status <- as.numeric(plot_df$Status)
  fit <- survfit(Surv(Survival, Status) ~ Group, data = plot_df)
  
  ### print survival plot
  fName <- "Survival_Plot_Somatic_Mutation"
  png(paste0(outputDir, fName, ".png"), width = 1500, height = 1200, res = 130)
  ggsurvplot(
    fit,
    data = plot_df,
    title = "Survival Differences Based On Somatic Mutation",
    legend.labs = levels(as.factor(plot_df$Group)),
    risk.table = TRUE,
    tables.col = "strata",
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Time in Months",
    break.time.by = round(max(plot_df$Survival, na.rm = TRUE)/5),
    ggtheme = theme_classic(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.25,
    conf.int.style = "ribbon"
  )
  dev.off()
  
  ### correlation plot
  ### the number of somatic mutations - survival outcomes
  
  ### MSI-H vs MSS
  new_plot_df <- plot_df[union(which(plot_df$MSI == "MSI-H"),
                               which(plot_df$MSI == "MSS")),]
  colors <- c("skyblue", "pink")
  names(colors) <- c("MSI-H", "MSS")
  
  ### print out
  fName <- "Correlation_Plot_MSI-H_vs_MSS"
  png(paste0(outputDir, fName, ".png"), width = 1500, height = 1200, res = 130)
  par(oma = c(0,0,3,0))
  plot(new_plot_df$The_Number_Of_Somatic_Mutations, new_plot_df$Survival, pch = 19,
       col = colors[new_plot_df$MSI],
       main = paste0(names(colors)[1], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[1])],
                                                             new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[1])], use = "pairwise.complete.obs"), 5),
                    ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[1])],
                                                    new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[1])])$p.value, 5),
                    "\n", names(colors)[2], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[2])],
                                                                   new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[2])], use = "pairwise.complete.obs"), 5),
                    ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[2])],
                                                    new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[2])])$p.value, 5)),
       xlab = "The Number of Somatic Mutations",
       ylab = "Overall Survival (Months)")
  abline(lm(new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[1])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[1])]), col=colors[1], lwd=2)
  abline(lm(new_plot_df$Survival[which(new_plot_df$MSI == names(colors)[2])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$MSI == names(colors)[2])]), col=colors[2], lwd=2)
  legend("topright", legend = names(colors),
         col = colors, pch = 19,
         title = "Sample Groups", cex = 0.8)
  mtext(fName, outer = TRUE, cex = 2)
  dev.off()
  
  ### MSS_AA vs MSS_CC
  new_plot_df <- plot_df[union(which(plot_df$Self_Reported_Race == "MSS_AA"),
                               which(plot_df$Self_Reported_Race == "MSS_CC")),]
  colors <- c("skyblue", "pink")
  names(colors) <- c("MSS_AA", "MSS_CC")
  
  ### print out
  fName <- "Correlation_Plot_MSS_AA_vs_CC"
  png(paste0(outputDir, fName, ".png"), width = 1500, height = 1200, res = 130)
  par(oma = c(0,0,3,0))
  plot(new_plot_df$The_Number_Of_Somatic_Mutations, new_plot_df$Survival, pch = 19,
       col = colors[new_plot_df$Self_Reported_Race],
       main = paste0(names(colors)[1], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[1])],
                                                              new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[1])], use = "pairwise.complete.obs"), 5),
                     ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[1])],
                                                     new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[1])])$p.value, 5),
                     "\n", names(colors)[2], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[2])],
                                                                    new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[2])], use = "pairwise.complete.obs"), 5),
                     ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[2])],
                                                     new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[2])])$p.value, 5)),
       xlab = "The Number of Somatic Mutations",
       ylab = "Overall Survival (Months)")
  abline(lm(new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[1])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[1])]), col=colors[1], lwd=2)
  abline(lm(new_plot_df$Survival[which(new_plot_df$Self_Reported_Race == names(colors)[2])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Self_Reported_Race == names(colors)[2])]), col=colors[2], lwd=2)
  legend("topright", legend = names(colors),
         col = colors, pch = 19,
         title = "Sample Groups", cex = 0.8)
  mtext(fName, outer = TRUE, cex = 2)
  dev.off()
  
  ### MSS_AA vs MSS_CC with the predicted race info
  new_plot_df <- plot_df[union(which(plot_df$Predicted_Race == "MSS_AA"),
                               which(plot_df$Predicted_Race == "MSS_CC")),]
  colors <- c("skyblue", "pink")
  names(colors) <- c("MSS_AA", "MSS_CC")
  
  ### print out
  fName <- "Correlation_Plot_Predicted_MSS_AA_vs_CC"
  png(paste0(outputDir, fName, ".png"), width = 1500, height = 1200, res = 130)
  par(oma = c(0,0,3,0))
  plot(new_plot_df$The_Number_Of_Somatic_Mutations, new_plot_df$Survival, pch = 19,
       col = colors[new_plot_df$Predicted_Race],
       main = paste0(names(colors)[1], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[1])],
                                                              new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[1])], use = "pairwise.complete.obs"), 5),
                     ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[1])],
                                                     new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[1])])$p.value, 5),
                     "\n", names(colors)[2], " P.Cor = ", round(cor(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[2])],
                                                                    new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[2])], use = "pairwise.complete.obs"), 5),
                     ", p-value = ", signif(cor.test(new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[2])],
                                                     new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[2])])$p.value, 5)),
       xlab = "The Number of Somatic Mutations",
       ylab = "Overall Survival (Months)")
  abline(lm(new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[1])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[1])]), col=colors[1], lwd=2)
  abline(lm(new_plot_df$Survival[which(new_plot_df$Predicted_Race == names(colors)[2])]~new_plot_df$The_Number_Of_Somatic_Mutations[which(new_plot_df$Predicted_Race == names(colors)[2])]), col=colors[2], lwd=2)
  legend("topright", legend = names(colors),
         col = colors, pch = 19,
         title = "Sample Groups", cex = 0.8)
  mtext(fName, outer = TRUE, cex = 2)
  dev.off()
  
}
