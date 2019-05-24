###
#   File name : DemographicAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Apr 15, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a degraphic table and plots with the TCGA-COAD-READ data
#
#   Instruction
#               1. Source("DemographicAnalysis.R")
#               2. Run the function "demoAnalysis" - specify clinical info file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DemographicAnalysis.R/DemographicAnalysis.R")
#               > demoAnalysis(clinInfoPath_640 = "./data/coadread_tcga_clinical_data.tsv",
#                              msiInfoPath = "./data/nationwidechildrens.org_auxiliary_coad_read.txt",
#                              predictedRaceInfoPath="./results/genotyping/Race_Prediction_Result.txt",
#                              outputDir="./results/demographic/")
###

demoAnalysis <- function(clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                         msiInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_auxiliary_coad_read.txt",
                         predictedRaceInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/genotyping/Race_Prediction_Result.txt",
                         outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/demographic/") {
  
  ### load necessary libraries
  if(!require(tidyverse, quietly = TRUE)) {
    install.packages("tidyverse")
    require(tidyverse, quietly = TRUE)
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
  
  ### load the data
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
  
  ### add group info (MSI/MSS - African American/Caucasian)
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
  
  ### make a data frame for plots
  plot_df <- clinicalInfo_640[,
                           c("Sample ID",
                             "Diagnosis Age",
                             "American Joint Committee on Cancer Tumor Stage Code",
                             "Cancer Type Detailed",
                             "Death from Initial Pathologic Diagnosis Date",
                             "Patient Height",
                             "Mutation Count",
                             "Overall Survival (Months)",
                             "Patient Primary Tumor Site",
                             "Race Category",
                             "Sample Initial Weight",
                             "Sex",
                             "Patient Vital Status",
                             "Patient Weight",
                             "MSI",
                             "MSI_AGE_Status",
                             "MSI_RACE_Status1",
                             "MSI_RACE_Status2")]
  
  ### change the colnames of the plot_df
  colnames(plot_df) <- c("Sample_ID",
                      "Age",
                      "Tumor_Stage",
                      "Cancer_Type",
                      "Death_From_Diagnosis",
                      "Height",
                      "Mutation_Count",
                      "Survival",
                      "Tumor_Site",
                      "Race",
                      "Initial_Weight",
                      "Gender",
                      "Status",
                      "Weight",
                      "MSI",
                      "MSI_AGE_Status",
                      "MSI_RACE_Status1",
                      "MSI_RACE_Status2")
  
  
  ### MSI/AGE
  
  ### remove rows with MSI_AGE_Status == NA
  plot_df1 <- plot_df
  if(length(which(is.na(plot_df1$MSI_AGE_Status))) > 0) {
    plot_df1 <- plot_df1[-which(is.na(plot_df1$MSI_AGE_Status)),]
  }
  
  ### make NA to "NA"
  plot_df1[is.na(plot_df1)] <- "NA"
  
  ### set parameters for the circular plot (stage)
  group <- unique(plot_df1$MSI_AGE_Status)
  separator_size <- 3
  tumor_stage_list <- unique(plot_df1$Tumor_Stage)[order(unique(plot_df1$Tumor_Stage))]
  
  ### make a data frame for the circular plot (stage)
  cdf <- data.frame(Stage=rep(NA, (length(tumor_stage_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Stage[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list))] <- tumor_stage_list
    cdf$Group[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list)+separator_size)] <- group[i]
    for(j in 1:length(tumor_stage_list)) {
      cdf$Number[(i-1)*(length(tumor_stage_list)+separator_size)+j] <- length(intersect(which(plot_df1$MSI_AGE_Status == group[i]), which(plot_df1$Tumor_Stage == tumor_stage_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Stage, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (stage)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g1 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                   top = "Stage Differences in Various MSI/Age Status")
  ggsave(file = paste0(outputDir, "demo_msi_age_stage.png"), g1, width = 20, height = 12)
  
  ### set parameters for the circular plot (gender)
  separator_size <- 2
  gender_list <- unique(plot_df1$Gender)[order(unique(plot_df1$Gender))]
  
  ### make a data frame for the circular plot (gender)
  cdf <- data.frame(Gender=rep(NA, (length(gender_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Gender[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list))] <- gender_list
    cdf$Group[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list)+separator_size)] <- group[i]
    for(j in 1:length(gender_list)) {
      cdf$Number[(i-1)*(length(gender_list)+separator_size)+j] <- length(intersect(which(plot_df1$MSI_AGE_Status == group[i]), which(plot_df1$Gender == gender_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Gender, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (gender)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g2 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                    top = "Gender Differences in Various MSI/Age Status")
  ggsave(file = paste0(outputDir, "demo_msi_age_gender.png"), g2, width = 20, height = 12)
  
  ### survival plot
  plot_df1$Survival <- as.numeric(plot_df1$Survival)
  plot_df1$Status[plot_df1$Status == "Alive"] <- 0
  plot_df1$Status[plot_df1$Status == "Dead"] <- 1
  plot_df1$Status <- as.numeric(plot_df1$Status)
  fit <- survfit(Surv(Survival, Status) ~ MSI_AGE_Status, data = plot_df1)
  
  p6 <- ggsurvplot(
    fit,
    data = plot_df1,
    title = "Survival Differences in Various MSI/Age Status",
    legend.labs = levels(as.factor(plot_df1$MSI_AGE_Status)),
    risk.table = TRUE,
    tables.col = "strata",
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Time in Months",
    break.time.by = round(max(plot_df1$Survival, na.rm = TRUE)/5),
    ggtheme = theme_classic(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.25,
    conf.int.style = "ribbon"
  )
  p6$table <- p6$table + theme(legend.position = "none")
  p6$ncensor.plot <- p6$ncensor.plot + theme(legend.position = "none")
  
  ### beeswarm plot
  p7 <- ggplot(plot_df1, aes(x=MSI_AGE_Status, y=Survival)) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=MSI_AGE_Status), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "Survival (Months)") +
    theme(legend.position="top")
  
  ### arrange the plots and print out
  g3 <- arrangeGrob(p6$plot, p6$table, p6$ncensor.plot, p7,
                    layout_matrix = rbind(c(1, 4), c(1, 4), c(1, 4), c(2, 4), c(3, 4)),
                    top = "Survival Differences in Various MSI/Age Status")
  ggsave(file = paste0(outputDir, "demo_msi_age_survival.png"), g3, width = 20, height = 12)
  
  
  ### MSI/RACE - SELF-REPORTED INFO
  
  ### remove rows with MSI_RACE_Status1 == NA
  plot_df2 <- plot_df
  if(length(which(is.na(plot_df2$MSI_RACE_Status1))) > 0) {
    plot_df2 <- plot_df2[-which(is.na(plot_df2$MSI_RACE_Status1)),]
  }
  
  ### make NA to "NA"
  plot_df2[is.na(plot_df2)] <- "NA"
  
  ### set parameters for the circular plot (stage)
  group <- unique(plot_df2$MSI_RACE_Status1)
  separator_size <- 3
  tumor_stage_list <- unique(plot_df2$Tumor_Stage)[order(unique(plot_df2$Tumor_Stage))]
  
  ### make a data frame for the circular plot (stage)
  cdf <- data.frame(Stage=rep(NA, (length(tumor_stage_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Stage[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list))] <- tumor_stage_list
    cdf$Group[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list)+separator_size)] <- group[i]
    for(j in 1:length(tumor_stage_list)) {
      cdf$Number[(i-1)*(length(tumor_stage_list)+separator_size)+j] <- length(intersect(which(plot_df2$MSI_RACE_Status1 == group[i]), which(plot_df2$Tumor_Stage == tumor_stage_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Stage, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (stage)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g1 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                    top = "Stage Differences in Various MSI/Race Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_stage.png"), g1, width = 20, height = 12)
  
  ### set parameters for the circular plot (gender)
  separator_size <- 2
  gender_list <- unique(plot_df2$Gender)[order(unique(plot_df2$Gender))]
  
  ### make a data frame for the circular plot (gender)
  cdf <- data.frame(Gender=rep(NA, (length(gender_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Gender[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list))] <- gender_list
    cdf$Group[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list)+separator_size)] <- group[i]
    for(j in 1:length(gender_list)) {
      cdf$Number[(i-1)*(length(gender_list)+separator_size)+j] <- length(intersect(which(plot_df2$MSI_RACE_Status1 == group[i]), which(plot_df2$Gender == gender_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Gender, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (gender)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g2 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                    top = "Gender Differences in Various MSI/RACE Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_gender.png"), g2, width = 20, height = 12)
  
  ### survival plot
  plot_df2$Survival <- as.numeric(plot_df2$Survival)
  plot_df2$Status[plot_df2$Status == "Alive"] <- 0
  plot_df2$Status[plot_df2$Status == "Dead"] <- 1
  plot_df2$Status <- as.numeric(plot_df2$Status)
  fit <- survfit(Surv(Survival, Status) ~ MSI_RACE_Status1, data = plot_df2)
  
  p6 <- ggsurvplot(
    fit,
    data = plot_df2,
    title = "Survival Differences in Various MSI/Race Status",
    legend.labs = levels(as.factor(plot_df2$MSI_RACE_Status1)),
    risk.table = TRUE,
    tables.col = "strata",
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Time in Months",
    break.time.by = round(max(plot_df2$Survival, na.rm = TRUE)/5),
    ggtheme = theme_classic(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.25,
    conf.int.style = "ribbon"
  )
  p6$table <- p6$table + theme(legend.position = "none")
  p6$ncensor.plot <- p6$ncensor.plot + theme(legend.position = "none")
  
  ### beeswarm plot
  p7 <- ggplot(plot_df2, aes(x=MSI_RACE_Status1, y=Survival)) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=MSI_RACE_Status1), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "Survival (Months)") +
    theme(legend.position="top")
  
  ### arrange the plots and print out
  g3 <- arrangeGrob(p6$plot, p6$table, p6$ncensor.plot, p7,
                    layout_matrix = rbind(c(1, 4), c(1, 4), c(1, 4), c(2, 4), c(3, 4)),
                    top = "Survival Differences in Various MSI/Race Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_survival.png"), g3, width = 20, height = 12)
  
  
  ### MSI/RACE - PREDICTED INFO
  
  ### remove rows with MSI_RACE_Status2 == NA
  plot_df3 <- plot_df
  if(length(which(is.na(plot_df3$MSI_RACE_Status2))) > 0) {
    plot_df3 <- plot_df3[-which(is.na(plot_df3$MSI_RACE_Status2)),]
  }
  
  ### make NA to "NA"
  plot_df3[is.na(plot_df3)] <- "NA"
  
  ### set parameters for the circular plot (stage)
  group <- unique(plot_df3$MSI_RACE_Status2)
  separator_size <- 3
  tumor_stage_list <- unique(plot_df3$Tumor_Stage)[order(unique(plot_df3$Tumor_Stage))]
  
  ### make a data frame for the circular plot (stage)
  cdf <- data.frame(Stage=rep(NA, (length(tumor_stage_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Stage[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list))] <- tumor_stage_list
    cdf$Group[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list)+separator_size)] <- group[i]
    for(j in 1:length(tumor_stage_list)) {
      cdf$Number[(i-1)*(length(tumor_stage_list)+separator_size)+j] <- length(intersect(which(plot_df3$MSI_RACE_Status2 == group[i]), which(plot_df3$Tumor_Stage == tumor_stage_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Stage, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (stage)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Stage)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g1 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                    top = "Stage Differences in Various MSI/Race Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_predicted_stage.png"), g1, width = 20, height = 12)
  
  ### set parameters for the circular plot (gender)
  separator_size <- 2
  gender_list <- unique(plot_df3$Gender)[order(unique(plot_df3$Gender))]
  
  ### make a data frame for the circular plot (gender)
  cdf <- data.frame(Gender=rep(NA, (length(gender_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Gender[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list))] <- gender_list
    cdf$Group[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list)+separator_size)] <- group[i]
    for(j in 1:length(gender_list)) {
      cdf$Number[(i-1)*(length(gender_list)+separator_size)+j] <- length(intersect(which(plot_df3$MSI_RACE_Status2 == group[i]), which(plot_df3$Gender == gender_list[j])))
    }
  }
  cdf$ID <- 1:nrow(cdf)
  
  ### get the name and the y position of each label
  label_data <- cdf
  number_of_bar <- nrow(label_data)
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle <- 90 - 360 * (label_data$ID - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  ### prepare a data frame for base lines
  base_data <- cdf %>% 
    group_by(Group) %>% 
    summarize(start=min(ID), end=max(ID) - separator_size) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data <- base_data[order(base_data$start),]
  
  ### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$start <- grid_data$start - 0.5
  grid_data$end <- grid_data$end + separator_size - 0.5
  
  ### make the circular plot
  quant <- round(quantile(1:max(cdf$Number, na.rm = TRUE)))
  p1 <- ggplot(cdf, aes(x = as.factor(ID), y = Number, fill = Group)) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    # add a val=25/50/75/100 lines. Do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = start, y = quant[2], xend = end, yend = quant[2]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[3], xend = end, yend = quant[3]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[4], xend = end, yend = quant[4]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    geom_segment(data = grid_data, aes(x = start, y = quant[5], xend = end, yend = quant[5]), colour = "grey", alpha = 1, size = 0.3 , inherit.aes = FALSE) +
    # add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(base_data$end[1]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[2]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[3]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    annotate("text", x = rep(base_data$end[4]+separator_size,4), y = quant[2:5], label = as.character(quant[2:5]) , color = "grey", size = 4 , angle = 0, fontface = "bold", hjust = 1) +
    geom_bar(aes(x = as.factor(ID), y = Number, fill = Group), stat = "identity", alpha = 0.5) +
    ylim(-1.1*quant[5],1.1*quant[5]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = ID, y = Number+1, label = Gender, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 4, angle = label_data$angle, inherit.aes = FALSE ) +
    # add base line information
    geom_segment(data = base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -10, label = group), hjust = c(1,1,0,0), colour = "black", alpha = 0.8, size = 4, fontface = "bold", inherit.aes = FALSE)
  
  ### make pie charts (gender)
  cdf <- data.frame(cdf, pcnt = NA, stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(group)) {
    df <- cdf[cdf$Group == group[i],]
    cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
  }
  if(length(which(cdf$pcnt < 3)) > 0) {
    cdf$pcnt[which(cdf$pcnt < 3)] <- ""
    cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
  }
  # for loop failed because of a bug of geom_text (it uses the fixed text input - not dynamic)
  # just manually perform it 4 times
  p2 <- ggplot(data = cdf[cdf$Group == group[1],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[1]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[1])
  p3 <- ggplot(data = cdf[cdf$Group == group[2],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[2]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[2])
  p4 <- ggplot(data = cdf[cdf$Group == group[3],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[3]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[3])
  p5 <- ggplot(data = cdf[cdf$Group == group[4],], aes(x = "", y = Number, fill = Gender)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, color = "black")) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    geom_text(aes(label = cdf$pcnt[cdf$Group == group[4]]), position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, title = group[4])
  
  ### arrange the plots and print out
  g2 <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                    top = "Gender Differences in Various MSI/RACE Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_predicted_gender.png"), g2, width = 20, height = 12)
  
  ### survival plot
  plot_df3$Survival <- as.numeric(plot_df3$Survival)
  plot_df3$Status[plot_df3$Status == "Alive"] <- 0
  plot_df3$Status[plot_df3$Status == "Dead"] <- 1
  plot_df3$Status <- as.numeric(plot_df3$Status)
  fit <- survfit(Surv(Survival, Status) ~ MSI_RACE_Status2, data = plot_df3)
  
  p6 <- ggsurvplot(
    fit,
    data = plot_df3,
    title = "Survival Differences in Various MSI/Race Status",
    legend.labs = levels(as.factor(plot_df3$MSI_RACE_Status2)),
    risk.table = TRUE,
    tables.col = "strata",
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Time in Months",
    break.time.by = round(max(plot_df3$Survival, na.rm = TRUE)/5),
    ggtheme = theme_classic(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.25,
    conf.int.style = "ribbon"
  )
  p6$table <- p6$table + theme(legend.position = "none")
  p6$ncensor.plot <- p6$ncensor.plot + theme(legend.position = "none")
  
  ### beeswarm plot
  p7 <- ggplot(plot_df3, aes(x=MSI_RACE_Status2, y=Survival)) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=MSI_RACE_Status2), na.rm = TRUE) +
    stat_compare_means() +
    labs(x = "", y = "Survival (Months)") +
    theme(legend.position="top")
  
  ### arrange the plots and print out
  g3 <- arrangeGrob(p6$plot, p6$table, p6$ncensor.plot, p7,
                    layout_matrix = rbind(c(1, 4), c(1, 4), c(1, 4), c(2, 4), c(3, 4)),
                    top = "Survival Differences in Various MSI/Race Status")
  ggsave(file = paste0(outputDir, "demo_msi_race_predicted_survival.png"), g3, width = 20, height = 12)
  
}
