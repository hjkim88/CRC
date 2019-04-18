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
#               > demoAnalysis(clinInfoPath_coad = "./data/nationwidechildrens.org_clinical_patient_coad.txt",
#                              clinInfoPath_read = "./data/nationwidechildrens.org_clinical_patient_read.txt",
#                              clinInfoPath_640 = "./data/coadread_tcga_clinical_data.tsv",
#                              outputDir="./results/demographic/")
###

demoAnalysis <- function(clinInfoPath_coad = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_clinical_patient_coad.txt",
                         clinInfoPath_read = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/nationwidechildrens.org_clinical_patient_read.txt",
                         clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
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
  clinicalInfo_COAD <- read.table(file = clinInfoPath_coad,
                                  header = TRUE, sep = "\t",stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_COAD <- clinicalInfo_COAD[-c(1,2),]
  clinicalInfo_READ <- read.table(file = clinInfoPath_read,
                                  header = TRUE, sep = "\t",stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_READ <- clinicalInfo_READ[-c(1,2),]
  clinicalInfo_640 <- read.table(file = clinInfoPath_640, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
  
  ### add msi status
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 MSI=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  msi_grp <- sample(nrow(clinicalInfo_640), round(nrow(clinicalInfo_640)/2))
  mss_grp <- setdiff(c(1:nrow(clinicalInfo_640)), msi_grp)
  clinicalInfo_640$MSI[msi_grp] <- "MSI"
  clinicalInfo_640$MSI[mss_grp] <- "MSS"
  
  ### add group info (MSI/MSS - Young/Old)
  clinicalInfo_640 <- data.frame(clinicalInfo_640,
                                 Group1=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  clinicalInfo_640$Group1[intersect(which(clinicalInfo_640$MSI == "MSI"),
                                    which(clinicalInfo_640$`Diagnosis Age` < 50))] <- "MSI_Young"
  clinicalInfo_640$Group1[intersect(which(clinicalInfo_640$MSI == "MSI"),
                                    which(clinicalInfo_640$`Diagnosis Age` >= 50))] <- "MSI_Old"
  clinicalInfo_640$Group1[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                    which(clinicalInfo_640$`Diagnosis Age` < 50))] <- "MSS_Young"
  clinicalInfo_640$Group1[intersect(which(clinicalInfo_640$MSI == "MSS"),
                                    which(clinicalInfo_640$`Diagnosis Age` >= 50))] <- "MSS_Old"
  
  
  ### make a data frame for beeswarm plot
  bsdf <- clinicalInfo_640[union(which(clinicalInfo_640$`Race Category` == "BLACK OR AFRICAN AMERICAN"),
                                 which(clinicalInfo_640$`Race Category` == "WHITE")),
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
                             "Group1")]
  
  ### change the colnames of the bsdf
  colnames(bsdf) <- c("Sample_ID",
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
                      "Group1")
  
  ### make NA to "NA"
  bsdf[is.na(bsdf)] <- "NA"
  
  ### draw a beeswarm plot
  ggplot(bsdf, aes(x=MSI, y=Age)) +
    ggtitle(paste("")) +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Tumor_Stage)) +
    stat_compare_means()
  ggsave(filename = paste0(""), width = 12, height = 10)
  
  ### set parameters for the circular plot (stage)
  group <- unique(bsdf$Group1)
  separator_size <- 3
  tumor_stage_list <- unique(bsdf$Tumor_Stage)[order(unique(bsdf$Tumor_Stage))]
  
  ### make a data frame for the circular plot (stage)
  cdf <- data.frame(Stage=rep(NA, (length(tumor_stage_list)+separator_size)*length(group)),
                    Group=NA, Number=NA, ID=NA)
  for(i in 1:length(group)) {
    cdf$Stage[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list))] <- tumor_stage_list
    cdf$Group[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list)+separator_size)] <- group[i]
    for(j in 1:length(tumor_stage_list)) {
      cdf$Number[(i-1)*(length(tumor_stage_list)+separator_size)+j] <- length(intersect(which(bsdf$Group1 == group[i]), which(bsdf$Tumor_Stage == tumor_stage_list[j])))
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
  grid_data$start <- grid_data$start
  grid_data$end <- grid_data$end + separator_size -1
  
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
    ylim(-100,120) +
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
    geom_segment(data = base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha = 1, size = 1, inherit.aes = FALSE) +
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
  g <- arrangeGrob(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 1, 5, 2), c(1, 1, 4, 3)),
                   top = "Stage Differences in Various MSI/Age Status")
  ggsave(file = paste0(outputDir, "demo_msi_age_stage.png"), g, width = 20, height = 10)
  
  ### survival plot
  bsdf$Survival <- as.numeric(bsdf$Survival)
  bsdf$Status[bsdf$Status == "Alive"] <- 0
  bsdf$Status[bsdf$Status == "Dead"] <- 1
  bsdf$Status <- as.numeric(bsdf$Status)
  fit <- survfit(Surv(Survival, Status) ~ Group1, data = bsdf)
  
  p6 <- ggsurvplot(
    fit,
    data = bsdf,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Time in Months",
    break.time.by = round(max(bsdf$Survival, na.rm = TRUE)/5),
    ggtheme = theme_classic(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = TRUE,
    ncensor.plot.height = 0.25,
    conf.int.style = "ribbon"
  )
  
  
  
  
  
  
}







