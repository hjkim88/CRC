###
#   File name : DemographicAnalysis2.R
#   Author    : Hyunjin Kim
#   Date      : July 11, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a degraphic table and plots with the TCGA-COAD-READ data
#
#   Instruction
#               1. Source("DemographicAnalysis2.R")
#               2. Run the function "demoAnalysis" - specify clinical info file paths and the output directory
#               3. The results will be generated under the output directory
#
#   * Now it uses updated clinical info other than original clinical info,
#     POLE-mutated samples are removed, and MSI-L samples are now treated as MSS.
#     colorectal location analysis has been also added.
#     The aboves are the differences from the DemographicAnalysis.R.
#
#   Example
#               > source("The_directory_of_DemographicAnalysis2.R/DemographicAnalysis2.R")
#               > demoAnalysis2(clinInfoPath_640 = "./data/coadread_tcga_clinical_data_updated2.txt",
#                               pie_chart_option = c("percentage", "count"),
#                               outputDir="./results/demographic/")
###

demoAnalysis2 <- function(clinInfoPath_640 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt",
                          pie_chart_option = c("percentage", "count"),
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
  rownames(clinicalInfo_640) <- clinicalInfo_640$`Sample ID`
  
  ### remove POLE-mutated samples
  clinicalInfo_640 <- clinicalInfo_640[-which(clinicalInfo_640$POLE_MUTANT == TRUE),]
  
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
                                "NEW_MSI",
                                "MSI_AGE_Status",
                                "MSI_RACE_Status1",
                                "MSI_RACE_Status2",
                                "TUMOR_LOCATION")]
  
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
                         "MSI_RACE_Status2",
                         "TUMOR_LOCATION")
  
  
  ### a function to perform demographic anlaysis with an input
  ### df = a data frame that contains all the sample info and will be used to generate plots
  ### colName = a column name of the df that will be the group info of the analysis
  demo_plots <- function(input_df, colName) {
    
    ### survfit() needs global variable
    assign("colName", colName, envir = globalenv())
    
    ### create a directory for the results
    dir.create(paste0(outputDir, colName))
    
    ### remove rows with grp == NA
    new_df <- input_df
    if(length(which(is.na(new_df[,colName]))) > 0) {
      new_df <- new_df[-which(is.na(new_df[,colName])),]
    }
    
    ### make NA to "NA"
    new_df[is.na(new_df)] <- "NA"
    
    ### set parameters for the circular plot (stage)
    group <- unique(new_df[,colName])
    separator_size <- 3
    tumor_stage_list <- unique(new_df$Tumor_Stage)[order(unique(new_df$Tumor_Stage))]
    
    ### make a data frame for the circular plot (stage)
    cdf <- data.frame(Stage=rep(NA, (length(tumor_stage_list)+separator_size)*length(group)),
                      Group=NA, Number=NA, ID=NA)
    for(i in 1:length(group)) {
      cdf$Stage[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list))] <- tumor_stage_list
      cdf$Group[((i-1)*(length(tumor_stage_list)+separator_size)+1):((i-1)*(length(tumor_stage_list)+separator_size)+length(tumor_stage_list)+separator_size)] <- group[i]
      for(j in 1:length(tumor_stage_list)) {
        cdf$Number[(i-1)*(length(tumor_stage_list)+separator_size)+j] <- length(intersect(which(new_df[,colName] == group[i]), which(new_df$Tumor_Stage == tumor_stage_list[j])))
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
    if(pie_chart_option[1] == "percentage") {
      for(i in 1:length(group)) {
        df <- cdf[cdf$Group == group[i],]
        cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
      }
      cdf$pcnt[which(cdf$pcnt < 3)] <- ""
      cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
    } else if(pie_chart_option[1] == "count") {
      cdf$pcnt <- cdf$Number
      cdf$pcnt[which(cdf$pcnt == 0)] <- ""
    } else {
      stop("The pie_chart_option should be either \"percentage\" or \"count\".")
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
                      top = paste0("Stage Differences in Various ", colName))
    ggsave(file = paste0(outputDir, colName, "/demo_", colName, "_stage.png"), g1, width = 20, height = 12)
    
    ### set parameters for the circular plot (gender)
    separator_size <- 2
    gender_list <- unique(new_df$Gender)[order(unique(new_df$Gender))]
    
    ### make a data frame for the circular plot (gender)
    cdf <- data.frame(Gender=rep(NA, (length(gender_list)+separator_size)*length(group)),
                      Group=NA, Number=NA, ID=NA)
    for(i in 1:length(group)) {
      cdf$Gender[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list))] <- gender_list
      cdf$Group[((i-1)*(length(gender_list)+separator_size)+1):((i-1)*(length(gender_list)+separator_size)+length(gender_list)+separator_size)] <- group[i]
      for(j in 1:length(gender_list)) {
        cdf$Number[(i-1)*(length(gender_list)+separator_size)+j] <- length(intersect(which(new_df[,colName] == group[i]), which(new_df$Gender == gender_list[j])))
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
    if(pie_chart_option[1] == "percentage") {
      for(i in 1:length(group)) {
        df <- cdf[cdf$Group == group[i],]
        cdf$pcnt[cdf$Group == group[i]] <- round((df$Number / sum(df$Number, na.rm = TRUE))*100)
      }
      cdf$pcnt[which(cdf$pcnt < 3)] <- ""
      cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)] <- paste0(cdf$pcnt[which(as.numeric(cdf$pcnt) >= 3)], "%")
    } else if(pie_chart_option[1] == "count") {
      cdf$pcnt <- cdf$Number
      cdf$pcnt[which(cdf$pcnt == 0)] <- ""
    } else {
      stop("The pie_chart_option should be either \"percentage\" or \"count\".")
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
                      top = paste0("Gender Differences in Various ", colName))
    ggsave(file = paste0(outputDir, colName, "/demo_", colName, "_gender.png"), g2, width = 20, height = 12)
    
    ### survival plot
    new_df$Survival <- as.numeric(new_df$Survival)
    new_df$Status[new_df$Status == "Alive"] <- 0
    new_df$Status[new_df$Status == "Dead"] <- 1
    new_df$Status <- as.numeric(new_df$Status)
    fit <- survfit(as.formula(paste("Surv(Survival, Status)", "~", colName)), data = new_df)
    
    ### additional
    # new_df <- new_df[union(which(new_df$MSI_RACE_Status2 == "MSS_AA"),
    #                        which(new_df$MSI_RACE_Status2 == "MSS_CC")),]
    # fit <- survfit(as.formula(paste("Surv(Survival, Status)", "~", colName)), data = new_df)
    
    if(length(which(fit$n < 2)) > 0) {
      p6 <- ggsurvplot(
        fit,
        data = new_df,
        title = paste0("Survival Differences in Various ", colName),
        legend.labs = levels(as.factor(new_df[,colName])),
        risk.table = TRUE,
        tables.col = "strata",
        pval = TRUE,
        conf.int = TRUE,
        xlab = "Time in Months",
        break.time.by = round(max(new_df$Survival, na.rm = TRUE)/5),
        ggtheme = theme_classic(),
        risk.table.y.text.col = TRUE,
        risk.table.height = 0.25,
        risk.table.y.text = FALSE,
        ncensor.plot = TRUE,
        ncensor.plot.height = 0.25,
        conf.int.style = "step"
      )
    } else {
      p6 <- ggsurvplot(
        fit,
        data = new_df,
        title = paste0("Survival Differences in Various ", colName),
        legend.labs = levels(as.factor(new_df[,colName])),
        risk.table = TRUE,
        tables.col = "strata",
        pval = TRUE,
        conf.int = TRUE,
        xlab = "Time in Months",
        break.time.by = round(max(new_df$Survival, na.rm = TRUE)/5),
        ggtheme = theme_classic(),
        risk.table.y.text.col = TRUE,
        risk.table.height = 0.25,
        risk.table.y.text = FALSE,
        ncensor.plot = TRUE,
        ncensor.plot.height = 0.25,
        conf.int.style = "ribbon"
      )
    }
    p6$table <- p6$table + theme(legend.position = "none")
    p6$ncensor.plot <- p6$ncensor.plot + theme(legend.position = "none")
    
    ### beeswarm plot
    p7 <- ggplot(new_df, aes_string(x=colName, y="Survival")) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes_string(color=colName), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Survival (Months)") +
      theme(legend.position="top")
    
    # ### additional
    # g3 <- arrangeGrob(p6$plot, p6$table, p6$ncensor.plot,
    #                   layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 2), c(3, 3)),
    #                   top = paste0("Survival Differences between MSS AA vs CC"))
    # ggsave(file = paste0(outputDir, "/demo_MSS_AA_vs_CC_survival.png"), g3, width = 10, height = 8, dpi = 300)
    
    ### arrange the plots and print out
    g3 <- arrangeGrob(p6$plot, p6$table, p6$ncensor.plot, p7,
                      layout_matrix = rbind(c(1, 4), c(1, 4), c(1, 4), c(2, 4), c(3, 4)),
                      top = paste0("Survival Differences in Various ", colName))
    ggsave(file = paste0(outputDir, colName, "/demo_", colName, "_survival.png"), g3, width = 20, height = 12)
    
  }
  
  
  ### add more group info to the plot_df
  plot_df$MSI_H_Age_Distance <- NA
  plot_df$MSI_H_Age_Distance <- paste0(plot_df$MSI_AGE_Status, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSI_H_Age_Distance[!grepl("MSI-H", plot_df$MSI_H_Age_Distance)] <- NA
  plot_df$MSI_H_Age_Distance[grep("NA", plot_df$MSI_H_Age_Distance)] <- NA
  
  plot_df$MSS_Age_Distance <- NA
  plot_df$MSS_Age_Distance <- paste0(plot_df$MSI_AGE_Status, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSS_Age_Distance[!grepl("MSS", plot_df$MSS_Age_Distance)] <- NA
  plot_df$MSS_Age_Distance[grep("NA", plot_df$MSS_Age_Distance)] <- NA
  
  plot_df$MSI_H_Race1_Distance <- NA
  plot_df$MSI_H_Race1_Distance <- paste0(plot_df$MSI_RACE_Status1, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSI_H_Race1_Distance[!grepl("MSI-H", plot_df$MSI_H_Race1_Distance)] <- NA
  plot_df$MSI_H_Race1_Distance[grep("NA", plot_df$MSI_H_Race1_Distance)] <- NA
  
  plot_df$MSS_Race1_Distance <- NA
  plot_df$MSS_Race1_Distance <- paste0(plot_df$MSI_RACE_Status1, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSS_Race1_Distance[!grepl("MSS", plot_df$MSS_Race1_Distance)] <- NA
  plot_df$MSS_Race1_Distance[grep("NA", plot_df$MSS_Race1_Distance)] <- NA
  
  plot_df$MSI_H_Race2_Distance <- NA
  plot_df$MSI_H_Race2_Distance <- paste0(plot_df$MSI_RACE_Status2, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSI_H_Race2_Distance[!grepl("MSI-H", plot_df$MSI_H_Race2_Distance)] <- NA
  plot_df$MSI_H_Race2_Distance[grep("NA", plot_df$MSI_H_Race2_Distance)] <- NA
  
  plot_df$MSS_Race2_Distance <- NA
  plot_df$MSS_Race2_Distance <- paste0(plot_df$MSI_RACE_Status2, "_", plot_df$TUMOR_LOCATION)
  plot_df$MSS_Race2_Distance[!grepl("MSS", plot_df$MSS_Race2_Distance)] <- NA
  plot_df$MSS_Race2_Distance[grep("NA", plot_df$MSS_Race2_Distance)] <- NA
  
  ### run demographic analysis for each group
  demo_plots(plot_df, "MSI_AGE_Status")
  demo_plots(plot_df, "MSI_RACE_Status1")
  demo_plots(plot_df, "MSI_RACE_Status2")
  demo_plots(plot_df, "MSI_H_Age_Distance")
  demo_plots(plot_df, "MSS_Age_Distance")
  demo_plots(plot_df, "MSI_H_Race1_Distance")
  demo_plots(plot_df, "MSS_Race1_Distance")
  demo_plots(plot_df, "MSI_H_Race2_Distance")
  demo_plots(plot_df, "MSS_Race2_Distance")
  
}
