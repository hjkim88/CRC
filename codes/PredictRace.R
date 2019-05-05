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
  
  
  
  
  
  
  
}




