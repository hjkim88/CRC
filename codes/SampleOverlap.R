###
#   File name : SampleOverlap.R
#   Author    : Hyunjin Kim
#   Date      : Apr 8, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are two CRC data in cBioPortal.
#               Determine if they are actually different patients or an overlap of the same patients
#               Draw Venn diagrams of how they are overlapped 
#
#   Instruction
#               1. Source("SampleOverlap.R")
#               2. Run the function "checkOverlap" - specify two clinical info files and output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_SampleOverlap.R/SampleOverlap.R")
#               > checkOverlap(clinInfoPath1 = "./data/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv",
#                              clinInfoPath2 = "./data/coadread_tcga_clinical_data.tsv",
#                              outputDir = "./results/")
###

checkOverlap <- function(clinInfoPath1 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv",
                         clinInfoPath2 = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data.tsv",
                         outputDir = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/") {
  
  ### load required library
  if(!require(VennDiagram, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("VennDiagram")
    require(VennDiagram, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  ### load the clinical info
  clinInfo1 <- read.table(file = clinInfoPath1, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
  clinInfo2 <- read.table(file = clinInfoPath2, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the venn diagram
  v <- venn.diagram(list(clinInfo1$`Sample ID`, clinInfo2$`Sample ID`),
                    category.names = c("PanCancer Atlas", "Provisional"),
                    cat.cex = 1.0, cex = 1.5,
                    filename = NULL)
  png(paste0(outputDir, "sample_overlap.png"),
      width = 1200, height = 1000, res = 180)
  grid.arrange(gTree(children=v),
               top=paste0("Sample Overlap of between PanCancer_Atlas & Provisional"),
               bottom="")
  dev.off()
  
}


