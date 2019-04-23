###
#   File name : ScriptForMSIsensor.R
#   Author    : Hyunjin Kim
#   Date      : Apr 22, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : This will generates a shell script that downloads appropriate Bam files
#               and runs MSIsensor to compute MSIsensor score. But the point is,
#               while downloading Bam files, the MSIsensor should be run simultaneously
#               on the already downloaded data.
#
#   Instruction
#               1. Source("ScriptForMSIsensor.R")
#               2. Run the function "createScript" - specify the file path of the TCGA WXS file info,
#                  other necessary GDC paths, output directory for the MSIsensor, and script output path
#                  The first one will be used on the current Windows OS, and the others will be
#                  in the script that will be run on linux.
#               3. The results will be generated in the script output path
#
#   Example
#               > source("The_directory_of_ScriptForMSIsensor.R/ScriptForMSIsensor.R")
#               > createScript(fileInfoPath = "C:/Research/CUMC/CRC/data/MSIsensor/WXS_matched_tumor_normal_sample_info.txt",
#                              dos2unixPath="C:/cygwin64/bin/dos2unix.exe",
#                              scriptPath="C:/Research/CUMC/CRC/data/MSIsensor/MSIsensor_script.sh",
#                              gdcToolPath="/mnt/c/Research/gdc-client_v1.4.0/gdc-client_v1.4.0_Ubuntu_x64",
#                              gdcTokenPath="/mnt/c/Research/CUMC/CRC/data/gdc-user-token.2019-04-22T15_26_09.730Z.txt",
#                              msListPath="/mnt/c/Research/CUMC/CRC/data/msi_scan_list/microsatellites_hg38.list",
#                              outputPath="/mnt/c/Research/CUMC/CRC/results/MSIsensor/")
###

createScript <- function(fileInfoPath = "C:/Research/CUMC/CRC/data/MSIsensor/WXS_matched_tumor_normal_sample_info.txt",
                         dos2unixPath="C:/cygwin64/bin/dos2unix.exe",
                         scriptPath="C:/Research/CUMC/CRC/data/MSIsensor/MSIsensor_script.sh",
                         gdcToolPath="/mnt/c/Research/gdc-client_v1.4.0/gdc-client_v1.4.0_Ubuntu_x64",
                         gdcTokenPath="/mnt/c/Research/CUMC/CRC/data/gdc-user-token.2019-04-22T15_26_09.730Z.txt",
                         msListPath="/mnt/c/Research/CUMC/CRC/data/msi_scan_list/microsatellites_hg38.list",
                         outputPath="/mnt/c/Research/CUMC/CRC/results/MSIsensor/") {
  
  ### load the file info
  fileInfo <- read.table(file = fileInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### script's first line
  script <- "#!/bin/bash\n\n"
  
  ### set environmental parameters
  script <- paste0(script, "GDCTOOL=", gdcToolPath, "\n")
  script <- paste0(script, "GDCTOKEN=", gdcTokenPath, "\n")
  script <- paste0(script, "WD=", outputPath, "\n")
  script <- paste0(script, "\n")
  
  ### for all the WXS Bam files, download and run MSIsensor
  for(i in 1:nrow(fileInfo)) {
    ### some patients have more than one Bam files
    ### just use the first one - in normal, they are always samples from blood
    normalFileName <- strsplit(fileInfo[i,"Normal_File_Name"], split = ",", fixed = TRUE)[[1]][1]
    normalFileID <- strsplit(fileInfo[i,"Normal_File_ID"], split = ",", fixed = TRUE)[[1]][1]
    tumorFileName <- strsplit(fileInfo[i,"Tumor_File_Name"], split = ",", fixed = TRUE)[[1]][1]
    tumorFileID <- strsplit(fileInfo[i,"Tumor_File_ID"], split = ",", fixed = TRUE)[[1]][1]
    
    ### add a comment to each loop
    script <- paste0(script, "### ", fileInfo[i,"Case_ID"], "\n")
    
    ### download using GDC transfer tool
    script <- paste0(script, "${GDCTOOL} download ", normalFileID, " ",
                     tumorFileID, " -t ${GDCTOKEN} -d ${WD}\n")
    
    ### move the Bam files out from the sub-directory
    ### and for *.bai files, rename them as *.bam.bai
    ### MSIsensor reconizes index files as *.bam.bai
    fileName <- substr(normalFileName, 1, nchar(normalFileName)-4)
    script <- paste0(script, "cp ${WD}", normalFileID, "/",
                     normalFileName, " ${WD}",
                     normalFileName, "\n")
    script <- paste0(script, "cp ${WD}", normalFileID, "/",
                     fileName, ".bai ${WD}", fileName, ".bam.bai\n")
    script <- paste0(script, "rm -rf ${WD}", normalFileID, "\n\n")
    fileName <- substr(tumorFileName, 1, nchar(tumorFileName)-4)
    script <- paste0(script, "cp ${WD}", tumorFileID, "/",
                     tumorFileName, " ${WD}",
                     tumorFileName, "\n")
    script <- paste0(script, "cp ${WD}", tumorFileID, "/",
                     fileName, ".bai ${WD}", fileName, ".bam.bai\n")
    script <- paste0(script, "rm -rf ${WD}", tumorFileID, "\n\n")
    
    ### run MSIsensor
    script <- paste0(script, "msisensor msi -d ", msListPath, " -n ${WD}",
                     normalFileName, " - t ${WD}",
                     tumorFileName, " -o ${WD}", fileInfo[i,"Case_ID"], "\n\n")
    
    ### delete the used Bam files
    script <- paste0(script, "rm ${WD}", normalFileName, "\n")
    script <- paste0(script, "rm ${WD}", tumorFileName, "\n")
    script <- paste0(script, "\n")
    
    ### print progress
    script <- paste0(script, "echo \"MSIsensor run on ", fileInfo[i,"Case_ID"], " is done - ",
                     i, "/", nrow(fileInfo), "\"\n\n\n")
  }
  
  ### save the script
  write.table(script, file=scriptPath, sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### dos2unix to the result
  system(paste(dos2unixPath, scriptPath))
  
}
