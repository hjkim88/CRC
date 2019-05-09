###
#   File name : MakeGenoRDA.R
#   Author    : Hyunjin Kim
#   Date      : Apr 30, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : The genotype data of TCGA COAD/READ are downloaded.
#               The reference genotype data are also downloaded from 1000 Genomes Project.
#               Combine the genotype data and make a RDA file.
#               But when combining the data, only keep the common SNP between them.
#
#   Instruction
#               1. Source("MakeGenoRDA.R")
#               2. Run the function "make_genotype_rda" - specify genotype file path, reference list, and output result path
#               3. The result will be generated in the output result path
#
#   Example
#               > source("The_directory_of_MakeGenoRDA.R/MakeGenoRDA.R")
#               > make_genotype_rda(genoPath="./data/genotyping/TCGA/",
#                                   metaDataPath="./data/genotyping/metadata.cart.2019-04-26.json",
#                                   refSampleDataPath="./data/genotyping/20130606_g1k.ped",
#                                   genoSNPListPath="./data/genotyping/TCGA.variants.txt",
#                                   referenceSNPListPath="./data/genotyping/Merge3.variants.txt",
#                                   referenceGenoPath="./data/genotyping/Merge3.012.txt",
#                                   outputResultPath="./data/genotyping/1000g_tcga_combined_genotype.RDA")
###

### THE FILE SIZE IS TOO LARGE, THEREFORE, I DID NOT UPLOAD THE DATA ON THE SERVER,
### AND JUST RAN IT ON MY LOCAL COMPUTER
make_genotype_rda <- function(genoPath="./data/genotyping/TCGA/",
                              metaDataPath="./data/genotyping/metadata.cart.2019-04-26.json",
                              refSampleDataPath="./data/genotyping/20130606_g1k.ped",
                              genoSNPListPath="./data/genotyping/TCGA.variants.txt",
                              referenceSNPListPath="./data/genotyping/Merge.variants.txt",
                              referenceGenoPath="./data/genotyping/Merge.012.txt",
                              outputResultPath="./data/genotyping/1000g_tcga_combined_genotype.RDA") {
  
  ### load library
  if(!require(rjson, quietly = TRUE)) {
    install.packages("rjson")
    require(rjson, quietly = TRUE)
  }
  if(!require(data.table, quietly = TRUE)) {
    install.packages("data.table")
    require(data.table, quietly = TRUE)
  }
  if(!require(pd.genomewidesnp.6, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("pd.genomewidesnp.6", version = "3.8")
    require(pd.genomewidesnp.6, quietly = TRUE)
  }
  
  ### load the TCGA genotype data SNP list (CHECK ALL THE GENO FILES HAVE THE SAME ROWS AS THIS ONE)
  genoList <- read.table(file = genoSNPListPath, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE, check.names = FALSE)[,1]
  
  ### load the reference SNP list
  refList <- read.table(file = referenceSNPListPath, sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE, check.names = FALSE)[,1]
  refList2 <- sapply(refList, function(x) paste0(strsplit(x, ":", fixed = TRUE)[[1]][1:2], collapse = ":"), USE.NAMES = FALSE)
  
  ### get the mapping info between Affymetrix SNP ID and CHR:POS
  con <- db(pd.genomewidesnp.6)
  db <- dbGetQuery(con, "SELECT man_fsetid, dbsnp_rs_id, chrom, physical_pos, strand, allele_a, allele_b from featureSet")
  db <- db[which(!is.na(db$dbsnp_rs_id)),]
  rownames(db) <- db[,"man_fsetid"]
  
  ### get CHR:POS version of the TCGA genotype data SNP list
  genoList_chr_pos <- paste0(db[genoList,"chrom"], ":", db[genoList,"physical_pos"])
  
  ### SNP list - TCGA data vs Reference - keep the commons only
  geno_retain_idx <- which(genoList_chr_pos %in% refList2)
  ref_retain_idx <- which(refList2 %in% genoList_chr_pos)
  
  ### remove duplicates - here, I got rid of all the duplicates including the original ones
  ### because I have no info which one I should retain if there are duplicates
  dups <- which(duplicated(genoList_chr_pos[geno_retain_idx]))
  if(length(dups) > 0) {
    rmv_obj <- enoList_chr_pos[geno_retain_idx][dups]
    geno_retain_idx <- geno_retain_idx[-which(genoList_chr_pos[geno_retain_idx] %in% rmv_obj)]
    ref_retain_idx <- ref_retain_idx[-which(refList2[ref_retain_idx] %in% rmv_obj)]
  }
  dups <- which(duplicated(refList2[ref_retain_idx]))
  if(length(dups) > 0) {
    rmv_obj <- refList2[ref_retain_idx][dups]
    ref_retain_idx <- ref_retain_idx[-which(refList2[ref_retain_idx] %in% rmv_obj)]
    geno_retain_idx <- geno_retain_idx[-which(genoList_chr_pos[geno_retain_idx] %in% rmv_obj)]
  }
  
  ### order the idx vectors based on the SNP list
  geno_retain_idx <- geno_retain_idx[order(genoList_chr_pos[geno_retain_idx])]
  ref_retain_idx <- ref_retain_idx[order(refList2[ref_retain_idx])]
  
  ### check if the two list are the same
  if(identical(genoList_chr_pos[geno_retain_idx], refList2[ref_retain_idx])) {
    writeLines("The two SNP list are the same. Now you can proceed.")
  } else {
    writeLines("The two SNP list are NOT the same. There is something wrong with them. Do not proceed further.")
  }
  
  ### get the sample info of the genotype data files from the TCGA metadata
  temp <- fromJSON(file = metaDataPath)
  tcga_meta <- data.frame(matrix(NA, length(temp), 3))
  colnames(tcga_meta) <- c("Sample_ID", "Gender", "Ethnicity")
  rownames(tcga_meta) <- sapply(temp, function(x) {
    if(is.null(x$file_name)) {
      return(NA)
    } else {
      return(x$file_name)
    }
  })
  tcga_meta[,"Sample_ID"] <- sapply(temp, function(x) {
    if(is.null(x$associated_entities[[1]]$entity_submitter_id)) {
      return(NA)
    } else {
      return(x$associated_entities[[1]]$entity_submitter_id)
    }
  })
  tcga_meta[,"Gender"] <- sapply(temp, function(x) {
    if(is.null(x$cases[[1]]$demographic$gender)) {
      return(NA)
    } else {
      return(x$cases[[1]]$demographic$gender)
    }
  })
  tcga_meta[,"Ethnicity"] <- sapply(temp, function(x) {
    if(is.null(x$cases[[1]]$demographic$race)) {
      return(NA)
    } else {
      return(toupper(x$cases[[1]]$demographic$race))
    }
  })
  
  ### get the sample info of the reference data
  ref_meta <- read.table(file = refSampleDataPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(ref_meta) <- ref_meta[,"Individual ID"]
  
  ### clean up some used variables that will not be used anymore
  rm(db)
  rm(genoList)
  rm(genoList_chr_pos)
  rm(refList)
  rm(refList2)
  rm(con)
  rm(rmv_obj)
  rm(temp)
  gc()
  
  ### load the 1000 Genomes Project reference genotype data
  refGeno <- fread(input = referenceGenoPath, sep = " ", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  refGeno <- transpose(refGeno)
  setDF(refGeno)
  
  ### get samples IDs
  ref_sample_ids <- as.character(refGeno[2,2:ncol(refGeno)])
  
  ### only keep the ones in the TCGA samples
  ### first 6 rows are sample information
  refGeno <- refGeno[ref_retain_idx+6,]
  combined_geno <- refGeno
  rm(refGeno)
  gc()
  
  ### make row names and column names
  rownames(combined_geno) <- sapply(combined_geno[,1], function(x) paste0(strsplit(x, ":", fixed = TRUE)[[1]][1:2], collapse = ":"), USE.NAMES = FALSE)
  combined_geno <- combined_geno[,-1]
  colnames(combined_geno) <- ref_sample_ids
  
  ### get genotype file list
  f <- list.files(genoPath)
  fNames <- sapply(f, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1], USE.NAMES = FALSE)
  
  ### load the genotype data
  ### we already know that all the files have the same rows
  for(i in 1:length(f)) {
    temp <- read.table(file = paste0(genoPath, f[i]), header = TRUE, sep = "\t",
                       row.names = 1, skip = 1, stringsAsFactors = FALSE, check.names = FALSE)[,1][geno_retain_idx]
    combined_geno <- data.frame(combined_geno, temp, stringsAsFactors = FALSE, check.names = FALSE)
    writeLines(paste(i, "/", length(f)))
    gc()
  }
  colnames(combined_geno)[(length(ref_sample_ids)+1):ncol(combined_geno)] <- tcga_meta[f,"Sample_ID"]
  temp <- rownames(combined_geno)
  combined_geno <- data.frame(lapply(combined_geno, as.numeric), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(combined_geno) <- temp
  
  ### create sample info for the combined data
  sample_info <- data.frame(matrix("", ncol(combined_geno), 5), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sample_info) <- colnames(combined_geno)
  colnames(sample_info) <- c("Source", "File_Name", "Sample_ID", "Gender", "Ethnicity")
  
  ### fill the sample info with tcga_meta and ref_meta
  sample_info[1:length(ref_sample_ids),"Source"] <- "1000 Genomes Project"
  sample_info[(length(ref_sample_ids)+1):nrow(sample_info),"Source"] <- "TCGA"
  sample_info[1:length(ref_sample_ids),"File_Name"] <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  sample_info[(length(ref_sample_ids)+1):nrow(sample_info),"File_Name"] <- rownames(tcga_meta[f,])
  sample_info[,"Sample_ID"] <- colnames(combined_geno)
  sample_info[1:length(ref_sample_ids),"Gender"] <- sapply(ref_meta[ref_sample_ids,"Gender"], function(x) {
    if(x == 1) {
      return("male")
    } else if(x == 2) {
      return("female")
    } else {
      return(NA)
    }
  })
  sample_info[(length(ref_sample_ids)+1):nrow(sample_info),"Gender"] <- tcga_meta[f,"Gender"]
  sample_info[1:length(ref_sample_ids),"Ethnicity"] <- ref_meta[ref_sample_ids,"Population"]
  sample_info[(length(ref_sample_ids)+1):nrow(sample_info),"Ethnicity"] <- tcga_meta[f,"Ethnicity"]
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"combined_geno\" is a data frame of the 1000 Genomes Project - TCGA combined")
    writeLines("genotype data. There are 36544 rows (SNPs) and 3814 columns (samples).")
    writeLines("Among the 3814 samples, the first 2504 samples are from 1000 Genome Project and")
    writeLines("can be used as reference while the later 1310 samples are from TCGA - the test samples.")
    writeLines("0: homozygous reference, 1: heterozygous, 2: homozygous alternate")
    writeLines("The \"sample_info\" is a data frame of sample information of the patients in \"combined_geno\".")
    writeLines("The rownames(sample_info) is same as colnames(combined_geno).")
    writeLines("It describes where the info come from, file names, gender, and ethnicity of the samples.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the distance matrices in RDA file
  save(list = c("combined_geno", "sample_info", "README"), file = outputResultPath)
  
}
