###
#   File name : MCP_Counter_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Jun 10, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : RUN MCP-Counter to get abundance of tissue-infiltrating immune and
#               stromal cell populations
#
#   Instruction
#               1. Source("MCP_Counter_Analysis.R")
#               2. Run the function "mcp_counter" - specify the input files output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MCP_Counter_Analysis.R/MCP_Counter_Analysis.R")
#               > mcp_counter(rCntPath = "./data/rnaseq/raw_count_tcga_coad_read.rda",
#                             clinInfoPath="./data/coadread_tcga_clinical_data_updated2.txt",
#                             outputDir="./results/MCPcounter/")
###

mcp_counter <- function(rCntPath = "//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/raw_count_tcga_coad_read.rda",
                        clinInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/data/coadread_tcga_clinical_data_updated2.txt",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Parvathi_Myer/results/mcp_counter/") {
  
  ### load libraries
  options(java.parameters = "-Xmx8000m")
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
    library(devtools, quietly = TRUE)
  }
  if(!require(curl, quietly = TRUE)) {
    install.packages("curl")
    library(curl, quietly = TRUE)
  }
  if(!require(MCPcounter, quietly = TRUE)) {
    install_github("ebecht/MCPcounter",ref="master", subdir="Source")
    library(MCPcounter, quietly = TRUE)
  }
  if(!require(gplots)) {
    install.packages("gplots")
    library(gplots)
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
  if(!require(ComplexHeatmap, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap", version = "3.8")
    require(ComplexHeatmap, quietly = TRUE)
  }
  
  ### load the data
  load(rCntPath)
  clinInfo <- read.table(clinInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinInfo) <- clinInfo$`Sample ID`
  
  ### rownames of clinInfo == colnames of norm_gep
  clinInfo <- clinInfo[colnames(raw_count_tcga_coad_read),]
  
  ### remove POLE-mutated samples
  raw_count_tcga_coad_read <- raw_count_tcga_coad_read[,-which(clinInfo$POLE_MUTANT == TRUE)]
  clinInfo <- clinInfo[-which(clinInfo$POLE_MUTANT == TRUE),]
  
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
  
  ### msi status Indeterminate -> NA
  clinInfo$NEW_MSI[which(clinInfo$NEW_MSI == "Indeterminate")] <- NA
  
  ### add age categorical info (Old/Young)
  clinInfo$Age <- NA
  clinInfo$Age[which(clinInfo$`Diagnosis Age` >= 50)] <- "Old"
  clinInfo$Age[which(clinInfo$`Diagnosis Age` < 50)] <- "Young"
  
  ### add race categorical info
  clinInfo$Self_Reported_Race <- NA
  clinInfo$Self_Reported_Race[which(clinInfo$`Race Category` == "BLACK OR AFRICAN AMERICAN")] <- "AA"
  clinInfo$Self_Reported_Race[which(clinInfo$`Race Category` == "WHITE")] <- "CC"
  clinInfo$Predicted_Race <- NA
  clinInfo$Predicted_Race[which(clinInfo$Prediction_Filtered == "African")] <- "AA"
  clinInfo$Predicted_Race[which(clinInfo$Prediction_Filtered == "Caucasian")] <- "CC"
  
  ### add more group info (location-related) to the clinInfo
  clinInfo$MSI_H_Age_Distance <- NA
  clinInfo$MSI_H_Age_Distance <- paste0(clinInfo$MSI_AGE_Status, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSI_H_Age_Distance[!grepl("MSI-H", clinInfo$MSI_H_Age_Distance)] <- NA
  clinInfo$MSI_H_Age_Distance[grep("NA", clinInfo$MSI_H_Age_Distance)] <- NA
  
  clinInfo$MSS_Age_Distance <- NA
  clinInfo$MSS_Age_Distance <- paste0(clinInfo$MSI_AGE_Status, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSS_Age_Distance[!grepl("MSS", clinInfo$MSS_Age_Distance)] <- NA
  clinInfo$MSS_Age_Distance[grep("NA", clinInfo$MSS_Age_Distance)] <- NA
  
  clinInfo$MSI_H_Race1_Distance <- NA
  clinInfo$MSI_H_Race1_Distance <- paste0(clinInfo$MSI_RACE_Status1, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSI_H_Race1_Distance[!grepl("MSI-H", clinInfo$MSI_H_Race1_Distance)] <- NA
  clinInfo$MSI_H_Race1_Distance[grep("NA", clinInfo$MSI_H_Race1_Distance)] <- NA
  
  clinInfo$MSS_Race1_Distance <- NA
  clinInfo$MSS_Race1_Distance <- paste0(clinInfo$MSI_RACE_Status1, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSS_Race1_Distance[!grepl("MSS", clinInfo$MSS_Race1_Distance)] <- NA
  clinInfo$MSS_Race1_Distance[grep("NA", clinInfo$MSS_Race1_Distance)] <- NA
  
  clinInfo$MSI_H_Race2_Distance <- NA
  clinInfo$MSI_H_Race2_Distance <- paste0(clinInfo$MSI_RACE_Status2, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSI_H_Race2_Distance[!grepl("MSI-H", clinInfo$MSI_H_Race2_Distance)] <- NA
  clinInfo$MSI_H_Race2_Distance[grep("NA", clinInfo$MSI_H_Race2_Distance)] <- NA
  
  clinInfo$MSS_Race2_Distance <- NA
  clinInfo$MSS_Race2_Distance <- paste0(clinInfo$MSI_RACE_Status2, "_", clinInfo$TUMOR_LOCATION)
  clinInfo$MSS_Race2_Distance[!grepl("MSS", clinInfo$MSS_Race2_Distance)] <- NA
  clinInfo$MSS_Race2_Distance[grep("NA", clinInfo$MSS_Race2_Distance)] <- NA
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
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
  
  ### normalize the raw counts with VST
  norm_gep <- normalizeRNASEQwithVST(raw_count_tcga_coad_read, filtering = FALSE)
  
  ### run MCP-Counter
  mcp_result <- MCPcounter.estimate(expression = norm_gep, featuresType = "HUGO_symbols")
  
  ### write out the result
  write.xlsx2(data.frame(Cell_Type=rownames(mcp_result), mcp_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "TCGA_COAD_READ_MCP-Counter_Result.xlsx"),
              sheetName = "MCP-Counter_Result", row.names = FALSE)
  
  ###
  #   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
  ###
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  
  ### set colside colors
  clinInfo[is.na(clinInfo)] <- "NA"
  uniqueV <- unique(c(clinInfo$NEW_MSI, clinInfo$Age, clinInfo$Self_Reported_Race, clinInfo$Predicted_Race, clinInfo$TUMOR_LOCATION))
  colors <- c("black", "white", "gray", "blue", "red", "purple", "green", "darkcyan", "orange")
  names(colors) <- uniqueV
  
  ### heatmap
  fName <- "Heatmap_MCP-Counter_Result"
  png(paste0(outputDir, fName, ".png"), width = 2000, height = 1400, res = 130)
  par(oma=c(0,0,0,8))
  heatmap.3(mcp_result, main = fName,
            distfun = function(x) dist(x, method = "euclidean"),
            hclustfun = function(x) hclust(x, method = "ward.D"),
            xlab = "", ylab = "", col=greenred(100),
            scale="row", key=T, keysize=0.8, dendrogram = 'none', trace = 'none',
            labRow = rownames(mcp_result), labCol = "",
            Rowv = TRUE, Colv = TRUE,
            ColSideColors = cbind(colors[clinInfo$NEW_MSI],
                                  colors[clinInfo$Age],
                                  colors[clinInfo$Self_Reported_Race],
                                  colors[clinInfo$Predicted_Race],
                                  colors[clinInfo$TUMOR_LOCATION]),
            ColSideColorsSize = 3,
            cexRow = 1.3, cexCol = 1, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Sample Annotation", legend = names(colors), fill = colors, cex = 0.7, box.lty = 1)
  legend("topright", inset = -0.02, xpd = TRUE, title = "The 5 Horizontal Bars", legend = c("Tumor Location", "Predicted Race", "Self-reported Race", "Age", "MSI Status"), fill = "white", cex = 0.7, box.lty = 1)
  dev.off()
  
  ### beeswarm plot
  mcp_beeswarm <- function(group) {
    ### make an empty plot list
    p <- vector("list", length = nrow(mcp_result))
    
    ### iteratively draw a plot for every cell types
    for(i in 1:length(p)) {
      ### create a data frame for the plot
      plot_df <- data.frame(MCP_Counter_Score=mcp_result[i,], Sample_Group=clinInfo[,group],
                            stringsAsFactors = FALSE, check.names = FALSE)
      
      ### remove NA rows in the group
      plot_df <- plot_df[which(plot_df$Sample_Group != "NA"),]
      
      ### additional
      # plot_df <- plot_df[union(which(plot_df$Sample_Group == "MSS_AA"),
      #                          which(plot_df$Sample_Group == "MSS_CC")),]
      
      ### draw a plot with the data frame
      p[[i]] <- ggplot(plot_df, aes(x=Sample_Group, y=MCP_Counter_Score)) +
        ggtitle(rownames(mcp_result)[i]) +
        theme_classic(base_size = 10) +
        geom_boxplot() +
        geom_beeswarm(aes(color=Sample_Group), na.rm = TRUE) +
        stat_compare_means() +
        labs(x = "", y = "MCP-Counter Score") +
        theme(legend.position="none", plot.title=element_text(hjust = 0.5))
    }
    
    ### arrange the plots and print out
    fName <- paste0("Beeswarm_Plot_MCP-Counter_Result_", group)
    g <- arrangeGrob(grobs = p,
                     layout_matrix = rbind(c(1, 2, 3, 4, 5), c(6, 7, 8, 9, 10)),
                     top = fName)
    ggsave(file = paste0(outputDir, fName, ".png"), g, width = 22, height = 12)
  }
  
  ### run the mcp_beeswarm() for each group
  mcp_beeswarm("MSI_AGE_Status")
  mcp_beeswarm("MSI_RACE_Status1")
  mcp_beeswarm("MSI_RACE_Status2")
  mcp_beeswarm("MSI_H_Age_Distance")
  mcp_beeswarm("MSS_Age_Distance")
  mcp_beeswarm("MSI_H_Race1_Distance")
  mcp_beeswarm("MSS_Race1_Distance")
  mcp_beeswarm("MSI_H_Race2_Distance")
  mcp_beeswarm("MSS_Race2_Distance")
  mcp_beeswarm("Age")
  mcp_beeswarm("Self_Reported_Race")
  mcp_beeswarm("Predicted_Race")
  
}
