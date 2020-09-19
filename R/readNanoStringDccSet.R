readNanoStringDccSet <-
function(dccFiles,
         pkcFile = NULL,
         phenoDataFile = NULL,
         phenoDataSheet = NULL,
         phenoDataDccColName = "Sample_ID",
         phenoDataColPrefix = "")
{
  # Read data rccFiles
  data <- structure(lapply(dccFiles, readDccFile), names = basename(dccFiles))

  # Create assayData
  assay <- lapply(data, function(x)
    structure(x[["Code_Summary"]][["Count"]],
              names = rownames(x[["Code_Summary"]])))

  # Create phenoData
  if (is.null(phenoDataFile)) {
    pheno <- annotatedDataFrameFrom(assay, byrow = FALSE)
  } else {
    pheno <- read_xlsx(phenoDataFile, col_names = TRUE, sheet = phenoDataSheet)
    pheno <- data.frame(pheno, stringsAsFactors = FALSE)
    j <- grep(phenoDataDccColName, colnames(pheno), ignore.case = TRUE)
    if (length(j) == 0L){
      stop("Column `phenoDataDccColName` not found in `phenoDataFile`")
    } else if (length(j) > 1L){
      stop("Multiple columns in `phenoDataFile` match `phenoDataDccColName`")
    }
    missingPhenoCount <- sum(!(colnames(assay) %in% pheno[[j]]))
    pheno[[j]] <- paste0(pheno[[j]], ".dcc")
    rownames(pheno) <- pheno[[j]]
    pheno[[j]] <- NULL
    pheno <- pheno[names(assay), , drop = FALSE]
    if (missingPhenoCount != 0L) {
      rownames(pheno) <- colnames(assay)
      warning(sprintf("Column `phenoDataDccColName` in `phenoDataFile` is missing %d of %d Samples",
                      missingPhenoCount, ncol(assay)))
    }
    if (phenoDataColPrefix != "") {
      colnames(pheno) <- paste0(phenoDataColPrefix, colnames(pheno))
    }
    pheno <- AnnotatedDataFrame(pheno,
                                dimLabels = c("sampleNames", "sampleColumns"))
  }

  #stopifnot(all(sapply(feature, function(x) identical(feature[[1L]], x))))
  if (is.null(pkcFile)) {
    pkcHeader <- list()
  } else if (!is.null(pkcFile)) {
    pkcData <- readPKCFile(pkcFile)

    pkcHeader <- metadata(pkcData)
    pkcHeader[["PKCFileDate"]] <- as.character(pkcHeader[["PKCFileDate"]])

    pkcData <- as.data.frame(pkcData)
    rownames(pkcData) <- pkcData[["RTS_ID"]]

  }
  
  for (index in seq_len(length(data))) {
    countMat <- data[[index]]$Code_Summary
    countMat <- countMat[which(rownames(countMat) %in% rownames(pkcData)), , drop = FALSE]
    countMat$Gene <- pkcData$Gene[match(countMat$RNAID, 
                                        pkcData$RTS_ID)]
    
    countMat$Pool <- pkcData$Module[match(countMat$RNAID, 
                                          pkcData$RTS_ID)]
    data[[index]]$Code_Summary <- countMat
  }
  
  probe_assay <- lapply(1:length(data), function(x)
    data.frame(data[[x]][["Code_Summary"]],
               Sample_ID = names(data)[x]))
  probe_assay <- do.call(rbind, probe_assay)
  
  gene_assay <- dcast(probe_assay, Gene + Pool ~ Sample_ID, 
                      value.var = 'Count', fun.aggregate = ngeoMean, fill = 1)
  
  if ( length(unique(gene_assay$Gene)) != nrow(gene_assay) ) {
    stop('Some genes are listed in multiple pools.')
  }
  
  rownames(gene_assay) <- gene_assay[, "Gene"]
  assay <- as.matrix(gene_assay[, -c(1:2)])
  
  # Create featureData
  feature <- gene_assay[, "Gene", drop = FALSE]
  rownames(feature) <- feature[["Gene"]]
  
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

  # Create experimentData
  construct <- unique(na.omit(pheno@data[["construct"]]))
  instrument_type <- unique(na.omit(pheno@data[["instrument_type"]]))
  read_pattern <- unique(na.omit(pheno@data[["read_pattern"]]))
  panel <- unique(na.omit(pheno@data[["panel"]]))
  pooling <- unique(na.omit(pheno@data[["pooling"]]))
  
  experiment <- MIAME(name = "", 
                      other = list(construct = construct, 
                                   instrument_type = instrument_type,
                                   read_pattern = read_pattern,
                                   panel = panel,
                                   pooling = pooling))

  # Create annotation
  annotation <- sort(sapply(strsplit(pkcFile, "/"), function(x) x[length(x)]))
  if( !identical(annotation, paste0(sort(unique(probe_assay[['Pool']])), ".pkc")) ) {
    stop("Name mismatch between pool and PKC files")
  }

  # Create protocolData
  protocol <-
    do.call(rbind,
            lapply(seq_along(dccFiles), function(i) {
              cbind(data[[i]][["Header"]], data[[i]][["Scan_Attributes"]])
            }))
  
  protocol <- data.frame(protocol, 
                         pheno@data[, which(colnames(pheno@data) %in% 
                            c("dsp_scan", "dsp_collection_plate", "dsp_collection_well", 
                              "pcr_primer_plate", "pcr_primer_well"))])
  
  annot_labelDescription <-  data.frame(labelDescription =
                                           c(NA_character_,
                                             NA_character_,
                                             NA_character_,
                                             NA_character_,
                                             NA_character_),
                                         row.names =
                                           c("dsp_scan", "dsp_collection_plate", "dsp_collection_well", 
                                             "pcr_primer_plate", "pcr_primer_well"),
                                         stringsAsFactors = FALSE)
  
  protocol <- AnnotatedDataFrame(protocol,
                                 rbind(.dccMetadata[["protocolData"]], 
                                       annot_labelDescription),
                                 dimLabels = c("sampleNames", "sampleColumns"))

  # Create NanoStringDccSet
  NanoStringDccSet(assayData = assay,
                   phenoData = pheno,
                   featureData = feature,
                   experimentData = experiment,
                   annotation = annotation,
                   protocolData = protocol,
                   check = FALSE)
}
