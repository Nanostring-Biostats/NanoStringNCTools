readNanoStringRccSet <-
function(rccFiles,
         rlfFile = NULL,
         phenoDataFile = NULL,
         phenoDataRccColName = "^RCC",
         phenoDataColPrefix = "")
{
  # Read data rccFiles
  data <- structure(lapply(rccFiles, readRccFile), names = basename(rccFiles))

  # Create assayData
  assay <-
    do.call(cbind,
            lapply(data, function(x)
              structure(x[["Code_Summary"]][["Count"]],
                        names = rownames(x[["Code_Summary"]]))))

  # Create phenoData
  if (is.null(phenoDataFile)) {
    pheno <- annotatedDataFrameFrom(assay, byrow = FALSE)
  } else {
    pheno <- read.csv(phenoDataFile, row.names = NULL, check.names = FALSE,
                      stringsAsFactors = FALSE)
    j <- grep(phenoDataRccColName, colnames(pheno), ignore.case = TRUE)
    if (length(j) == 0L)
      stop("Column `phenoDataRccColName` not found in `phenoDataFile`")
    else if (length(j) > 1L)
      stop("Multiple columns in `phenoDataFile` match `phenoDataRccColName`")
    missingPhenoCount <- sum(!(colnames(assay) %in% pheno[[j]]))
    rownames(pheno) <- pheno[[j]]
    pheno[[j]] <- NULL
    pheno <- pheno[colnames(assay), , drop = FALSE]
    if (missingPhenoCount != 0L) {
      rownames(pheno) <- colnames(assay)
      warning(sprintf("Column `phenoDataRccColName` in `phenoDataFile` is missing %d of %d Samples",
                      missingPhenoCount, ncol(assay)))
    }
    if (phenoDataColPrefix != "") {
      colnames(pheno) <- paste0(phenoDataColPrefix, colnames(pheno))
    }
    pheno <- AnnotatedDataFrame(pheno,
                                dimLabels = c("sampleNames", "sampleColumns"))
  }

  # Create featureData
  feature <- lapply(data, function(x) {
    x[["Code_Summary"]][, c("BarcodeClass", "GeneName", "Accession")]
  })
  stopifnot(all(sapply(feature, function(x) identical(feature[[1L]], x))))
  feature <- feature[[1L]]
  feature[["IsControl"]] <- NA
  feature[["IsControl"]][feature[["BarcodeClass"]] %in%
                           c("Negative", "Positive")] <- TRUE
  feature[["IsControl"]][feature[["BarcodeClass"]] %in%
                           c("Endogenous", "Housekeeping")] <- FALSE
  feature[["ControlConc"]] <- NA_real_
  feature[["ControlConc"]][feature[["BarcodeClass"]] == "Negative"] <- 0
  feature[["ControlConc"]][feature[["BarcodeClass"]] == "Positive"] <-
    as.numeric(sub("POS_\\w\\((.*)\\)", "\\1",
                   feature[["GeneName"]][feature[["BarcodeClass"]] ==
                                           "Positive"]))
  if (is.null(rlfFile)) {
    rlfHeader <- list()
  } else if (!is.null(rlfFile)) {
    rlfData <- readRlfFile(rlfFile)

    rlfHeader <- metadata(rlfData)
    rlfHeader[["RlfFileDate"]] <- as.character(rlfHeader[["RlfFileDate"]])

    rlfData <- as.data.frame(rlfData)
    rlfData <- rlfData[rlfData[["BarcodeClassActive"]] %in% c(2L, 3L), ,
                       drop = FALSE]
    rownames(rlfData) <-
      sprintf("%s_%s_%s", rlfData[["BarcodeClass"]], rlfData[["GeneName"]],
              rlfData[["Accession"]])

    if (!identical(sort(rownames(feature)), sort(rownames(rlfData))))
      stop("featureData mismatch between RLF and RCC files")

    rlfData <- rlfData[rownames(feature), , drop = FALSE]
    feature[["IsControl"]] <- rlfData[["BarcodeClassActive"]] == 3L

    for (j in c("BarcodeClass", "GeneName", "Accession", "BarcodeClassActive")) {
      rlfData[[j]] <- NULL
    }
    feature <- cbind(feature, rlfData)
  }
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

  # Create experimentData
  name <- sapply(data, function(x) x[["Sample_Attributes"]][["SampleOwner"]])
  name <- unique(na.omit(name))
  experiment <- MIAME(name = name, other = rlfHeader)

  # Create annotation
  annotation <- sapply(data, function(x) x[["Sample_Attributes"]][["GeneRLF"]])
  annotation <- unique(annotation)
  stopifnot(length(annotation) == 1L)

  # Create protocolData
  protocol <-
    do.call(rbind,
            lapply(seq_along(rccFiles), function(i) {
              x <- data[[i]][["Sample_Attributes"]]
              x <- x[, setdiff(names(x), "GeneRLF")]
              cbind(data[[i]][["Header"]], x, data[[i]][["Lane_Attributes"]])
            }))
  protocol <- AnnotatedDataFrame(protocol,
                                 .rccMetadata[["protocolData"]],
                                 dimLabels = c("sampleNames", "sampleColumns"))

  # Create NanoStringRccSet
  NanoStringRccSet(assayData = assay,
                   phenoData = pheno,
                   featureData = feature,
                   experimentData = experiment,
                   annotation = annotation,
                   protocolData = protocol,
                   check = FALSE)
}
