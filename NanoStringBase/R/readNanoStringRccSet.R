readNanoStringRccSet <-
function(rccFiles,
         rlfFile = NULL,
         phenoDataFile = NULL,
         phenoDataRccColName = "^RCC")
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
    if (!all(colnames(assay) %in% pheno[[j]]))
      stop("Column `phenoDataRccColName` in `phenoDataFile` does not match `rccFiles`")
    rownames(pheno) <- pheno[[j]]
    pheno[[j]] <- NULL
    pheno <- pheno[colnames(assay), , drop = FALSE]
    pheno <- AnnotatedDataFrame(pheno,
                                dimLabels = c("sampleNames", "sampleColumns"))
  }

  # Create featureData
  feature <- lapply(data, function(x) {
    x[["Code_Summary"]][, c("CodeClass", "GeneName", "Accession")]
  })
  stopifnot(all(sapply(feature, function(x) identical(feature[[1L]], x))))
  feature <- feature[[1L]]
  if (!is.null(rlfFile)) {
    rlfData <- as.data.frame(readRlfFile(rlfFile))
    rlfData <- rlfData[rlfData[["GeneName"]] %in% feature[["GeneName"]] &
                       rlfData[["Accession"]] %in% feature[["Accession"]], ,
                       drop = FALSE]
    rownames(rlfData) <-
      sprintf("%s_%s_%s", rlfData[["CodeClass"]], rlfData[["GeneName"]],
              rlfData[["Accession"]])
    for (j in c("CodeClass", "GeneName", "Accession")) {
      rlfData[[j]] <- NULL
    }
    rlfData <- rlfData[rownames(feature), , drop = FALSE]
    feature <- cbind(feature, rlfData)
  }
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

  # Create experimentData
  name <- sapply(data, function(x) x[["Sample_Attributes"]][["Owner"]])
  name <- unique(na.omit(name))
  experiment <- MIAME(name = name)

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
