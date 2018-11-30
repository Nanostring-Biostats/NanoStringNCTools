readNanoStringRccSet <-
function(rccFiles, rlfFile = NULL)
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
  pheno <-
    do.call(rbind,
            lapply(seq_along(rccFiles), function(i) {
              x <- data[[i]][["Sample_Attributes"]]
              x[, setdiff(names(x), "GeneRLF")]
            }))
  pheno <- AnnotatedDataFrame(pheno,
                              dimLabels = c("sampleNames", "sampleColumns"))

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
              cbind(data[[i]][["Lane_Attributes"]], data[[i]][["Header"]])
            }))
  protocol <- AnnotatedDataFrame(protocol,
                                 dimLabels = c("sampleNames", "sampleColumns"))

  # Create NanoStringRccSet
  NanoStringRccSet(assayData = assay,
                   phenoData = pheno,
                   featureData = feature,
                   experimentData = experiment,
                   annotation = annotation,
                   protocolData = protocol)
}
