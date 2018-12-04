setClass("NanoStringRccSet",
         contains = "ExpressionSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(classVersion("ExpressionSet"),
                            NanoStringRccSet = "1.0.0"))))

setValidity("NanoStringRccSet",
function(object)
{
  msg <- NULL
  if (dim(object)[["Samples"]] > 0L) {
    # sampleNames
    if (!all(grepl("\\.rcc$", sampleNames(object), ignore.case = TRUE))) {
      msg <- c(msg, "'sampleNames' must all have an \".RCC\" file extension")
    }
    # phenoData
    phenoDataColNames <- c("SampleID", "Owner", "Comments", "Date", "SystemAPF")
    if (!all(phenoDataColNames %in% varLabels(phenoData(object)))) {
      msg <-
        c(msg,
          sprintf("'phenoData' must contain columns %s",
                  paste0("\"", phenoDataColNames, "\"", collapse = ", ")))
    }
    # protocolData
    protocolDataColNames <-
      c("LaneID", "FovCount", "FovCounted", "ScannerID", "StagePosition",
        "BindingDensity", "CartridgeID", "CartridgeBarcode", "FileVersion",
        "SoftwareVersion")
    if (!all(protocolDataColNames %in% varLabels(protocolData(object)))) {
      msg <-
        c(msg,
          sprintf("'protocolData' must contain columns %s",
                  paste0("\"", protocolDataColNames, "\"", collapse = ", ")))
    }
  }
  if (dim(object)[["Features"]] > 0L) {
    # featureData
    featureDataColNames <- c("CodeClass", "GeneName", "Accession")
    if (!all(featureDataColNames %in% varLabels(featureData(object)))) {
      msg <-
        c(msg,
          sprintf("'featureData' must contain columns %s",
                  paste0("\"", featureDataColNames, "\"", collapse = ", ")))
    }
  }
  if (sum(dim(object)) > 0L) {
    # annotation
    if (length(annotation(object)) != 1L || is.na(annotation(object)) ||
        !nzchar(annotation(object))) {
      msg <- c(msg, "'annotation' must contain the GeneRLF")
    }
  }
  if (prod(dim(object)) > 0L) {
    # assayData
    if (!is.integer(exprs(object)) || anyNA(exprs(object)) ||
        min(exprs(object)) < 0L) {
      msg <- c(msg, "'exprs' must be a non-negative integer matrix")
    }
  }
  if (is.null(msg)) TRUE else msg
})

setGeneric("NanoStringRccSet",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         ...)
  standardGeneric("NanoStringRccSet"),
signature = "assayData")


setMethod("NanoStringRccSet", "missing",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData, ...)
})

setMethod("NanoStringRccSet", "environment",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         ...)
{
  new("NanoStringRccSet",
      assayData = assayData,
      phenoData = phenoData,
      featureData = featureData,
      experimentData = experimentData,
      annotation = annotation,
      protocolData = protocolData,
      ...)
})

setMethod("NanoStringRccSet", "matrix",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData, ...)
})
