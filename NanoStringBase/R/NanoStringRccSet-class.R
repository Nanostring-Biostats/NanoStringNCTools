setClass("NanoStringRccSet",
         contains = "ExpressionSet",
         slots = c(signatureWeights = "NumericList"),
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(classVersion("ExpressionSet"),
                            NanoStringRccSet = "1.0.0")),
           signatureWeights = NumericList()))

setMethod("initialize", "NanoStringRccSet",
function(.Object, signatureWeights = NumericList(), ...)
{
  callNextMethod(.Object,
                 signatureWeights = as(signatureWeights, "NumericList"),
                 ...)
})

setValidity2("NanoStringRccSet",
function(object)
{
  msg <- NULL
  if (dim(object)[["Samples"]] > 0L) {
    # sampleNames
    if (!all(grepl("\\.rcc$", sampleNames(object), ignore.case = TRUE))) {
      msg <- c(msg, "'sampleNames' must all have an \".RCC\" file extension")
    }
    # protocolData
    protocolDataColNames <- rownames(.rccMetadata[["protocolData"]])
    if (!all(protocolDataColNames %in% varLabels(protocolData(object)))) {
      msg <-
        c(msg,
          sprintf("'protocolData' must contain columns %s",
                  paste0("\"", protocolDataColNames, "\"", collapse = ", ")))
    }
    # protocolData - FileVersion
    if (!all(protocolData(object)[["FileVersion"]] %in%
             numeric_version(c("1.7", "2.0")))) {
      msg <-
        c(msg, "'protocolData' \"FileVersion\" must all be either 1.7 or 2.0")
    }
    # protocolData - LaneID
    if (!all(protocolData(object)[["LaneID"]] %in% 1L:12L)) {
      msg <-
        c(msg, "'protocolData' \"LaneID\" must all be integers from 1 to 12")
    }
    # protocolData - FovCount
    if (!.validNonNegativeInteger(protocolData(object)[["FovCount"]])) {
      msg <-
        c(msg, "'protocolData' \"FovCount\" must all be non-negative integers")
    }
    # protocolData - FovCounted
    if (!.validNonNegativeInteger(protocolData(object)[["FovCounted"]])) {
      msg <-
        c(msg,
          "'protocolData' \"FovCounted\" must all be non-negative integers")
    } else if (!all(protocolData(object)[["FovCounted"]] <=
                    protocolData(object)[["FovCount"]])) {
      msg <-
        c(msg,
          "'protocolData' \"FovCounted\" must all be less than or equal to \"FovCount\"")
    }
    # protocolData - StagePosition
    if (!all(protocolData(object)[["StagePosition"]] %in% 1L:6L)) {
      msg <-
        c(msg,
          "'protocolData' \"StagePosition\" must all be integers from 1 to 6")
    }
    # protocolData - BindingDensity
    if (!.validNonNegativeNumber(protocolData(object)[["BindingDensity"]])) {
      msg <-
        c(msg,
          "'protocolData' \"BindingDensity\" must all be non-negative numbers")
    }
  }
  if (dim(object)[["Features"]] > 0L) {
    # featureData
    featureDataColNames <- c("BarcodeClass", "GeneName", "Accession",
                             "IsControl", "ControlConc")
    if (!all(featureDataColNames %in% varLabels(featureData(object)))) {
      msg <-
        c(msg,
          sprintf("'featureData' must contain columns %s",
                  paste0("\"", featureDataColNames, "\"", collapse = ", ")))
    } else {
      barcodeClass <- featureData(object)[["BarcodeClass"]]
      isControl <- featureData(object)[["IsControl"]]
      controlConc <- featureData(object)[["ControlConc"]]
      endogen <- which(barcodeClass == "Endogenous")
      posCtrl <- which(barcodeClass == "Positive")
      negCtrl <- which(barcodeClass == "Negative")
      housekp <- which(barcodeClass == "Housekeeping")
      if (length(endogen) > 0L) {
        if (!.allFALSE(isControl[endogen])) {
          msg <-
            c(msg,
              "'featureData': \"IsControl\" must all be FALSE when BarcodeClass == \"Endogenous\"")
        }
        if (!.allNA(controlConc[endogen])) {
          msg <-
            c(msg,
              "'featureData': \"ControlConc\" must all be NA when BarcodeClass == \"Endogenous\"")
        }
      }
      if (length(posCtrl) > 0L) {
        if (!.allTRUE(isControl[posCtrl])) {
          msg <-
            c(msg,
              "'featureData': \"IsControl\" must all be TRUE when BarcodeClass == \"Positive\"")
        }
        if (!.validPositiveNumber(controlConc[posCtrl])) {
          msg <-
            c(msg,
              "'featureData': \"ControlConc\" must be positive numbers when BarcodeClass == \"Positive\"")
        }
      }
      if (length(negCtrl) > 0L) {
        if (!.allTRUE(isControl[negCtrl])) {
          msg <-
            c(msg,
              "'featureData': \"IsControl\" must all be TRUE when BarcodeClass == \"Negative\"")
        }
        if (!.allZero(controlConc[negCtrl])) {
          msg <-
            c(msg,
              "'featureData': \"ControlConc\" must all be zero when BarcodeClass == \"Negative\"")
        }
      }
      if (length(housekp) > 0L) {
        if (!.allFALSE(isControl[housekp])) {
          msg <-
            c(msg,
              "'featureData': \"IsControl\" must all be FALSE when BarcodeClass == \"Housekeeping\"")
        }
        if (!.allNA(controlConc[housekp])) {
          msg <-
            c(msg,
              "'featureData': \"ControlConc\" must all be NA when BarcodeClass == \"Housekeeping\"")
        }
      }
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
    if (!.validNonNegativeInteger(exprs(object))) {
      msg <- c(msg, "'exprs' must be a non-negative integer matrix")
    }
  }
  if (any(duplicated(c(fvarLabels(object), svarLabels(object),
                       assayDataElementNames(object))))) {
    msg <-
      c(msg,
        "'fvarLabels', 'svarLabels', and 'assayDataElementNames' must be unique")
  }
  if (length(signatureWeights(object)) > 0L) {
    numGenes <- lengths(signatureWeights(object))
    if (is.null(names(numGenes)) || any(nchar(names(numGenes)) == 0L)) {
      msg <- c(msg, "'signatureWeights' must be a named NumericList")
    }
    if (any(numGenes == 0L)) {
      msg <- c(msg, "'signatureWeights' vectors must be non-empty")
    } else {
      genes <- names(unlist(unname(signatureWeights(object))))
      if (is.null(genes) || any(nchar(genes) == 0L)) {
        msg <- c(msg, "'signatureWeights' vectors must be named")
      } else if(!all(unique(genes) %in% featureData(object)[["GeneName"]])) {
        msg <-
          c(msg,
            "'signatureWeights' vectors must be named with values from 'featureData' \"GeneName\"")
      }
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
         signatureWeights = NumericList(),
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
         signatureWeights = NumericList(),
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              signatureWeights = signatureWeights, ...)
})

setMethod("NanoStringRccSet", "environment",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         ...)
{
  new2("NanoStringRccSet",
       assayData = assayData,
       phenoData = phenoData,
       featureData = featureData,
       experimentData = experimentData,
       annotation = annotation,
       protocolData = protocolData,
       signatureWeights = signatureWeights,
       ...)
})

setMethod("NanoStringRccSet", "matrix",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              signatureWeights = signatureWeights, ...)
})

setMethod("NanoStringRccSet", "NanoStringRccSet",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         ...)
{
  callGeneric(assayData = copyEnv(assayData(assayData)),
              phenoData = Biobase::phenoData(assayData),
              featureData = Biobase::featureData(assayData),
              experimentData = Biobase::experimentData(assayData),
              annotation = Biobase::annotation(assayData),
              protocolData = Biobase::protocolData(assayData),
              signatureWeights = signatureWeights,
              ...)
})
