test_NanoStringRccSet_constructor_empty <- function() {
  checkTrue(validObject(NanoStringRccSet()))
}

rcc <-
  list(assayData =
         matrix(0:11, 4L, 3L,
                dimnames = list(letters[1:4], sprintf("%s.RCC", LETTERS[1:3]))),
       phenoData =
         AnnotatedDataFrame(data.frame(Treatment = c("A", "A", "B"),
                                       Age = c(58L, 42L, 27L),
                                       row.names = sprintf("%s.RCC", LETTERS[1:3]),
                                       stringsAsFactors = FALSE),
                            dimLabels = c("sampleNames", "sampleColumns")),
       featureData =
         AnnotatedDataFrame(data.frame(CodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
                                       GeneName = letters[1:4],
                                       Accession = letters[1:4],
                                       IsControl = c(FALSE, TRUE, TRUE, TRUE),
                                       ControlConc = c(NA_real_, 0.125, 0, NA_real_),
                                       row.names = letters[1:4],
                                       stringsAsFactors = FALSE),
                            dimLabels = c("featureNames", "featureColumns")),
       annotation = "rlffile",
       protocolData =
         AnnotatedDataFrame(data.frame(FileVersion = numeric_version(rep("1.7", 3L)),
                                       SoftwareVersion = numeric_version(rep("4.0.0.3", 3L)),
                                       SystemType = rep("Gen2", 3L),
                                       SampleID = letters[1:3],
                                       SampleOwner = rep("", 3L),
                                       SampleComments = rep("DNA-RNA-Protein", 3L),
                                       SampleDate = as.Date(rep("1999-12-31", 3L)),
                                       SystemAPF = rep("n6_vDV1", 3L),
                                       AssayType = rep(NA_character_, 3L),
                                       LaneID = 1:3,
                                       FovCount = rep(280L, 3L),
                                       FovCounted = 1:3,
                                       ScannerID = rep("a", 3L),
                                       StagePosition = 1:3,
                                       BindingDensity = c(0.75, 1, 1.25),
                                       CartridgeID = letters[1:3],
                                       CartridgeBarcode = rep("", 3L),
                                       row.names = sprintf("%s.RCC", LETTERS[1:3]),
                                       stringsAsFactors = FALSE),
                            NanoStringNCTools:::.rccMetadata[["protocolData"]],
                            dimLabels = c("sampleNames", "sampleColumns")))

test_NanoStringRccSet_constructor_simple <- function() {
  checkTrue(validObject(NanoStringRccSet(rcc$assayData,
                                         phenoData = rcc$phenoData,
                                         featureData = rcc$featureData,
                                         annotation = rcc$annotation,
                                         protocolData = rcc$protocolData)))
}

test_NanoStringRccSet_exception_sample_name <- function() {
  x <- rcc$assayData
  colnames(x) <- letters[1:3]
  y <- rcc$protocolData
  rownames(y) <- letters[1:3]
  checkException(validObject(NanoStringRccSet(x,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = y)))
}

test_NanoStringRccSet_exception_featureData <- function() {
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = rcc$featureData[,1:2],
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))
}

test_NanoStringRccSet_exception_featureData_IsControl <- function() {
  x <- rcc$featureData
  x[["IsControl"]] <- TRUE
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))

  x <- rcc$featureData
  x[["IsControl"]] <- FALSE
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))

  x <- rcc$featureData
  x[["IsControl"]] <- NA
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))
}

test_NanoStringRccSet_exception_featureData_ControlConc <- function() {
  x <- rcc$featureData
  x[["ControlConc"]] <- 0.125
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))

  x <- rcc$featureData
  x[["ControlConc"]] <- 0
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))

  x <- rcc$featureData
  x[["ControlConc"]] <- NA_real_
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = x,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))
}

test_NanoStringRccSet_exception_annotation <- function() {
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              phenoData = rcc$phenoData,
                                              featureData = rcc$featureData,
                                              protocolData = rcc$protocolData)))
}

test_NanoStringRccSet_exception_protocolData <- function() {
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData[,1:2])))
}

test_NanoStringRccSet_exception_protocolData_FileVersion <- function() {
  x <- rcc$protocolData
  x[["FileVersion"]] <- numeric_version("1.6")
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_protocolData_LaneID <- function() {
  x <- rcc$protocolData
  x[["LaneID"]] <- 0L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
  x[["LaneID"]] <- 13L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_protocolData_FovCount <- function() {
  x <- rcc$protocolData
  x[["FovCount"]] <- -1L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_protocolData_FovCounted <- function() {
  x <- rcc$protocolData
  x[["FovCounted"]] <- -1L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_protocolData_StagePosition <- function() {
  x <- rcc$protocolData
  x[["StagePosition"]] <- 0L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
  x[["StagePosition"]] <- 7L
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_protocolData_BindingDensity <- function() {
  x <- rcc$protocolData
  x[["BindingDensity"]] <- -1.2
  checkException(validObject(NanoStringRccSet(rcc$assayData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = x)))
}

test_NanoStringRccSet_exception_assayData <- function() {
  # exprs matrix
  x <- rcc$assayData
  x[] <- -1L
  checkException(validObject(NanoStringRccSet(x,
                                              phenoData = rcc$phenoData,
                                              featureData = rcc$featureData,
                                              annotation = rcc$annotation,
                                              protocolData = rcc$protocolData)))

  # additional assayData matrices
  x <- NanoStringRccSet(rcc$assayData,
                        phenoData = rcc$phenoData,
                        featureData = rcc$featureData,
                        annotation = rcc$annotation,
                        protocolData = rcc$protocolData)
  checkTrue(validObject(x))
  assayDataElement(x, "shifted") <- rcc$assayData - 1L
  checkTrue(validObject(x))
}

test_NanoStringRccSet_exception_duplicate_names <- function() {
  x <- NanoStringRccSet(rcc$assayData,
                        phenoData = rcc$phenoData,
                        featureData = rcc$featureData,
                        annotation = rcc$annotation,
                        protocolData = rcc$protocolData)

  # Duplicate fvarLabels
  y <- x
  fData(y)[["exprs"]] <- 1L
  checkException(validObject(y))

  y <- x
  fData(y)[["SampleID"]] <- 1L
  checkException(validObject(y))

  y <- x
  fData(y)[["signatures"]] <- 1L
  checkException(validObject(y))

  y <- x
  fData(y)[["design"]] <- 1L
  checkException(validObject(y))

  # Duplicate pvarLabels
  y <- x
  pData(y)[["exprs"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(y)[["CodeClass"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(y)[["signatures"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(y)[["design"]] <- 1L
  checkException(validObject(y))

  # Duplicate pvarLables(protocolData)
  y <- x
  pData(protocolData(y))[["exprs"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(protocolData(y))[["CodeClass"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(protocolData(y))[["signatures"]] <- 1L
  checkException(validObject(y))

  y <- x
  pData(protocolData(y))[["design"]] <- 1L
  checkException(validObject(y))
}


test_NanoStringRccSet_exception_signatures_unnamed <- function() {
  checkException(NanoStringRccSet(rcc$assayData,
                                  phenoData = rcc$phenoData,
                                  featureData = rcc$featureData,
                                  annotation = rcc$annotation,
                                  protocolData = rcc$protocolData,
                                  signatures =
                                    SignatureSet(weights =
                                                   list(c(b = 1, d = 3)))))
}

test_NanoStringRccSet_exception_signatures_empty <- function() {
  checkException(NanoStringRccSet(rcc$assayData,
                                  phenoData = rcc$phenoData,
                                  featureData = rcc$featureData,
                                  annotation = rcc$annotation,
                                  protocolData = rcc$protocolData,
                                  signatures =
                                    SignatureSet(weights =
                                                   list(a = numeric()))))
}

test_NanoStringRccSet_exception_signatures_GeneName <- function() {
  checkException(NanoStringRccSet(rcc$assayData,
                                  phenoData = rcc$phenoData,
                                  featureData = rcc$featureData,
                                  annotation = rcc$annotation,
                                  protocolData = rcc$protocolData,
                                  signatures =
                                    SignatureSet(weights =
                                                   list(a = c(x = 1, y = 3)))))
}

test_NanoStringRccSet_valid_design <- function() {
  x <- NanoStringRccSet(rcc$assayData,
                        phenoData = rcc$phenoData,
                        featureData = rcc$featureData,
                        annotation = rcc$annotation,
                        protocolData = rcc$protocolData,
                        design = ~ Treatment)
  checkTrue(validObject(x))
  checkIdentical(design(x), ~ Treatment)

  design(x) <- ~ Age + Treatment
  checkTrue(validObject(x))
  checkIdentical(design(x), ~ Age + Treatment)
}

test_NanoStringRccSet_exception_design <- function() {
  checkException(NanoStringRccSet(rcc$assayData,
                                  phenoData = rcc$phenoData,
                                  featureData = rcc$featureData,
                                  annotation = rcc$annotation,
                                  protocolData = rcc$protocolData,
                                  design = Age ~ Treatment))

  checkException(NanoStringRccSet(rcc$assayData,
                                  phenoData = rcc$phenoData,
                                  featureData = rcc$featureData,
                                  annotation = rcc$annotation,
                                  protocolData = rcc$protocolData,
                                  design = ~ CodeClass))
}
