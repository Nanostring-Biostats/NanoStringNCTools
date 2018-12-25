rcc <-
  NanoStringRccSet(assayData =
         matrix(0:11, 4L, 3L,
                dimnames = list(letters[1:4], sprintf("%s.RCC", LETTERS[1:3]))),
       phenoData =
         AnnotatedDataFrame(data.frame(Treatment = c("A", "A", "B"),
                                       Age = c(58L, 42L, 27L),
                                       row.names = sprintf("%s.RCC", LETTERS[1:3]),
                                       stringsAsFactors = FALSE),
                            dimLabels = c("sampleNames", "sampleColumns")),
       featureData =
         AnnotatedDataFrame(data.frame(BarcodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
                                       GeneName = letters[1:4],
                                       Accession = letters[1:4],
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
                            NanoStringBase:::.rccMetadata[["protocolData"]],
                            dimLabels = c("sampleNames", "sampleColumns")))

# Accessing
test_NanoStringRccSet_utils_sData <- function() {
  checkIdentical(cbind(pData(rcc), pData(protocolData(rcc))), sData(rcc))
}

test_NanoStringRccSet_utils_svarLabels <- function() {
  checkIdentical(c(varLabels(rcc), varLabels(protocolData(rcc))), svarLabels(rcc))
}

# Sumarizing
test_NanoStringRccSet_utils_summary <- function() {
  # Marginal summaries by Feature
  checkIdentical(cbind(Mean = structure(4:7, names = letters[1:4]),
                       SD = 4,
                       Skewness = 0,
                       Kurtosis = NaN,
                       Min = 0:3,
                       Q1 = 2:5,
                       Median = 4:7,
                       Q3 = 6:9,
                       Max = 8:11,
                       N = 3,
                       NMiss = 0),
                 summary(rcc, 1L))

  # Marginal summaries by Sample
  checkEquals(cbind(Mean = structure(c(1.5, 5.5, 9.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    Min = c(0, 4, 8),
                    Q1 = c(0.75, 4.75, 8.75),
                    Median = c(1.5, 5.5, 9.5),
                    Q3 = c(2.25, 6.25, 10.25),
                    Max = c(3, 7, 11),
                    N = 4,
                    NMiss = 0),
              summary(rcc, 2L))
}

# Subsetting
test_NanoStringRccSet_utils_controlSet <- function() {
  checkException(controlSet(rcc))

  x <- rcc
  featureData(x)[["IsControl"]] <- c(TRUE, FALSE, TRUE, FALSE)
  checkEquals(x[featureData(x)[["IsControl"]], ], controlSet(x))
}

test_NanoStringRccSet_utils_nonControlSet <- function() {
  checkException(nonControlSet(rcc))

  x <- rcc
  featureData(x)[["IsControl"]] <- c(TRUE, FALSE, TRUE, FALSE)
  checkEquals(x[!featureData(x)[["IsControl"]], ], nonControlSet(x))
}

test_NanoStringRccSet_utils_endogenousSet <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous", ],
              endogenousSet(rcc))
}

test_NanoStringRccSet_utils_housekeepingSet <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Housekeeping", ],
              housekeepingSet(rcc))
}

test_NanoStringRccSet_utils_negativeControlSet <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Negative", ],
              negativeControlSet(rcc))
}

test_NanoStringRccSet_utils_positiveControlSet <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Positive", ],
              positiveControlSet(rcc))
}

# transforming
test_NanoStringRccSet_utils_sweep <- function() {
  # scale across Features
  rcc2 <- sweep(rcc, 1L, rowMeans(exprs(rcc)), toElt = "centered")
  rcc2 <- sweep(rcc2, 1L, apply(exprs(rcc2), 1L, sd), FUN = "/",
                fromElt = "centered", toElt = "scaled")
  target <- t(scale(t(exprs(rcc))))
  attr(target, "scaled:center") <- attr(target, "scaled:scale") <- NULL
  checkTrue(validObject(rcc2))
  checkIdentical(target, assayDataElement(rcc2, "scaled"))

  # scale across Samples
  rcc2 <- sweep(rcc, 2L, colMeans(exprs(rcc)), toElt = "centered")
  rcc2 <- sweep(rcc2, 2L, apply(exprs(rcc2), 2L, sd), FUN = "/",
                fromElt = "centered", toElt = "scaled")
  target <- scale(exprs(rcc))
  attr(target, "scaled:center") <- attr(target, "scaled:scale") <- NULL
  checkTrue(validObject(rcc2))
  checkIdentical(target, assayDataElement(rcc2, "scaled"))
}

test_NanoStringRccSet_utils_transform <- function() {
  rcc2 <- transform(rcc,
                    exprs_thresh = pmax(exprs - 2, 0),
                    log1p_exprs_thresh = log1p(exprs_thresh))
  checkTrue(validObject(rcc2))
  checkIdentical(pmax(exprs(rcc) - 2, 0),
                 assayDataElement(rcc2, "exprs_thresh"))
  checkIdentical(log1p(pmax(exprs(rcc) - 2, 0)),
                 assayDataElement(rcc2, "log1p_exprs_thresh"))
}

# evaluating
test_NanoStringRccSet_utils_with <- function() {
  nms <- sort(c(assayDataElementNames(rcc), fvarLabels(rcc), svarLabels(rcc)))
  checkIdentical(nms, with(rcc, ls()))

  # calculate means across Features
  checkIdentical(esApply(rcc, 1L, mean),
                 with(rcc, apply(exprs, 1L, mean)))

  # calculate means across Samples
  checkIdentical(esApply(rcc, 2L, mean),
                 with(rcc, apply(exprs, 2L, mean)))
}
