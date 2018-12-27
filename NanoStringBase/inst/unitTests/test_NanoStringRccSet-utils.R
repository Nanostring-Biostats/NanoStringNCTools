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
  rcc2 <- transform(rcc, exprsp1 = exprs + 1L)

  # Marginal summaries by Feature
  checkIdentical(cbind(N = 3,
                       NMiss = 0,
                       Mean = structure(4:7, names = letters[1:4]),
                       SD = 4,
                       Skewness = 0,
                       Kurtosis = NaN,
                       Min = 0:3,
                       Q1 = 2:5,
                       Median = 4:7,
                       Q3 = 6:9,
                       Max = 8:11,
                       MAD = 5.9304,
                       MedPolEffect = c(-1.5, -0.5, 0.5, 1.5)),
                 summary(rcc2, 1L))
  checkIdentical(cbind(N = 3,
                       NMiss = 0,
                       Mean = structure(5:8, names = letters[1:4]),
                       SD = 4,
                       Skewness = 0,
                       Kurtosis = NaN,
                       Min = 1:4,
                       Q1 = 3:6,
                       Median = 5:8,
                       Q3 = 7:10,
                       Max = 9:12,
                       MAD = 5.9304,
                       MedPolEffect = c(-1.5, -0.5, 0.5, 1.5)),
                 summary(rcc2, 1L, elt = "exprsp1"))

  # Marginal summaries by Sample
  checkEquals(cbind(N = 4,
                    NMiss = 0,
                    Mean = structure(c(1.5, 5.5, 9.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    Min = c(0, 4, 8),
                    Q1 = c(0.75, 4.75, 8.75),
                    Median = c(1.5, 5.5, 9.5),
                    Q3 = c(2.25, 6.25, 10.25),
                    Max = c(3, 7, 11),
                    MAD = 1.4826,
                    MedPolEffect = c(-4, 0, 4)),
              summary(rcc2, 2L))
  checkEquals(cbind(N = 4,
                    NMiss = 0,
                    Mean = structure(c(2.5, 6.5, 10.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    Min = c(1, 5, 9),
                    Q1 = c(1.75, 5.75, 9.75),
                    Median = c(2.5, 6.5, 10.5),
                    Q3 = c(3.25, 7.25, 11.25),
                    Max = c(4, 8, 12),
                    MAD = 1.4826,
                    MedPolEffect = c(-4, 0, 4)),
              summary(rcc2, 2L, elt = "exprsp1"))
}

# Subsetting
test_NanoStringRccSet_utils_subset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous", ],
              subset(rcc, BarcodeClass == "Endogenous"))

  checkEquals(rcc[, phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, select = Treatment == "A"))

  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous",
                  phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, BarcodeClass == "Endogenous", Treatment == "A"))
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

# looping
test_NanoStringRccSet_utils_esApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(rcc2), 1L, mean), esApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(rcc2), 2L, mean), esApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(rcc2, "log1p_exprs"), 1L, mean),
                 esApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(rcc2, "log1p_exprs"), 2L, mean),
                 esApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

test_NanoStringRccSet_utils_endogenousApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(endogenousSet(rcc2)), 1L, mean),
                 endogenousApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(endogenousSet(rcc2)), 2L, mean),
                 endogenousApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(endogenousSet(rcc2), "log1p_exprs"), 1L, mean),
                 endogenousApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(endogenousSet(rcc2), "log1p_exprs"), 2L, mean),
                 endogenousApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

test_NanoStringRccSet_utils_housekeepingApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(housekeepingSet(rcc2)), 1L, mean),
                 housekeepingApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(housekeepingSet(rcc2)), 2L, mean),
                 housekeepingApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(housekeepingSet(rcc2), "log1p_exprs"), 1L, mean),
                 housekeepingApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(housekeepingSet(rcc2), "log1p_exprs"), 2L, mean),
                 housekeepingApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

test_NanoStringRccSet_utils_negativeControlApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(negativeControlSet(rcc2)), 1L, mean),
                 negativeControlApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(negativeControlSet(rcc2)), 2L, mean),
                 negativeControlApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(negativeControlSet(rcc2), "log1p_exprs"), 1L, mean),
                 negativeControlApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(negativeControlSet(rcc2), "log1p_exprs"), 2L, mean),
                 negativeControlApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

test_NanoStringRccSet_utils_positiveControlApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(positiveControlSet(rcc2)), 1L, mean),
                 positiveControlApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(positiveControlSet(rcc2)), 2L, mean),
                 positiveControlApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(positiveControlSet(rcc2), "log1p_exprs"), 1L, mean),
                 positiveControlApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(positiveControlSet(rcc2), "log1p_exprs"), 2L, mean),
                 positiveControlApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

# transforming
test_NanoStringRccSet_utils_transform <- function() {
  rcc2 <- transform(rcc,
                    exprs_thresh = pmax(exprs - 2, 0),
                    log1p_exprs_thresh = log1p(exprs_thresh))
  checkTrue(validObject(rcc2))
  checkIdentical(pmax(exprs(rcc) - 2, 0),
                 assayDataElement(rcc2, "exprs_thresh"))
  checkIdentical(log1p(pmax(exprs(rcc) - 2, 0)),
                 assayDataElement(rcc2, "log1p_exprs_thresh"))
  checkIdentical(list(exprs_thresh = substitute(pmax(exprs - 2, 0)),
                      log1p_exprs_thresh = substitute(log1p(exprs_thresh))),
                 preproc(rcc2))
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
