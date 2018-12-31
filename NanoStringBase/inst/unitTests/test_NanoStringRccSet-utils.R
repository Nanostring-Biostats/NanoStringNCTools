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
                                       IsControl = c(FALSE, TRUE, TRUE, FALSE),
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
                            NanoStringBase:::.rccMetadata[["protocolData"]],
                            dimLabels = c("sampleNames", "sampleColumns")),
       signatureWeights =
         list(x = c(a = 1), y = c(b = 1/3, d = 2/3), z = c(a = 2, c = 4)))

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
  checkEquals(cbind(N = 3,
                    Mean = structure(4:7, names = letters[1:4]),
                    SD = 4,
                    Skewness = 0,
                    Kurtosis = NA_real_,
                    GM_LCL95 = c(NA_real_, 0.2108001028, 0.6394432507, 1.1919981586),
                    GeomMean = c(NA_real_, 3.556893304, 4.932424149, 6.135792440),
                    GM_UCL95 = c(NA_real_, 60.01652661, 38.04686023, 31.58389851),
                    Min = 0:3,
                    Q1 = 2:5,
                    Median = 4:7,
                    Q3 = 6:9,
                    Max = 8:11,
                    MAD = 5.9304,
                    MedPolEff = c(-1.5, -0.5, 0.5, 1.5)),
              summary(rcc2, 1L))
  checkEquals(cbind(N = 3,
                    Mean = structure(5:8, names = letters[1:4]),
                    SD = 4,
                    Skewness = 0,
                    Kurtosis = NA_real_,
                    GM_LCL95 = c(0.2108001028, 0.6394432507, 1.1919981586, 1.8284870354),
                    GeomMean = c(3.556893304, 4.932424149, 6.135792440, 7.268482371),
                    GM_UCL95 = c(60.01652661, 38.04686023, 31.58389851, 28.89319692),
                    Min = 1:4,
                    Q1 = 3:6,
                    Median = 5:8,
                    Q3 = 7:10,
                    Max = 9:12,
                    MAD = 5.9304,
                    MedPolEff = c(-1.5, -0.5, 0.5, 1.5)),
              summary(rcc2, 1L, elt = "exprsp1"))

  # Marginal summaries by Sample
  checkEquals(cbind(N = 4,
                    Mean = structure(c(1.5, 5.5, 9.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    GM_LCL95 = c(NA_real_, 3.668188416, 7.584766676),
                    GeomMean = c(NA_real_, 5.383563271, 9.433683366),
                    GM_UCL95 = c(NA_real_, 7.901108178, 11.733305143),
                    Min = c(0, 4, 8),
                    Q1 = c(0.75, 4.75, 8.75),
                    Median = c(1.5, 5.5, 9.5),
                    Q3 = c(2.25, 6.25, 10.25),
                    Max = c(3, 7, 11),
                    MAD = 1.4826,
                    MedPolEff = c(-4, 0, 4)),
              summary(rcc2, 2L))
  checkEquals(cbind(N = 4,
                    Mean = structure(c(2.5, 6.5, 10.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    GM_LCL95 = c(0.8503745621, 4.6391603807, 8.5728555505),
                    GeomMean = c(2.213363839, 6.402171746, 10.440086817),
                    GM_UCL95 = c(5.760966642, 8.835176993, 12.714014847),
                    Min = c(1, 5, 9),
                    Q1 = c(1.75, 5.75, 9.75),
                    Median = c(2.5, 6.5, 10.5),
                    Q3 = c(3.25, 7.25, 11.25),
                    Max = c(4, 8, 12),
                    MAD = 1.4826,
                    MedPolEff = c(-4, 0, 4)),
              summary(rcc2, 2L, elt = "exprsp1"))
}

test_NanoStringRccSet_utils_summary_GROUP <- function() {
  rcc2 <- transform(rcc, exprsp1 = exprs + 1L)

  # Marginal summaries by Feature
  checkEquals(list(A =
                     cbind(N = 2,
                           Mean = structure(2:5, names = letters[1:4]),
                           SD = sqrt(8),
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = c(NA_real_, 8.106941001e-05, 3.223965976e-03, 2.105307130e-02),
                           GeomMean = c(NA_real_, 2.236067977, 3.464101615, 4.582575695),
                           GM_UCL95 = c(NA_real_, 61675.5444443, 3722.1236485, 997.4791658),
                           Min = 0:3,
                           Q1 = 1:4,
                           Median = 2:5,
                           Q3 = 3:6,
                           Max = 4:7,
                           MAD = 2.9652,
                           MedPolEff = c(-1.5, -0.5, 0.5, 1.5)),
                   B =
                     cbind(N = 1,
                           Mean = structure(8:11, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = 8:11,
                           GM_UCL95 = NA_real_,
                           Min = 8:11,
                           Q1 = 8:11,
                           Median = 8:11,
                           Q3 = 8:11,
                           Max = 8:11,
                           MAD = 0,
                           MedPolEff = c(-1.5, -0.5, 0.5, 1.5))),
              summary(rcc2, 1L, "Treatment"))
  checkEquals(list("1" =
                     cbind(N = 1,
                           Mean = structure(1:4, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = 1:4,
                           GM_UCL95 = NA_real_,
                           Min = 1:4,
                           Q1 = 1:4,
                           Median = 1:4,
                           Q3 = 1:4,
                           Max = 1:4,
                           MAD = 0,
                           MedPolEff = c(-1.5, -0.5, 0.5, 1.5)),
                   "2" =
                     cbind(N = 1,
                           Mean = structure(5:8, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = 5:8,
                           GM_UCL95 = NA_real_,
                           Min = 5:8,
                           Q1 = 5:8,
                           Median = 5:8,
                           Q3 = 5:8,
                           Max = 5:8,
                           MAD = 0,
                           MedPolEff = c(-1.5, -0.5, 0.5, 1.5)),
                   "3" =
                     cbind(N = 1,
                           Mean = structure(9:12, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = 9:12,
                           GM_UCL95 = NA_real_,
                           Min = 9:12,
                           Q1 = 9:12,
                           Median = 9:12,
                           Q3 = 9:12,
                           Max = 9:12,
                           MAD = 0,
                           MedPolEff = c(-1.5, -0.5, 0.5, 1.5))),
              summary(rcc2, 1L, "LaneID", elt = "exprsp1"))

  # Marginal summaries by Sample
  checkEquals(list(Endogenous =
                     cbind(N = 1,
                           Mean = structure(c(1, 5, 9), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = c(1, 5, 9),
                           GM_UCL95 = NA_real_,
                           Min = c(1, 5, 9),
                           Q1 = c(1, 5, 9),
                           Median = c(1, 5, 9),
                           Q3 = c(1, 5, 9),
                           Max = c(1, 5, 9),
                           MAD = 0,
                           MedPolEff = c(-4, 0, 4)),
                   Housekeeping =
                     cbind(N = 1,
                           Mean = structure(c(4, 8, 12), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = c(4, 8, 12),
                           GM_UCL95 = NA_real_,
                           Min = c(4, 8, 12),
                           Q1 = c(4, 8, 12),
                           Median = c(4, 8, 12),
                           Q3 = c(4, 8, 12),
                           Max = c(4, 8, 12),
                           MAD = 0,
                           MedPolEff = c(-4, 0, 4)),
                   Negative =
                     cbind(N = 1,
                           Mean = structure(c(3, 7, 11), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = c(3, 7, 11),
                           GM_UCL95 = NA_real_,
                           Min = c(3, 7, 11),
                           Q1 = c(3, 7, 11),
                           Median = c(3, 7, 11),
                           Q3 = c(3, 7, 11),
                           Max = c(3, 7, 11),
                           MAD = 0,
                           MedPolEff = c(-4, 0, 4)),
                   Positive =
                     cbind(N = 1,
                           Mean = structure(c(2, 6, 10), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           GM_LCL95 = NA_real_,
                           GeomMean = c(2, 6, 10),
                           GM_UCL95 = NA_real_,
                           Min = c(2, 6, 10),
                           Q1 = c(2, 6, 10),
                           Median = c(2, 6, 10),
                           Q3 = c(2, 6, 10),
                           Max = c(2, 6, 10),
                           MAD = 0,
                           MedPolEff = c(-4, 0, 4))),
              summary(rcc2, 2L, "BarcodeClass", elt = "exprsp1"))
}

# Subsetting
test_NanoStringRccSet_utils_subset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous", ],
              subset(rcc, BarcodeClass == "Endogenous"))
  checkEquals(NumericList(x = c(a = 1), compress = FALSE),
              signatureWeights(subset(rcc, BarcodeClass == "Endogenous")))

  checkEquals(rcc[, phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, select = Treatment == "A"))
  checkEquals(signatureWeights(rcc),
              signatureWeights(subset(rcc, select = Treatment == "A")))

  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous",
                  phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, BarcodeClass == "Endogenous", Treatment == "A"))
  checkEquals(NumericList(x = c(a = 1), compress = FALSE),
              signatureWeights(subset(rcc, BarcodeClass == "Endogenous",
                                      Treatment == "A")))
}

test_NanoStringRccSet_utils_endogenousSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Endogenous", ],
              endogenousSubset(rcc))
}

test_NanoStringRccSet_utils_housekeepingSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Housekeeping", ],
              housekeepingSubset(rcc))
}

test_NanoStringRccSet_utils_negativeControlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Negative", ],
              negativeControlSubset(rcc))
}

test_NanoStringRccSet_utils_positiveControlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Positive", ],
              positiveControlSubset(rcc))

  checkEquals(rcc[featureData(rcc)[["BarcodeClass"]] == "Positive" &
                  featureData(rcc)[["ControlConc"]] >= 0.5, ],
              positiveControlSubset(rcc, subset = ControlConc >= 0.5))
}

test_NanoStringRccSet_utils_controlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["IsControl"]], ], controlSubset(rcc))
}

test_NanoStringRccSet_utils_nonControlSubset <- function() {
  checkEquals(rcc[!featureData(rcc)[["IsControl"]], ], nonControlSubset(rcc))
}

test_NanoStringRccSet_utils_signatureSubset <- function() {
  x <- rcc
  signatureWeights(x) <- signatureWeights(x)[1L]
  checkEquals(rcc[featureData(rcc)[["GeneName"]] == "a", ],
              signatureSubset(x))
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
  nms <- sort(c(assayDataElementNames(rcc), fvarLabels(rcc), svarLabels(rcc),
                "signatureWeights", "design"))
  checkIdentical(nms, with(rcc, ls()))

  # calculate means across Features
  checkIdentical(esApply(rcc, 1L, mean),
                 with(rcc, apply(exprs, 1L, mean)))

  # calculate means across Samples
  checkIdentical(esApply(rcc, 2L, mean),
                 with(rcc, apply(exprs, 2L, mean)))
}
