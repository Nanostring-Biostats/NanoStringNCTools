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
  checkEquals(cbind(GeomMean = c(2.519842100, 3.556893304, 4.932424149, 6.135792440),
                    SizeFactor = c(0.74300232, 0.92875290, 1.11450348, 1.30025406),
                    MedPolSF = c(0.7302967433, 0.9128709292, 1.0954451150, 1.2780193008),
                    Mean = structure(4:7, names = letters[1:4]),
                    SD = 4,
                    Skewness = 0,
                    Kurtosis = NA_real_,
                    Min = 0:3,
                    Q1 = 2:5,
                    Median = 4:7,
                    Q3 = 6:9,
                    Max = 8:11,
                    MAD = 5.9304),
              summary(rcc2, 1L))
  checkEquals(cbind(GeomMean = c(3.556893304, 4.932424149, 6.135792440, 7.268482371),
                    SizeFactor = c(0.7809849842, 0.9371819811, 1.0933789779, 1.2495759748),
                    MedPolSF = c(0.7715167498, 0.9258200998, 1.0801234497, 1.2344267997),
                    Mean = structure(5:8, names = letters[1:4]),
                    SD = 4,
                    Skewness = 0,
                    Kurtosis = NA_real_,
                    Min = 1:4,
                    Q1 = 3:6,
                    Median = 5:8,
                    Q3 = 7:10,
                    Max = 9:12,
                    MAD = 5.9304),
              summary(rcc2, 1L, elt = "exprsp1"))

  # Marginal summaries by Sample
  checkEquals(cbind(GeomMean = c(1.316074013, 5.383563271, 9.433683366),
                    SizeFactor = c(0.3376364857, 1.3076604860, 2.2649344008),
                    MedPolSF = c(0.2581988897, 1, 1.7320508076),
                    Mean = structure(c(1.5, 5.5, 9.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    Min = c(0, 4, 8),
                    Q1 = c(0.75, 4.75, 8.75),
                    Median = c(1.5, 5.5, 9.5),
                    Q3 = c(2.25, 6.25, 10.25),
                    Max = c(3, 7, 11),
                    MAD = 1.4826),
              summary(rcc2, 2L))
  checkEquals(cbind(GeomMean = c(2.213363839, 6.402171746, 10.440086817),
                    SizeFactor = c(0.4452563148, 1.1780374787, 1.9064736403),
                    MedPolSF = c(0.377964473, 1, 1.618347187),
                    Mean = structure(c(2.5, 6.5, 10.5), names = sampleNames(rcc)),
                    SD = sqrt(15/9),
                    Skewness = 0,
                    Kurtosis = -1.2,
                    Min = c(1, 5, 9),
                    Q1 = c(1.75, 5.75, 9.75),
                    Median = c(2.5, 6.5, 10.5),
                    Q3 = c(3.25, 7.25, 11.25),
                    Max = c(4, 8, 12),
                    MAD = 1.4826),
              summary(rcc2, 2L, elt = "exprsp1"))
}

test_NanoStringRccSet_utils_summary_GROUP <- function() {
  rcc2 <- transform(rcc, exprsp1 = exprs + 1L)

  # Marginal summaries by Feature
  checkEquals(list(A =
                     cbind(GeomMean = c(1.414213562, 2.236067977, 3.464101615, 4.582575695),
                           SizeFactor = c(0.5313001399, 0.8400592816, 1.3014142429, 1.7216092197),
                           MedPolSF = c(0.5081327482, 0.8034284189, 1.2446659546, 1.6465382906),
                           Mean = structure(2:5, names = letters[1:4]),
                           SD = sqrt(8),
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = 0:3,
                           Q1 = 1:4,
                           Median = 2:5,
                           Q3 = 3:6,
                           Max = 4:7,
                           MAD = 2.9652),
                   B =
                     cbind(GeomMean = 8:11,
                           SizeFactor = c(0.8480250704, 0.9540282042, 1.0600313380, 1.1660344717),
                           MedPolSF = c(0.8432740427, 0.9486832981, 1.0540925534, 1.1595018087),
                           Mean = structure(8:11, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = 8:11,
                           Q1 = 8:11,
                           Median = 8:11,
                           Q3 = 8:11,
                           Max = 8:11,
                           MAD = 0)),
              summary(rcc2, 1L, "Treatment"))
  checkEquals(list("1" =
                     cbind(GeomMean = 1:4,
                           SizeFactor = c(0.4518010018, 0.9036020036, 1.3554030054, 1.8072040072),
                           MedPolSF = c(0.4082482905, 0.8164965809, 1.2247448714, 1.6329931619),
                           Mean = structure(1:4, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = 1:4,
                           Q1 = 1:4,
                           Median = 1:4,
                           Q3 = 1:4,
                           Max = 1:4,
                           MAD = 0),
                   "2" =
                     cbind(GeomMean = 5:8,
                           SizeFactor = c(0.7809849842, 0.9371819811, 1.0933789779, 1.2495759748),
                           MedPolSF = c(0.7715167498, 0.9258200998, 1.0801234497, 1.2344267997),
                           Mean = structure(5:8, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = 5:8,
                           Q1 = 5:8,
                           Median = 5:8,
                           Q3 = 5:8,
                           Max = 5:8,
                           MAD = 0),
                   "3" =
                     cbind(GeomMean = 9:12,
                           SizeFactor = c(0.8620617968, 0.9578464409, 1.0536310849, 1.1494157290),
                           MedPolSF = c(0.8581163303, 0.9534625892, 1.0488088482, 1.1441551071),
                           Mean = structure(9:12, names = letters[1:4]),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = 9:12,
                           Q1 = 9:12,
                           Median = 9:12,
                           Q3 = 9:12,
                           Max = 9:12,
                           MAD = 0)),
              summary(rcc2, 1L, "LaneID", elt = "exprsp1"))

  # Marginal summaries by Sample
  checkEquals(list(Endogenous =
                     cbind(GeomMean = c(1, 5, 9),
                           SizeFactor = c(0.2811442218, 1.4057211088, 2.5302979959),
                           MedPolSF = c(0.2, 1, 1.8),
                           Mean = structure(c(1, 5, 9), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = c(1, 5, 9),
                           Q1 = c(1, 5, 9),
                           Median = c(1, 5, 9),
                           Q3 = c(1, 5, 9),
                           Max = c(1, 5, 9),
                           MAD = 0),
                   Housekeeping =
                     cbind(GeomMean = c(4, 8, 12),
                           SizeFactor = c(0.5503212081, 1.1006424163, 1.6509636244),
                           MedPolSF = c(0.5, 1.0, 1.5),
                           Mean = structure(c(4, 8, 12), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = c(4, 8, 12),
                           Q1 = c(4, 8, 12),
                           Median = c(4, 8, 12),
                           Q3 = c(4, 8, 12),
                           Max = c(4, 8, 12),
                           MAD = 0),
                   Negative =
                     cbind(GeomMean = c(3, 7, 11),
                           SizeFactor = c(0.4889344008, 1.1408469352, 1.7927594696),
                           MedPolSF = c(0.4285714286, 1, 1.5714285714),
                           Mean = structure(c(3, 7, 11), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = c(3, 7, 11),
                           Q1 = c(3, 7, 11),
                           Median = c(3, 7, 11),
                           Q3 = c(3, 7, 11),
                           Max = c(3, 7, 11),
                           MAD = 0),
                   Positive =
                     cbind(GeomMean = c(2, 6, 10),
                           SizeFactor = c(0.405480133, 1.216440399, 2.027400665),
                           MedPolSF = c(1/3, 1, 5/3),
                           Mean = structure(c(2, 6, 10), names = sampleNames(rcc)),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = c(2, 6, 10),
                           Q1 = c(2, 6, 10),
                           Median = c(2, 6, 10),
                           Q3 = c(2, 6, 10),
                           Max = c(2, 6, 10),
                           MAD = 0)),
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
