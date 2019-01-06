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

test_NanoStringRccSet_utils_modelData <- function() {
  # Exceptions
  checkException(modelData(rcc))

  rcc2 <- rcc
  design(rcc2) <- ~ BarcodeClass + BindingDensity
  checkException(modelData(rcc2))

  checkException(modelData(rcc, ~ GNP, longley))

  # Valid results
  design(rcc2) <- ~ Treatment + Age
  target <- pData(rcc2)
  attr(target, "terms") <- terms(design(rcc2))
  checkEquals(target, modelData(rcc2))

  target <- fData(rcc)[,"BarcodeClass", drop = FALSE]
  attr(target, "terms") <- terms(~ BarcodeClass)
  checkEquals(target, modelData(rcc, ~ BarcodeClass))

  target <- sData(rcc)[, c("Treatment", "LaneID")]
  attr(target, "terms") <- terms(~ Treatment + LaneID)
  checkEquals(target, modelData(rcc, ~ Treatment + LaneID))

  target <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  attr(target, "terms") <- terms(~ V1)
  checkEquals(target, modelData(rcc, ~ V1, target))

  newdata <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  target <- cbind(pData(rcc)[, "Treatment", drop = FALSE], newdata)
  attr(target, "terms") <- terms(~ Treatment + V1)
  checkEquals(target, modelData(rcc, ~ Treatment + V1, newdata))
}

# Sumarizing
test_NanoStringRccSet_utils_summary <- function() {
  rcc2 <- transform(rcc, log2_exprs = log2t(exprs))

  # Marginal summaries by Feature
  checkEquals(cbind(GeomMean = c(2.519842100, 3.556893304, 4.932424149, 6.135792440),
                    SizeFactor = c(0.6209112338, 0.8764497626, 1.2153926486, 1.5119131688),
                    MeanLog2 = c(1.333333333, 1.830617699, 2.302296865, 2.617249680),
                    SDLog2 = c(2.0816659995, 1.6410806067, 1.1864916416, 0.9515847943),
                    Min = structure(0:3, names = letters[1:4]),
                    Q1 = 2:5,
                    Median = 4:7,
                    Q3 = 6:9,
                    Max = 8:11),
              summary(rcc2, 1L))
  checkEquals(cbind(Mean = c(1.333333333, 1.830617699, 2.302296865, 2.617249680),
                    SD = c(2.0816659995, 1.6410806067, 1.1864916416, 0.9515847943),
                    Skewness = c(-1.2933427807, -1.2264691300, -1.0112175576, -0.8631188095),
                    Kurtosis = NA_real_,
                    Min = structure(c(-1, 0, 1, 1.584962501), names = letters[1:4]),
                    Q1 = c(0.5, 1.160964047, 1.792481250, 2.196158711),
                    Median = c(2, 2.321928095, 2.584962501, 2.807354922),
                    Q3 = c(2.5, 2.745926548, 2.953445298, 3.133393270),
                    Max = c(3, 3.169925001, 3.321928095, 3.459431619)),
              summary(rcc2, 1L, elt = "log2_exprs", log2scale = FALSE))

  # Marginal summaries by Sample
  checkEquals(cbind(GeomMean = c(1.316074013, 5.383563271, 9.433683366),
                    SizeFactor = c(0.3242922004, 1.3265572922, 2.3245424696),
                    MeanLog2 = c(0.3962406252, 2.4285613794, 3.2378211787),
                    SDLog2 = c(1.1378458590, 0.3478416172, 0.1977826790),
                    Min = structure(c(0, 4, 8), names = sampleNames(rcc)),
                    Q1 = c(0.75, 4.75, 8.75),
                    Median = c(1.5, 5.5, 9.5),
                    Q3 = c(2.25, 6.25, 10.25),
                    Max = c(3, 7, 11)),
              summary(rcc2, 2L))
  checkEquals(cbind(Mean = c(0.3962406252, 2.4285613794, 3.2378211787),
                    SD = c(1.1378458590, 0.3478416172, 0.1977826790),
                    Skewness = c(-0.4002032906, -0.3444852394, -0.1969256558),
                    Kurtosis = c(-1.6584093446, -0.9658168584, -1.1224275064),
                    Min = structure(c(-1, 2, 3), names = sampleNames(rcc)),
                    Q1 = c(-0.25, 2.241446071, 3.127443751),
                    Median = c(0.5, 2.453445298, 3.245926548),
                    Q3 = c(1.146240625, 2.640560606, 3.356303976),
                    Max = c(1.584962501, 2.807354922, 3.459431619)),
              summary(rcc2, 2L, elt = "log2_exprs", log2scale = FALSE))
}

test_NanoStringRccSet_utils_summary_GROUP <- function() {
  rcc2 <- transform(rcc, log2_exprs = log2t(exprs))

  # Marginal summaries by Feature
  checkEquals(list(A =
                     cbind(GeomMean = c(1.414213562, 2.236067977, 3.464101615, 4.582575695),
                           SizeFactor = c(0.5313001399, 0.8400592816, 1.3014142429, 1.7216092197),
                           MeanLog2 = c(0.5, 1.160964047, 1.792481250, 2.196158711),
                           SDLog2 = c(2.1213203436, 1.6418511013, 1.1207377322, 0.8643619704),
                           Min = structure(0:3, names = letters[1:4]),
                           Q1 = 1:4,
                           Median = 2:5,
                           Q3 = 3:6,
                           Max = 4:7),
                   B =
                     cbind(GeomMean = 8:11,
                           SizeFactor = c(0.8480250704, 0.9540282042, 1.0600313380, 1.1660344717),
                           MeanLog2 = c(3, 3.169925001, 3.321928095, 3.459431619),
                           SDLog2 = NA_real_,
                           Min = structure(8:11, names = letters[1:4]),
                           Q1 = 8:11,
                           Median = 8:11,
                           Q3 = 8:11,
                           Max = 8:11)),
              summary(rcc2, 1L, "Treatment"))
  checkEquals(list("1" =
                     cbind(Mean = c(-1, 0, 1, 1.584962501),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = structure(c(-1, 0, 1, 1.584962501), names = letters[1:4]),
                           Q1 = c(-1, 0, 1, 1.584962501),
                           Median = c(-1, 0, 1, 1.584962501),
                           Q3 = c(-1, 0, 1, 1.584962501),
                           Max = c(-1, 0, 1, 1.584962501)),
                   "2" =
                     cbind(Mean = c(2, 2.321928095, 2.584962501, 2.807354922),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = structure(c(2, 2.321928095, 2.584962501, 2.807354922), names = letters[1:4]),
                           Q1 = c(2, 2.321928095, 2.584962501, 2.807354922),
                           Median = c(2, 2.321928095, 2.584962501, 2.807354922),
                           Q3 = c(2, 2.321928095, 2.584962501, 2.807354922),
                           Max = c(2, 2.321928095, 2.584962501, 2.807354922)),
                   "3" =
                     cbind(Mean = c(3, 3.169925001, 3.321928095, 3.459431619),
                           SD = NA_real_,
                           Skewness = NA_real_,
                           Kurtosis = NA_real_,
                           Min = structure(c(3, 3.169925001, 3.321928095, 3.459431619), names = letters[1:4]),
                           Q1 = c(3, 3.169925001, 3.321928095, 3.459431619),
                           Median = c(3, 3.169925001, 3.321928095, 3.459431619),
                           Q3 = c(3, 3.169925001, 3.321928095, 3.459431619),
                           Max = c(3, 3.169925001, 3.321928095, 3.459431619))),
              summary(rcc2, 1L, "LaneID", elt = "log2_exprs", log2scale = FALSE))

  # Marginal summaries by Sample
  checkEquals(list(Endogenous =
                     cbind(GeomMean = c(0.5, 4, 8),
                           SizeFactor = c(0.1984251315, 1.5874010520, 3.1748021039),
                           MeanLog2 = c(-1, 2, 3),
                           SDLog2 = NA_real_,
                           Min = structure(c(0, 4, 8), names = sampleNames(rcc)),
                           Q1 = c(0, 4, 8),
                           Median = c(0, 4, 8),
                           Q3 = c(0, 4, 8),
                           Max = c(0, 4, 8)),
                   Housekeeping =
                     cbind(GeomMean = c(3, 7, 11),
                           SizeFactor = c(0.4889344008, 1.1408469352, 1.7927594696),
                           MeanLog2 = c(1.584962501, 2.807354922, 3.459431619),
                           SDLog2 = NA_real_,
                           Min = structure(c(3, 7, 11), names = sampleNames(rcc)),
                           Q1 = c(3, 7, 11),
                           Median = c(3, 7, 11),
                           Q3 = c(3, 7, 11),
                           Max = c(3, 7, 11)),
                   Negative =
                     cbind(GeomMean = c(2, 6, 10),
                           SizeFactor = c(0.405480133, 1.216440399, 2.027400665),
                           MeanLog2 = c(1, 2.584962501, 3.321928095),
                           SDLog2 = NA_real_,
                           Min = structure(c(2, 6, 10), names = sampleNames(rcc)),
                           Q1 = c(2, 6, 10),
                           Median = c(2, 6, 10),
                           Q3 = c(2, 6, 10),
                           Max = c(2, 6, 10)),
                   Positive =
                     cbind(GeomMean = c(1, 5, 9),
                           SizeFactor = c(0.2811442218, 1.4057211088, 2.5302979959),
                           MeanLog2 = c(0, 2.321928095, 3.169925001),
                           SDLog2 = NA_real_,
                           Min = structure(c(1, 5, 9), names = sampleNames(rcc)),
                           Q1 = c(1, 5, 9),
                           Median = c(1, 5, 9),
                           Q3 = c(1, 5, 9),
                           Max = c(1, 5, 9))),
              summary(rcc2, 2L, "BarcodeClass"))
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
                    exprs_scaled = sweep(exprs, 2L, c(2, 1, 0.5), FUN = "*"),
                    exprs_thresh = pmax(exprs_scaled - 2L, 0L))
  checkTrue(validObject(rcc2))
  checkEquals(sweep(exprs(rcc), 2L, c(2, 1, 0.5), FUN = "*"),
              assayDataElement(rcc2, "exprs_scaled"))
  checkEquals(pmax(sweep(exprs(rcc), 2L, c(2, 1, 0.5), FUN = "*") - 2L, 0L),
                 assayDataElement(rcc2, "exprs_thresh"))
  checkIdentical(list(exprs_scaled = substitute(sweep(exprs, 2L, c(2, 1, 0.5), FUN = "*")),
                      exprs_thresh = substitute(pmax(exprs_scaled - 2L, 0L))),
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
