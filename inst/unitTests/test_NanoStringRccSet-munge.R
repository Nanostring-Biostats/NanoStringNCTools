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
                            dimLabels = c("sampleNames", "sampleColumns")),
       signatureWeights =
         list(x = c(a = 1), y = c(b = 1/3, d = 2/3), z = c(a = 2, c = 4)))

# Munging data for plotting
test_NanoStringRccSet_munge_exception_missing_mapping <- function() {
  checkException(munge(rcc))
}

test_NanoStringRccSet_munge_exception_feature_and_sample_vars <- function() {
  checkException(munge(rcc, ~ CodeClass + BindingDensity))
}

test_NanoStringRccSet_munge_exception_mismatch_summaries <- function() {
  checkException(munge(rcc, ~ MeanLog2 + Mean + BindingDensity))
}

test_NanoStringRccSet_munge_exception_mismatch_extradata <- function() {
  checkException(munge(rcc, ~ GNP, longley))
}

test_NanoStringRccSet_munge_exception_ambiguous_summaries <- function() {
  checkException(munge(rcc, ~ MeanLog2 + SDLog2))
}

test_NanoStringRccSet_munge_featureData <- function() {
  target <- data.frame(FeatureName = featureNames(rcc),
                       CodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
                       stringsAsFactors = FALSE)
  checkIdentical(target, munge(rcc, ~ CodeClass))

  target <-
    data.frame(FeatureName = featureNames(rcc),
               MeanLog2 = c(1.333333333, 1.830617699, 2.302296865, 2.617249680),
               CodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ MeanLog2 + CodeClass))

  target <-
    data.frame(FeatureName = featureNames(rcc),
               Mean = 4:7,
               CodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ Mean + CodeClass))

  target <-
    data.frame(FeatureName = featureNames(rcc),
               Median = 4:7,
               CodeClass = c("Endogenous", "Positive", "Negative", "Housekeeping"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ Median + CodeClass))
}

test_NanoStringRccSet_munge_sampleData <- function() {
  target <- data.frame(SampleName = sampleNames(rcc),
                       Treatment = c("A", "A", "B"),
                       Age = c(58L, 42L, 27L),
                       stringsAsFactors = FALSE)
  checkIdentical(target, munge(rcc, ~ Treatment + Age))

  target <- data.frame(SampleName = sampleNames(rcc),
                       Treatment = c("A", "A", "B"),
                       LaneID = 1:3,
                       stringsAsFactors = FALSE)
  checkIdentical(target, munge(rcc, ~ Treatment + LaneID))

  target <-
    data.frame(SampleName = sampleNames(rcc),
               MeanLog2 = c(0.3962406252, 2.4285613794, 3.2378211787),
               Treatment = c("A", "A", "B"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ MeanLog2 + Treatment))

  target <-
    data.frame(SampleName = sampleNames(rcc),
               Mean = c(1.5, 5.5, 9.5),
               Treatment = c("A", "A", "B"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ Mean + Treatment))

  target <-
    data.frame(SampleName = sampleNames(rcc),
               Median = c(1.5, 5.5, 9.5),
               Treatment = c("A", "A", "B"),
               stringsAsFactors = FALSE)
  checkEquals(target, munge(rcc, ~ Median + Treatment))

  newdata <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  target <- data.frame(SampleName = sampleNames(rcc),
                       V1 = 11:13,
                       stringsAsFactors = FALSE)
  checkIdentical(target, munge(rcc, ~ V1, newdata))

  newdata <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  target <- data.frame(SampleName = sampleNames(rcc),
                       Treatment = c("A", "A", "B"),
                       V1 = 11:13,
                       stringsAsFactors = FALSE)
  checkIdentical(target, munge(rcc, ~ Treatment + V1, newdata))
}
