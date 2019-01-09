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

# Molding data for plotting
test_NanoStringRccSet_mold <- function() {
  # Exceptions
  checkException(mold(rcc))

  rcc2 <- rcc
  design(rcc2) <- ~ BarcodeClass + BindingDensity
  checkException(mold(rcc2))

  checkException(mold(rcc, ~ GNP, longley))

  # Valid results
  design(rcc2) <- ~ Treatment + Age
  target <- pData(rcc2)
  checkIdentical(target, mold(rcc2))

  target <- fData(rcc)[,"BarcodeClass", drop = FALSE]
  checkIdentical(target, mold(rcc, ~ BarcodeClass))

  target <- sData(rcc)[, c("Treatment", "LaneID")]
  checkIdentical(target, mold(rcc, ~ Treatment + LaneID))

  target <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  checkIdentical(target, mold(rcc, ~ V1, target))

  newdata <- data.frame(V1 = 11:13, row.names = sampleNames(rcc))
  target <- cbind(pData(rcc)[, "Treatment", drop = FALSE], newdata)
  checkIdentical(target, mold(rcc, ~ Treatment + V1, newdata))
}
