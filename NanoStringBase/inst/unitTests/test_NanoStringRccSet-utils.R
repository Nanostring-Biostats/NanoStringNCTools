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
                                       Owner = rep("", 3L),
                                       Comments = rep("DNA-RNA-Protein", 3L),
                                       Date = as.Date(rep("1999-12-31", 3L)),
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

test_NanoStringRccSet_utils_sData <- function() {
  checkIdentical(cbind(pData(rcc), pData(protocolData(rcc))), sData(rcc))
}

test_NanoStringRccSet_utils_svarLabels <- function() {
  checkIdentical(c(varLabels(rcc), varLabels(protocolData(rcc))), svarLabels(rcc))
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
