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
                            dimLabels = c("sampleNames", "sampleColumns"))
       )

# Subsetting
test_NanoStringRccSet_subset <- function() {
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Endogenous", ],
              subset(rcc, CodeClass == "Endogenous"))
  #checkEquals(IRanges::NumericList(x = c(a = 1), compress = FALSE),
  #            weights(signatures(subset(rcc, CodeClass == "Endogenous"))))

  checkEquals(rcc[, phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, select = Treatment == "A"))
  checkEquals(weights(signatures(rcc)),
              weights(signatures(subset(rcc, select = Treatment == "A"))))

  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Endogenous",
                  phenoData(rcc)[["Treatment"]] == "A"],
              subset(rcc, CodeClass == "Endogenous", Treatment == "A"))
  #checkEquals(IRanges::NumericList(x = c(a = 1), compress = FALSE),
  #            weights(signatures(subset(rcc, CodeClass == "Endogenous",
  #                                      Treatment == "A"))))
}

test_NanoStringRccSet_subset_in_function <- function() {
  subsetFUN1 <- function(object, x) subset(object, CodeClass == x)
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Endogenous", ],
              subsetFUN1(rcc, "Endogenous"))

  subsetFUN2 <- function(object, x) subset(object, select = Treatment == x)
  checkEquals(rcc[, phenoData(rcc)[["Treatment"]] == "A"],
              subsetFUN2(rcc, "A"))

  subsetFUN3 <-
    function(object, x, y) subset(object, CodeClass == x, Treatment == y)
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Endogenous",
                  phenoData(rcc)[["Treatment"]] == "A"],
              subsetFUN3(rcc, "Endogenous", "A"))
}

test_NanoStringRccSet_endogenousSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Endogenous", ],
              endogenousSubset(rcc))
}

test_NanoStringRccSet_housekeepingSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Housekeeping", ],
              housekeepingSubset(rcc))
}

test_NanoStringRccSet_negativeControlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Negative", ],
              negativeControlSubset(rcc))
}

test_NanoStringRccSet_positiveControlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Positive", ],
              positiveControlSubset(rcc))

  checkEquals(rcc[featureData(rcc)[["CodeClass"]] == "Positive" &
                  featureData(rcc)[["ControlConc"]] >= 0.5, ],
              positiveControlSubset(rcc, subset = ControlConc >= 0.5))
}

test_NanoStringRccSet_controlSubset <- function() {
  checkEquals(rcc[featureData(rcc)[["IsControl"]], ], controlSubset(rcc))
}

test_NanoStringRccSet_nonControlSubset <- function() {
  checkEquals(rcc[!featureData(rcc)[["IsControl"]], ], nonControlSubset(rcc))
}
