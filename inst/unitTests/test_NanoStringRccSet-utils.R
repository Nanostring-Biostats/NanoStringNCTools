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

# Looping
test_NanoStringRccSet_utils_esApply <- function() {
  rcc2 <- transform(rcc, log1p_exprs = log1p(exprs))

  checkIdentical(apply(exprs(rcc2), 1L, mean), esApply(rcc2, 1L, mean))
  checkIdentical(apply(exprs(rcc2), 2L, mean), esApply(rcc2, 2L, mean))

  checkIdentical(apply(assayDataElement(rcc2, "log1p_exprs"), 1L, mean),
                 esApply(rcc2, 1L, mean, elt = "log1p_exprs"))
  checkIdentical(apply(assayDataElement(rcc2, "log1p_exprs"), 2L, mean),
                 esApply(rcc2, 2L, mean, elt = "log1p_exprs"))
}

# Transforming
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

# Evaluating
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
