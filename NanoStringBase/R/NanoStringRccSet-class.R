setClassUnion("formulaOrNULL", c("formula", "NULL"))

# Class definition
setClass("NanoStringRccSet",
         contains = "ExpressionSet",
         slots = c(signatureWeights = "NumericList",
                   design = "formulaOrNULL"),
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(classVersion("ExpressionSet"),
                            NanoStringRccSet = "1.0.0")),
           signatureWeights = NumericList(),
           design = NULL))

# Initialization method
setMethod("initialize", "NanoStringRccSet",
function(.Object, signatureWeights = NumericList(), ...)
{
  callNextMethod(.Object,
                 signatureWeights = as(signatureWeights, "NumericList"),
                 ...)
})

# Show method
setMethod("show", signature = "NanoStringRccSet",
function(object) {
  callNextMethod(object)
  cat("signatureWeights: ")
  if (length(signatureWeights(object)) == 0L)
    cat("none\n")
  else
    cat("use 'signatureWeights(object)'")
})

# Constructors
setGeneric("NanoStringRccSet",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         design = NULL,
         ...)
  standardGeneric("NanoStringRccSet"),
signature = "assayData")

setMethod("NanoStringRccSet", "missing",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         design = NULL,
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              signatureWeights = signatureWeights, design = design, ...)
})

setMethod("NanoStringRccSet", "environment",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         design = NULL,
         ...)
{
  new2("NanoStringRccSet",
       assayData = assayData,
       phenoData = phenoData,
       featureData = featureData,
       experimentData = experimentData,
       annotation = annotation,
       protocolData = protocolData,
       signatureWeights = signatureWeights,
       design = design,
       ...)
})

setMethod("NanoStringRccSet", "matrix",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         design = NULL,
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              signatureWeights = signatureWeights, design = design, ...)
})

setMethod("NanoStringRccSet", "NanoStringRccSet",
function(assayData,
         phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
         signatureWeights = NumericList(),
         design = NULL,
         ...)
{
  callGeneric(assayData = copyEnv(assayData(assayData)),
              phenoData = Biobase::phenoData(assayData),
              featureData = Biobase::featureData(assayData),
              experimentData = Biobase::experimentData(assayData),
              annotation = Biobase::annotation(assayData),
              protocolData = Biobase::protocolData(assayData),
              signatureWeights = signatureWeights,
              design = design,
              ...)
})
