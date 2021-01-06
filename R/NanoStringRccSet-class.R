setClassUnion("formulaOrNULL", c("formula", "NULL"))
setClass("NanoStringRccSet", contains = "ExpressionSet", slots = c(dimLabels = "character", 
    signatures = "SignatureSet", design = "formulaOrNULL"), prototype = prototype(new("VersionedBiobase", 
    versions = c(classVersion("ExpressionSet"), NanoStringRccSet = "1.0.0")), signatures = SignatureSet(), 
    design = NULL))
setMethod("show", signature = "NanoStringRccSet", function(object) {
    callNextMethod(object)
    cat("signature: ")
    if (length(signatures(object)) == 0L) 
        cat("none\n")
    else cat("use 'signatures(object)'")
})
setGeneric("NanoStringRccSet", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) standardGeneric("NanoStringRccSet"), signature = "assayData")
setMethod("NanoStringRccSet", "missing", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) {
    assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
    callGeneric(assayData = assayData, phenoData = phenoData, featureData = featureData, 
        experimentData = experimentData, annotation = annotation, protocolData = protocolData, 
        dimLabels = dimLabels, signatures = signatures, design = design, ...)
})
setMethod("NanoStringRccSet", "environment", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) {
    new2("NanoStringRccSet", assayData = assayData, phenoData = phenoData, featureData = featureData, 
        experimentData = experimentData, annotation = annotation, protocolData = protocolData, 
        dimLabels = dimLabels, signatures = signatures, design = design, ...)
})
setMethod("NanoStringRccSet", "matrix", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) {
    assayData <- assayDataNew(exprs = assayData)
    callGeneric(assayData = assayData, phenoData = phenoData, featureData = featureData, 
        experimentData = experimentData, annotation = annotation, protocolData = protocolData, 
        dimLabels = dimLabels, signatures = signatures, design = design, ...)
})
setMethod("NanoStringRccSet", "ExpressionSet", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) {
    callGeneric(assayData = copyEnv(assayData(assayData)), phenoData = Biobase::phenoData(assayData), 
        featureData = Biobase::featureData(assayData), experimentData = Biobase::experimentData(assayData), 
        annotation = Biobase::annotation(assayData), protocolData = Biobase::protocolData(assayData), 
        dimLabels = dimLabels, signatures = signatures, design = design, ...)
})
setMethod("NanoStringRccSet", "NanoStringRccSet", function(assayData, phenoData = annotatedDataFrameFrom(assayData, 
    byrow = FALSE), featureData = annotatedDataFrameFrom(assayData, byrow = TRUE), experimentData = MIAME(), 
    annotation = character(), protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE), 
    dimLabels = c("GeneName", "SampleID"), signatures = SignatureSet(), design = NULL, 
    ...) {
    callGeneric(assayData = copyEnv(assayData(assayData)), phenoData = Biobase::phenoData(assayData), 
        featureData = Biobase::featureData(assayData), experimentData = Biobase::experimentData(assayData), 
        annotation = Biobase::annotation(assayData), protocolData = Biobase::protocolData(assayData), 
        dimLabels = dimLabels(assayData), signatures = signatures(assayData), design = design(assayData), 
        ...)
})
setAs("ExpressionSet", "NanoStringRccSet", function(from) NanoStringRccSet(from))
