assayDataElement2 <- function(object, elt) {
    if (elt %in% assayDataElementNames(object)) 
        assayDataElement(object, elt)
    else stop("'elt' not present in assayData(object)")
}
setGeneric("sData", signature = "object", function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringRccSet", function(object) cbind(pData(object), pData(protocolData(object))))
setGeneric("svarLabels", signature = "object", function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringRccSet", function(object) c(varLabels(object), varLabels(protocolData(object))))
setMethod("dimLabels", "NanoStringRccSet", function(object) object@dimLabels)
setReplaceMethod("dimLabels", c("NanoStringRccSet", "character"), function(object, value) {
    object@dimLabels <- value
    object
})
setMethod("design", "NanoStringRccSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringRccSet", "formula"), function(object, value) {
    object@design <- value
    object
})
setReplaceMethod("design", c("NanoStringRccSet", "ANY"), function(object, value) {
    object@design <- as.formula(value)
    object
})
setReplaceMethod("design", c("NanoStringRccSet", "NULL"), function(object, value) {
    object@design <- NULL
    object
})
