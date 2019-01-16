assayDataElement2 <- function(object, elt)
{
  if (elt %in% assayDataElementNames(object))
    assayDataElement(object, elt)
  else
    stop("'elt' not present in assayData(object)")
}

# sData Accessor
setGeneric("sData", signature = "object",
           function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringRccSet",
          function(object) cbind(pData(object), pData(protocolData(object))))

# svarLabels Accessor
setGeneric("svarLabels", signature = "object",
           function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringRccSet",
          function(object) c(varLabels(object), varLabels(protocolData(object))))

# signatureWeights Accessor and Replacer
setGeneric("signatureWeights", signature = "object",
           function(object) standardGeneric("signatureWeights"))
setMethod("signatureWeights", "NanoStringRccSet",
          function(object) object@signatureWeights)

setGeneric("signatureWeights<-", signature = c("object", "value"),
           function(object, value) standardGeneric("signatureWeights<-"))
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "NumericList"),
                 function(object, value) {
                   object@signatureWeights <- value
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "ANY"),
                 function(object, value) {
                   object@signatureWeights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "list"),
                 function(object, value) {
                   value <- lapply(value, function(x) {
                     if(is.matrix(x) && ncol(x) == 1L)
                       structure(x[, 1L, drop = TRUE],
                                 names = rownames(x))
                     else
                       x
                   })
                   object@signatureWeights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "NULL"),
                 function(object, value) {
                   object@signatureWeights <- NumericList()
                   object
                 })

# signatureScores Accessor and Replacer
setGeneric("signatureScores", signature = "object",
           function(object, ...) standardGeneric("signatureScores"))
setMethod("signatureScores", "NanoStringRccSet",
          function(object, elt = "exprs") {
            exprs <- assayDataElement2(object, elt)
            rownames(exprs) <- featureData(object)[["GeneName"]]
            t(sapply(signatureWeights(object),
                     function(wts){
                       if ("(Intercept)" %in% names(wts))
                         exprs <- rbind("(Intercept)" = 1, exprs)
                       exprs <- exprs[names(wts), , drop = FALSE]
                       colSums(wts * exprs)
                     }))
          })

# design Accessor and Replacer
setMethod("design", "NanoStringRccSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringRccSet", "formula"),
                 function(object, value) {
                   object@design <- value
                   object
                 })
setReplaceMethod("design", c("NanoStringRccSet", "ANY"),
                 function(object, value) {
                   object@design <- as.formula(value)
                   object
                 })
setReplaceMethod("design", c("NanoStringRccSet", "NULL"),
                 function(object, value) {
                   object@design <- NULL
                   object
                 })
