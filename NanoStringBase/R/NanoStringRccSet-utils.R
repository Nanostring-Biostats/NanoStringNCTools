# Subsetting
setGeneric("controlSet", signature = "object",
           function(object) standardGeneric("controlSet"))
setMethod("controlSet", "NanoStringRccSet",
          function(object)
            object[featureData(object)[["IsControl"]], ])

setGeneric("nonControlSet", signature = "object",
           function(object) standardGeneric("nonControlSet"))
setMethod("nonControlSet", "NanoStringRccSet",
          function(object)
            object[!featureData(object)[["IsControl"]], ])

setGeneric("endogenousSet", signature = "object",
           function(object) standardGeneric("endogenousSet"))
setMethod("endogenousSet", "NanoStringRccSet",
          function(object)
            object[featureData(object)[["BarcodeClass"]] == "Endogenous", ])

setGeneric("housekeepingSet", signature = "object",
           function(object) standardGeneric("housekeepingSet"))
setMethod("housekeepingSet", "NanoStringRccSet",
          function(object)
            object[featureData(object)[["BarcodeClass"]] == "Housekeeping", ])

setGeneric("negativeControlSet", signature = "object",
           function(object) standardGeneric("negativeControlSet"))
setMethod("negativeControlSet", "NanoStringRccSet",
          function(object)
            object[featureData(object)[["BarcodeClass"]] == "Negative", ])

setGeneric("positiveControlSet", signature = "object",
           function(object) standardGeneric("positiveControlSet"))
setMethod("positiveControlSet", "NanoStringRccSet",
          function(object)
            object[featureData(object)[["BarcodeClass"]] == "Positive", ])


# Utilities
setMethod("esApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ...)
{
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  kvs <- c(pData(X), pData(protocolData(X)), fData(X))
  multiassign(names(kvs), kvs, envir = e1)
  environment(FUN) <- e1
  apply(exprs(X), MARGIN, FUN, ...)
})

setMethod("subset", "NanoStringRccSet",
function(x, subset, select, ...)
{
  kvs <- c(pData(x), pData(protocolData(x)), fData(x))
  eval(substitute(with(kvs, x[subset, select])))
})
