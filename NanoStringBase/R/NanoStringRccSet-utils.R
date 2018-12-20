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


# Apply
setGeneric("endogenousApply", signature = "X",
           function(X, MARGIN, FUN, ...) standardGeneric("endogenousApply"))
setMethod("endogenousApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ...)
            esApply(endogenousSet(X), MARGIN = MARGIN, FUN = FUN, ...))

setGeneric("housekeepingApply", signature = "X",
           function(X, MARGIN, FUN, ...) standardGeneric("housekeepingApply"))
setMethod("housekeepingApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ...)
            esApply(housekeepingSet(X), MARGIN = MARGIN, FUN = FUN, ...))

setGeneric("negativeControlApply", signature = "X",
           function(X, MARGIN, FUN, ...) standardGeneric("negativeControlApply"))
setMethod("negativeControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ...)
            esApply(negativeControlSet(X), MARGIN = MARGIN, FUN = FUN, ...))

setGeneric("positiveControlApply", signature = "X",
           function(X, MARGIN, FUN, ...) standardGeneric("positiveControlApply"))
setMethod("positiveControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ...)
            esApply(positiveControlSet(X), MARGIN = MARGIN, FUN = FUN, ...))


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
