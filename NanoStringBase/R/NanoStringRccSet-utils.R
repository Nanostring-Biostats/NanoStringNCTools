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
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("endogenousApply"))
setMethod("endogenousApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(endogenousSet(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("housekeepingApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("housekeepingApply"))
setMethod("housekeepingApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(housekeepingSet(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("negativeControlApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("negativeControlApply"))
setMethod("negativeControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(negativeControlSet(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("positiveControlApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("positiveControlApply"))
setMethod("positiveControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(positiveControlSet(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))


# Utilities
setMethod("esApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ..., elt = "exprs")
{
  stopifnot(MARGIN %in% c(1L, 2L))
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  if (MARGIN == 1L)
    kvs <- cbind(pData(X), pData(protocolData(X)))
  else
    kvs <- fData(X)
  multiassign(names(kvs), kvs, envir = e1)
  environment(FUN) <- e1
  apply(assayDataElement(X, elt), MARGIN, FUN, ...)
})

setGeneric("esSweep", signature = "x",
           function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...,
                    fromElt = "exprs", toElt, validate = TRUE)
             standardGeneric("esSweep"))
setMethod("esSweep", "NanoStringRccSet",
function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...,
         fromElt = "exprs", toElt, validate = TRUE)
{
  stopifnot(MARGIN %in% c(1L, 2L))
  if (missing(toElt))
    stop("argument \"toElt\" is missing, with no default")

  # Construct FUN argument
  FUN <- match.fun(FUN)
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  if (MARGIN == 1L)
    kvs <- fData(x)
  else
    kvs <- cbind(pData(x), pData(protocolData(x)))
  multiassign(names(kvs), kvs, envir = e1)
  environment(FUN) <- e1

  # Evaluate STATS argument
  stats <- try(eval(substitute(STATS), e1), silent = TRUE)
  if (inherits(stats, "try-error"))
    STATS <- eval(substitute(STATS), parent.frame())
  else
    STATS <- stats

  # Calculate matrix
  value <- sweep(assayDataElement(x, fromElt), MARGIN = MARGIN, STATS = STATS,
                 FUN = FUN, check.margin = check.margin, ...)

  # Modify return value
  assayDataElement(x, toElt, validate = validate) <- value
  preproc(x)[[toElt]] <- match.call()

  x
})

setMethod("subset", "NanoStringRccSet",
function(x, subset, select, ...)
{
  kvs <- c(pData(x), pData(protocolData(x)), fData(x))
  eval(substitute(with(kvs, x[subset, select])))
})
