# Coercion
setAs("NanoStringRccSet", "list",
      function(from) c(as.list(assayData(from)), fData(from), sData(from)))
setMethod("as.list", "NanoStringRccSet", function(x, ...) as(x, "list"))


# Accessing
setGeneric("sData", signature = "object",
           function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringRccSet",
          function(object) cbind(pData(object), pData(protocolData(object))))

setGeneric("svarLabels", signature = "object",
           function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringRccSet",
          function(object) c(varLabels(object), varLabels(protocolData(object))))

assayDataElement2 <- function(object, elt)
{
  if (elt %in% assayDataElementNames(object))
    assayDataElement(object, elt)
  else
    stop("'elt' not present in assayData(object)")
}


# Summarizing
.marginal.summary <- function(x)
{
  # Handle missing data
  if (anyNA(x))
    x <- x[!is.na(x)]

  # Calculate statistics
  quartiles <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(quartiles) <- c("Min", "Q1", "Median", "Q3", "Max")
  c("N"        = length(x),
    "Mean"     = mean(x),
    "SD"       = sd(x),
    "Skewness" = skewness(x),
    "Kurtosis" = kurtosis(x),
    geomMean(x),
    quartiles,
    "MAD"      = mad(x))
}

setMethod("summary", "NanoStringRccSet",
function(object, MARGIN, GROUP = NULL, elt = "exprs", ...)
{
  stopifnot(MARGIN %in% c(1L, 2L))
  FUN <- function(x) {
    mp <- medpolish(assayDataElement2(x, elt), eps = 1e-8, maxiter = 100L,
                    trace.iter = FALSE, na.rm = TRUE)
    cbind(t(esApply(x, MARGIN = MARGIN, FUN = .marginal.summary, elt = elt)),
          MedPolEff = mp[[ifelse(MARGIN == 1L, "row", "col")]])
  }
  if (is.null(GROUP)) {
    FUN(object)
  } else {
    esBy(object, GROUP = GROUP, FUN = FUN, simplify = FALSE)
  }
})


# Subsetting
setMethod("subset", "NanoStringRccSet",
function(x, subset, select, ...)
{
  kvs <- c(sData(x), fData(x))
  eval(substitute(with(kvs, x[subset, select])))
})

setGeneric("endogenousSubset", signature = "x",
           function(x, subset, select) standardGeneric("endogenousSubset"))
setMethod("endogenousSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(featureData(x)[["BarcodeClass"]] == "Endogenous"), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })

setGeneric("housekeepingSubset", signature = "x",
           function(x, subset, select) standardGeneric("housekeepingSubset"))
setMethod("housekeepingSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(featureData(x)[["BarcodeClass"]] == "Housekeeping"), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })

setGeneric("negativeControlSubset", signature = "x",
           function(x, subset, select) standardGeneric("negativeControlSubset"))
setMethod("negativeControlSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(featureData(x)[["BarcodeClass"]] == "Negative"), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })

setGeneric("positiveControlSubset", signature = "x",
           function(x, subset, select) standardGeneric("positiveControlSubset"))
setMethod("positiveControlSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(featureData(x)[["BarcodeClass"]] == "Positive"), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })

setGeneric("controlSubset", signature = "x",
           function(x, subset, select) standardGeneric("controlSubset"))
setMethod("controlSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(featureData(x)[["IsControl"]]), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })

setGeneric("nonControlSubset", signature = "x",
           function(x, subset, select) standardGeneric("nonControlSubset"))
setMethod("nonControlSubset", "NanoStringRccSet",
          function(x, subset, select) {
            call <- match.call()
            call$x <- x[which(!featureData(x)[["IsControl"]]), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })


# Looping
setMethod("esApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ..., elt = "exprs")
{
  stopifnot(MARGIN %in% c(1L, 2L))
  if (MARGIN == 1L)
    kvs <- sData(X)
  else
    kvs <- fData(X)

  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  multiassign(names(kvs), kvs, envir = e1)
  environment(FUN) <- e1

  apply(assayDataElement2(X, elt), MARGIN, FUN, ...)
})

setGeneric("endogenousApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("endogenousApply"))
setMethod("endogenousApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(endogenousSubset(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("housekeepingApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("housekeepingApply"))
setMethod("housekeepingApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(housekeepingSubset(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("negativeControlApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("negativeControlApply"))
setMethod("negativeControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(negativeControlSubset(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("positiveControlApply", signature = "X",
           function(X, MARGIN, FUN, ..., elt = "exprs")
             standardGeneric("positiveControlApply"))
setMethod("positiveControlApply", "NanoStringRccSet",
          function(X, MARGIN, FUN, ..., elt = "exprs")
            esApply(positiveControlSubset(X), MARGIN = MARGIN, FUN = FUN, ...,
                    elt = elt))

setGeneric("esBy", signature = "X",
           function(X, GROUP, FUN, ...) standardGeneric("esBy"))
setMethod("esBy", "NanoStringRccSet",
function(X, GROUP, FUN, ..., simplify = TRUE)
{
  featureNames <- fvarLabels(X)
  phenoNames <- varLabels(X)
  protocolNames <- varLabels(protocolData(X))
  choices <- c(structure(rep.int("featureData", length(featureNames)),
                         names = featureNames),
               structure(rep.int("phenoData", length(phenoNames)),
                         names = phenoNames),
               structure(rep.int("protocolData", length(protocolNames)),
                         names = protocolNames))
  GROUP <- choices[match.arg(GROUP, names(choices))]
  values <- do.call(GROUP, list(X))[[names(GROUP)]]
  keys <- sort(unique(values))
  names(keys) <- as.character(keys)
  if (GROUP == "featureData") {
    sapply(keys, function(k) FUN(X[values == k, ], ...),
           simplify = simplify)
  } else {
    sapply(keys, function(k) FUN(X[, values == k], ...),
           simplify = simplify)
  }
})


# Transforming
setMethod("transform", "NanoStringRccSet",
function(`_data`, ...)
{
  exprs <- as.list(substitute(list(...))[-1L])
  if (any(names(exprs) == "")) {
    stop("all arguments in '...' must be named")
  }
  aData <- assayData(`_data`)
  isLocked <- environmentIsLocked(aData)
  if (isLocked) {
    aData <- copyEnv(aData)
  }
  for (elt in names(exprs)) {
    assign(elt, eval(exprs[[elt]], as.list(aData), parent.frame()), aData)
  }
  if (isLocked) {
    lockEnvironment(aData)
    assayData(`_data`) <- aData
  }
  preproc(`_data`)[names(exprs)] <- exprs
  `_data`
})


# Evaluating
setMethod("with", "NanoStringRccSet",
          function(data, expr, ...)
            eval(substitute(expr), as(data, "list"), parent.frame()))
