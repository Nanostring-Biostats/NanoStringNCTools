# Show
setMethod("show", signature = "NanoStringRccSet",
          function(object) {
            callNextMethod(object)
            cat("signatureWeights: ")
            if (length(signatureWeights(object)) == 0L)
              cat("none\n")
            else
              cat("use 'signatureWeights(object)'")
          })

# Coercion
setAs("NanoStringRccSet", "list",
      function(from) c(as.list(assayData(from)), fData(from), sData(from),
                       list(signatureWeights = signatureWeights(from),
                            design = design(from))))
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

setGeneric("modelData", signature = "object",
           function(object, formula = design(object), ...)
             standardGeneric("modelData"))
setMethod("modelData", "NanoStringRccSet",
function(object, formula = design(object), ...)
{
  if (is.null(formula))
    stop("\"formula\" argument is missing")
  vars <- all.vars(formula)
  hasFeatureVars <- any(vars %in% fvarLabels(object))
  hasSampleVars  <- any(vars %in% svarLabels(object))
  if (hasFeatureVars && hasSampleVars)
    stop("\"formula\" argument cannot use both feature and sample variables")
  if (hasFeatureVars)
    data <- fData(object)
  else if (hasSampleVars)
    data <- sData(object)
  else
    data <- NULL
  if (nargs() > 2L) {
    if (is.null(data)) {
      data <- cbind.data.frame(...)
      if (!identical(rownames(data), featureNames(object)) &&
          !identical(rownames(data), sampleNames(object))) {
        stop("data in \"...\" do not match 'featureNames' or 'sampleNames'")
      }
    } else {
      data <- cbind(data, ...)
    }
  }
  model.frame(formula, data)
})

assayDataElement2 <- function(object, elt)
{
  if (elt %in% assayDataElementNames(object))
    assayDataElement(object, elt)
  else
    stop("'elt' not present in assayData(object)")
}

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
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "NULL"),
                 function(object, value) {
                   object@signatureWeights <- NumericList()
                   object
                 })

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


# Summarizing
.marginal.summary <- function(x, log2scale = TRUE)
{
  # Handle missing data
  if (anyNA(x))
    x <- x[!is.na(x)]

  # Calculate statistics
  quartiles <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(quartiles) <- c("Min", "Q1", "Median", "Q3", "Max")
  if (log2scale) {
    log2X <- log2t(x, thresh = 0.5)
    stats <- c("GeomMean"   = geomMean(x),
               "SizeFactor" = NA_real_,
               "MeanLog2"   = mean(log2X),
               "SDLog2"     = sd(log2X))
  } else {
    stats <- c("Mean"       = mean(x),
               "SD"         = sd(x),
               "Skewness"   = skewness(x),
               "Kurtosis"   = kurtosis(x))
  }
  c(stats, quartiles)
}

setMethod("summary", "NanoStringRccSet",
function(object, MARGIN = 2L, GROUP = NULL, log2scale = TRUE, elt = "exprs", ...)
{
  stopifnot(MARGIN %in% c(1L, 2L))
  FUN <- function(x) {
    stats <- t(esApply(x, MARGIN = MARGIN, FUN = .marginal.summary,
                       log2scale = log2scale, elt = elt))
    if (log2scale) {
      stats[,"SizeFactor"] <- 2^(stats[,"MeanLog2"] - mean(stats[,"MeanLog2"]))
    }
    stats
  }
  if (is.null(GROUP)) {
    FUN(object)
  } else {
    esBy(object, GROUP = GROUP, FUN = FUN, simplify = FALSE)
  }
})


# Subsetting
setMethod("[", "NanoStringRccSet",
function(x, i, j, ..., drop = FALSE)
{
  x <- callNextMethod()
  weights <- signatureWeights(x)
  if (length(weights) > 0L) {
    genes <- featureData(x)[["GeneName"]]
    keep <- unlist(lapply(weights, function(y) all(names(y) %in% genes)))
    if (!all(keep)) {
      signatureWeights(x) <- weights[keep]
    }
  }
  x
})

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

setGeneric("signatureSubset", signature = "x",
           function(x, subset, select) standardGeneric("signatureSubset"))
setMethod("signatureSubset", "NanoStringRccSet",
          function(x, subset, select) {
            genes <- unique(names(unlist(unname(signatureWeights(x)))))
            call <- match.call()
            call$x <- x[which(featureData(x)[["GeneName"]] %in% genes), ]
            call[[1L]] <- as.name("subset")
            eval(call, parent.frame())
          })


# Looping
setMethod("esApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ..., elt = "exprs")
{
  stopifnot(MARGIN %in% c(1L, 2L))
  if (MARGIN == 1L)
    kvs <- c(sData(X), list(design = design(X)))
  else
    kvs <- c(fData(X), list(signatureWeights = signatureWeights(X)))

  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  multiassign(names(kvs), kvs, envir = e1)
  environment(FUN) <- e1

  apply(assayDataElement2(X, elt), MARGIN, FUN, ...)
})

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
