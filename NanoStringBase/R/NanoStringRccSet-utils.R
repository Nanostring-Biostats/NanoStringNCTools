# Coercion
setAs("NanoStringRccSet", "list",
      function(from) c(as.list(assayData(from)), fData(from), sData(from),
                       list(signatureWeights = signatureWeights(from),
                            design = design(from))))
setMethod("as.list", "NanoStringRccSet", function(x, ...) as(x, "list"))


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
