# Coercion
setAs("NanoStringRccSet", "list",
      function(from) c(as.list(assayData(from)), fData(from), sData(from),
                       list(signatures = signatures(from),
                            design = design(from))))
setMethod("as.list", "NanoStringRccSet", function(x, ...) as(x, "list"))


# Looping
setGeneric("assayDataApply", signature = "X",
           function(X, MARGIN, FUN, ...)
             standardGeneric("assayDataApply"))
setMethod("assayDataApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ..., elt = "exprs")
{
  stopifnot(MARGIN %in% c(1L, 2L))
  if (MARGIN == 1L) {
    df <- fData(X)
    kvs <- c(sData(X), list(design = design(X)))
  } else {
    df <- sData(X)
    kvs <- fData(X)
  }
  mat <- assayDataElement2(X, elt)
  .apply(X = mat, MARGIN = MARGIN, FUN = FUN, ..., .df = df, .kvs = kvs)
})

setGeneric("signatureScoresApply", signature = "X",
           function(X, MARGIN, FUN, ...)
             standardGeneric("signatureScoresApply"))
setMethod("signatureScoresApply", "NanoStringRccSet",
function(X, MARGIN, FUN, ..., elt = "exprs")
{
  stopifnot(MARGIN %in% c(1L, 2L))
  if (MARGIN == 1L) {
    df <- data.frame()
    kvs <- c(sData(X), list(design = design(X)))
  } else {
    df <- sData(X)
    kvs <- list()
  }
  mat <- signatureScores(X, elt)
  .apply(X = mat, MARGIN = MARGIN, FUN = FUN, ..., .df = df, .kvs = kvs)
})

.apply <- function(X, MARGIN, FUN, ..., .df, .kvs)
{
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  environment(FUN) <- new.env(parent = parent)

  if (length(.kvs) > 0L) {
    multiassign(names(.kvs), .kvs, environment(FUN))
  }

  if (length(.df) == 0L) {
    apply(X, MARGIN = MARGIN, FUN = FUN, ...)
  } else {
    if (MARGIN == 1L) {
      output <- vector("list", nrow(X))
      for (i in seq_along(output)) {
        multiassign(colnames(.df), .df[i, ], environment(FUN))
        output[[i]] <- FUN(X[i, ], ...)
      }
      names(output) <- rownames(X)
    } else {
      output <- vector("list", ncol(X))
      for (j in seq_along(output)) {
        multiassign(colnames(.df), .df[j, ], environment(FUN))
        output[[j]] <- FUN(X[, j], ...)
      }
      names(output) <- colnames(X)
    }
    simplify2array(output, higher = FALSE)
  }
}


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
  keys <- sort(unique(values), na.last = TRUE)
  names(keys) <- as.character(keys)
  vapply(keys,
         function(k) {
           if (is.na(k))
             keep <- which(is.na(values))
           else
             keep <- which(!is.na(values) & values == k)
           if (GROUP == "featureData")
             FUN(X[keep, ], ...)
           else
             FUN(X[, keep], ...)
         }, simplify = simplify)
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
