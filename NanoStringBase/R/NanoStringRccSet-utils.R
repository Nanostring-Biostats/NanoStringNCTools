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


# Plotting
setMethod("mold", "NanoStringRccSet",
function(data, mapping = design(data), extradata = NULL, ...)
{
  if (is.null(mapping))
    stop("\"mapping\" argument is missing")
  if (inherits(mapping, "formula"))
    vars <- all.vars(mapping)
  else if (is.list(mapping))
    vars <- unique(unlist(lapply(mapping, all.vars), use.names = FALSE))
  hasFeatureVars <- any(vars %in% fvarLabels(data))
  hasSampleVars  <- any(vars %in% svarLabels(data))
  if (hasFeatureVars && hasSampleVars)
    stop("\"mapping\" argument cannot use both feature and sample variables")
  if (hasFeatureVars)
    df <- fData(data)
  else if (hasSampleVars)
    df <- sData(data)
  else
    df <- NULL
  if (!is.null(extradata)) {
    if (!identical(rownames(extradata), featureNames(data)) &&
        !identical(rownames(extradata), sampleNames(data))) {
      stop("\"extradata\" 'rownames' do not match 'featureNames' or 'sampleNames'")
    }
    if (is.null(df))
      df <- extradata
    else
      df <- cbind(df, extradata)
  }
  assayDataElts <- intersect(vars, assayDataElementNames(data))
  if (length(assayDataElts) == 0L) {
    df[, vars, drop = FALSE]
  } else {
    df <- df[, setdiff(vars, assayDataElts), drop = FALSE]
    transpose <- identical(rownames(df), sampleNames(data))
    stackedData <-
      sapply(assayDataElts, function(elt) {
        mat <- assayDataElement(data, elt)
        if (transpose)
          mat <- t(mat)
        as.vector(mat)
      })
    if (transpose) {
      stackedData <-
        data.frame(FeatureName = rep(featureNames(data), each = ncol(data)),
                   SampleName = rep.int(sampleNames(data), nrow(data)),
                   stackedData)
      df <- df[stackedData[["SampleName"]], , drop = FALSE]
    } else {
      stackedData <-
        data.frame(FeatureName = rep.int(featureNames(data), ncol(data)),
                   SampleName = rep(sampleNames(data), each = nrow(data)),
                   stackedData)
      df <- df[stackedData[["FeatureName"]], , drop = FALSE]
    }
    cbind(stackedData, df)
  }
})

ggplot.NanoStringRccSet <-
function(data, mapping = aes(), extradata = NULL, ...,
         environment = parent.frame())
{
  if (length(mapping) == 0L) {
    mapping <- design(data)
    if (is.null(mapping))
      stop("\"mapping\" argument is missing")
  }
  df <- mold(data, mapping = mapping, extradata = extradata)
  g <- ggplot(df, mapping, ..., environment = environment)
  GGbio(g, data = data)
}

setMethod("autoplot", "NanoStringRccSet",
function(object, ...,
         type = c("meanlog2-sdlog2-features",
                  "meanlog2-sdlog2-samples",
                  "mean-sd-features",
                  "mean-sd-samples"),
         elt = "exprs",
         tooltip_digits = 6L)
{
  args <- list(...)
  type <- match.arg(type)
  switch(type,
         "meanlog2-sdlog2-features" =,
         "meanlog2-sdlog2-samples" =,
         "mean-sd-features" =,
         "mean-sd-samples" = {
           MARGIN <- 1L + (type %in% c("meanlog2-sdlog2-samples",
                                       "mean-sd-samples"))
           log2scale <- type %in% c("meanlog2-sdlog2-features",
                                    "meanlog2-sdlog2-samples")
           stats <- summary(object, MARGIN = MARGIN, log2scale = log2scale,
                            elt = elt)
           if (log2scale) {
             mapping <- aes_string(x = "MeanLog2", y = "SDLog2",
                                   tooltip = "ToolTip")
             df <- as.data.frame(stats[,c("MeanLog2", "SDLog2")])
           } else {
             mapping <- aes_string(x = "Mean", y = "SD", tooltip = "ToolTip")
             df <- as.data.frame(stats[,c("Mean", "SD")])
           }
           df[["ToolTip"]] <-
             sprintf("%s<br>%s = %s<br>%s = %s", rownames(df),
                     colnames(df)[1L], signif(df[,1L], tooltip_digits),
                     colnames(df)[2L], signif(df[,2L], tooltip_digits))
           p <- ggplot(df, mapping, ...) + geom_point_interactive(...)
         })

  p
})
