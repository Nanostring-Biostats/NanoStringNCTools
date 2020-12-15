setGeneric("munge", signature = "data",
           function(data, ...) standardGeneric("munge"))

setMethod("munge", "NanoStringRccSet",
function(data, mapping = update(design(data), exprs ~ .), extradata = NULL,
         elt = "exprs", ...)
{
  # Get list of variables
  mapping <- try(mapping, silent = TRUE)
  if (inherits(mapping, "try-error"))
    stop("\"mapping\" argument is missing")
  if (inherits(mapping, "formula"))
    vars <- all.vars(mapping)
  else if (is.list(mapping))
    vars <- unique(unlist(lapply(mapping, all.vars), use.names = FALSE))

  # Handle GeneMatrix / SignatureMatrix special case
  hasGeneMatrix <- "GeneMatrix" %in% vars
  hasSignatureMatrix <- "SignatureMatrix" %in% vars
  if (hasGeneMatrix || hasSignatureMatrix) {
    df <- data.frame(SampleName = sampleNames(data),
                     sData(data),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
    if (!is.null(extradata) && any(vars %in% colnames(extradata))) {
      df <- cbind(df, extradata)
    }
    sampleLabels <- sData(data)[[dimLabels(data)[2L]]]
    if (hasGeneMatrix) {
      mat <- assayDataElement2(data, elt)
      rownames(mat) <- featureData(data)[["GeneName"]]
      colnames(mat) <- sampleLabels
      df[["GeneMatrix"]] <- t(mat)
    }
    if (hasSignatureMatrix) {
      mat <- signatureScores(data, elt)
      colnames(mat) <- sampleLabels
      df[["SignatureMatrix"]] <- t(mat)
    }
    if (!all(vars %in% colnames(df))) {
      stop("\"mapping\" contains undefined variables")
    }
    return(as(df[vars], "DataFrame"))
  }

  # Determine the types of variables
  hasFeatureVars <- any(vars %in% fvarLabels(data))
  hasSampleVars  <- any(vars %in% svarLabels(data))
  hasLog2Summaries <- any(vars %in% rownames(.summaryMetadata[["log2"]]))
  hasSummaries <- any(vars %in% rownames(.summaryMetadata[["moments"]]))
  hasQuantiles <- any(vars %in% rownames(.summaryMetadata[["quantiles"]]))
  if (hasQuantiles && !hasLog2Summaries)
    hasSummaries <- TRUE
  hasAggregates <- hasLog2Summaries || hasSummaries
  hasAssayDataElts <- any(vars %in% assayDataElementNames(data))
  useSignatures <- "SignatureName" %in% vars
  if (hasAggregates) {
    hasFeatureVars <- hasFeatureVars || ("FeatureName" %in% vars)
    hasSampleVars  <- hasSampleVars  || ("SampleName"  %in% vars)
    if (hasAssayDataElts)
      stop("\"mapping\" argument cannot contain both aggregates and disaggregates")
    if (hasFeatureVars && hasSampleVars)
      stop("\"mapping\" argument cannot aggregate using both feature and sample variables")
    if (useSignatures && hasFeatureVars)
      stop("\"mapping\" argument cannot aggregate using both signatures and feature variables")
    if (useSignatures && hasSampleVars)
      stop("\"mapping\" argument cannot aggregate using both signatures and sample variables")
    if (!hasFeatureVars && !hasSampleVars && !useSignatures)
      stop("\"mapping\" argument contains an ambiguous aggregation")
  }
  if (hasLog2Summaries && hasSummaries)
    stop("\"mapping\" argument cannot use both log2 and linear aggregations")
  if (useSignatures && hasFeatureVars)
    stop("\"mapping\" argument cannot use both signatures and feature variables")
  if ((!hasAggregates && !hasAssayDataElts) &&
      (hasFeatureVars || useSignatures) && hasSampleVars)
    stop("\"mapping\" argument contains an ambiguous variable selection")

  copyRowNames <- function(df, key) {
    rn <- data.frame(rownames(df), stringsAsFactors = FALSE)
    colnames(rn) <- key
    cbind(rn, df)
  }

  # Produce either disaggregate or aggregate results
  if (hasAssayDataElts) {
    assayDataElts <- intersect(vars, assayDataElementNames(data))
    if (useSignatures) {
      df <-
        vapply(assayDataElts, function(elt) {
          as.vector(signatureScores(data, elt))
        })
      df <-
        data.frame(SignatureName = rep.int(names(signatures(data)), ncol(data)),
                   SampleName = rep(sampleNames(data),
                                    each = length(signatures(data))),
                   df,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
    } else {
      df <-
        vapply(assayDataElts, function(elt) {
          as.vector(assayDataElement2(data, elt))
        })
      df <-
        data.frame(FeatureName = rep.int(featureNames(data), ncol(data)),
                   SampleName = rep(sampleNames(data), each = nrow(data)),
                   df,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
    }
  } else if (hasAggregates) {
    # Calculate marginal summaries
    if (useSignatures) {
      df <- summary(data, MARGIN = 1L, log2scale = hasLog2Summaries,
                    elt = elt, signatureScores = TRUE)
      df <- df[, intersect(vars, colnames(df)), drop = FALSE]
      df <- copyRowNames(df, "SignatureName")
    } else {
      MARGIN <- 1L + hasSampleVars
      df <- summary(data, MARGIN = MARGIN, log2scale = hasLog2Summaries,
                    elt = elt)
      df <- df[, intersect(vars, colnames(df)), drop = FALSE]
      if (MARGIN == 1L) {
        df <- copyRowNames(df, "FeatureName")
      } else {
        df <- copyRowNames(df, "SampleName")
      }
    }
  } else {
    df <- NULL
  }

  # Add feature variables, if requested
  if (hasFeatureVars) {
    fvars <- intersect(vars, fvarLabels(data))
    fdf <- fData(data)[, fvars, drop = FALSE]
    if (is.null(df)) {
      df <- copyRowNames(fdf, "FeatureName")
    } else {
      df <- cbind(df, fdf[df[["FeatureName"]], , drop = FALSE])
    }
  }

  if (useSignatures && is.null(df)) {
    df <- data.frame(SignatureName = names(signatures(data)),
                     stringsAsFactors = FALSE)
  }

  # Add sample variables, if requested
  if (hasSampleVars) {
    svars <- intersect(vars, svarLabels(data))
    sdf <- sData(data)[, svars, drop = FALSE]
    if (is.null(df)) {
      df <- copyRowNames(sdf, "SampleName")
    } else {
      df <- cbind(df, sdf[df[["SampleName"]], , drop = FALSE])
    }
  }

  # Add extra data, if supplied
  if (!is.null(extradata)) {
    matchFeatureNames <- identical(rownames(extradata), featureNames(data))
    matchSampleNames  <- identical(rownames(extradata), sampleNames(data))
    if (!matchFeatureNames && !matchSampleNames) {
      stop("\"extradata\" 'rownames' do not match 'featureNames' or 'sampleNames'")
    }
    evars <- intersect(vars, colnames(extradata))
    edf <- extradata[, evars, drop = FALSE]
    if (matchFeatureNames) {
      if (is.null(df)) {
        df <- copyRowNames(edf, "FeatureName")
      } else {
        df <- cbind(df, edf[df[["FeatureName"]], , drop = FALSE])
      }
    } else {
      if (is.null(df)) {
        df <- copyRowNames(edf, "SampleName")
      } else {
        df <- cbind(df, edf[df[["SampleName"]], , drop = FALSE])
      }
    }
  }

  # Check that all columns are present
  if (!all(vars %in% colnames(df))) {
    stop("\"mapping\" contains undefined variables")
  }

  # Return result
  rownames(df) <- NULL
  df
})
