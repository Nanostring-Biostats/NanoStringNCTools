setMethod("mold", "NanoStringRccSet",
function(data, mapping = design(data), extradata = NULL, elt = "exprs", ...)
{
  # Get list of variables
  if (is.null(mapping))
    stop("\"mapping\" argument is missing")
  if (inherits(mapping, "formula"))
    vars <- all.vars(mapping)
  else if (is.list(mapping))
    vars <- unique(unlist(lapply(mapping, all.vars), use.names = FALSE))

  # Determine the types of variables
  hasFeatureVars <- any(vars %in% fvarLabels(data))
  hasSampleVars  <- any(vars %in% svarLabels(data))
  hasLog2Summaries <- any(vars %in% c("GeomMean", "MeanLog2", "SDLog2"))
  hasSummaries <- any(vars %in% c("Mean", "SD", "Skewness", "Kurtosis"))
  hasQuantiles <- any(vars %in% c("Min", "Q1", "Median", "Q3", "Max"))
  if (hasQuantiles && !hasSummaries)
    hasLog2Summaries <- TRUE
  if (hasFeatureVars && hasSampleVars)
    stop("\"mapping\" argument cannot use both feature and sample variables")
  if (hasLog2Summaries && hasSummaries)
    stop("\"mapping\" argument cannot use both log2 and linear summary statistics")

  # Get marginal data if needed
  if (hasFeatureVars)
    df <- fData(data)
  else if (hasSampleVars)
    df <- sData(data)
  else
    df <- NULL

  # Add extra data if supplied
  if (!is.null(extradata)) {
    matchFeatureNames <- identical(rownames(extradata), featureNames(data))
    matchSampleNames  <- identical(rownames(extradata), sampleNames(data))
    if (!matchFeatureNames && !matchSampleNames) {
      stop("\"extradata\" 'rownames' do not match 'featureNames' or 'sampleNames'")
    }
    if (is.null(df)) {
      hasFeatureVars <- matchFeatureNames
      hasSampleVars  <- matchSampleNames
      df <- extradata
    } else
      df <- cbind(df, extradata)
  }

  # Add marginal summaries if needed
  if (hasLog2Summaries || hasSummaries) {
    if (!hasFeatureVars && !hasSampleVars)
      stop("\"mapping\" argument contains ambiguous summary statistics")
    MARGIN <- 1L + hasSampleVars
    if (hasLog2Summaries)
      df <- cbind(df, summary(data, MARGIN = MARGIN, elt = elt))
    else
      df <- cbind(df, summary(data, MARGIN = MARGIN, log2scale = FALSE,
                              elt = elt))
  }

  # Determine if assay data elements are needed
  assayDataElts <- intersect(vars, assayDataElementNames(data))
  if (length(assayDataElts) > 0L) {
    vars <- setdiff(vars, c(assayDataElts, "FeatureName", "SampleName"))
  }
  if (!all(vars %in% colnames(df))) {
    stop("\"mapping\" contains undefined variables")
  }
  df <- df[, vars, drop = FALSE]

  # Get assay data elements if needed
  if (length(assayDataElts) > 0L) {
    df <- df[, setdiff(vars, c(assayDataElts, "FeatureName", "SampleName")),
             drop = FALSE]
    transpose <- identical(rownames(df), sampleNames(data))
    stackedData <-
      sapply(assayDataElts, function(elt) {
        mat <- assayDataElement2(data, elt)
        if (transpose)
          mat <- t(mat)
        as.vector(mat)
      })
    if (transpose) {
      stackedData <-
        data.frame(FeatureName = rep(featureNames(data), each = ncol(data)),
                   SampleName = rep.int(sampleNames(data), nrow(data)),
                   stackedData,
                   stringsAsFactors = FALSE)
      df <- df[stackedData[["SampleName"]], , drop = FALSE]
    } else {
      stackedData <-
        data.frame(FeatureName = rep.int(featureNames(data), ncol(data)),
                   SampleName = rep(sampleNames(data), each = nrow(data)),
                   stackedData,
                   stringsAsFactors = FALSE)
      df <- df[stackedData[["FeatureName"]], , drop = FALSE]
    }
    rownames(df) <- NULL
    df <- cbind(stackedData, df)
  }

  # Return result
  df
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
  if ("tooltip" %in% names(mapping)) {
    tooltip <- as.character(mapping[["tooltip"]][[2L]])
    for (j in c("x", "y")) {
      if (j %in% names(mapping)) {
        mf <- model.frame(mapping[[j]], df)
        df[[tooltip]] <-
          sprintf("%s<br>%s = %s", df[[tooltip]], names(mf)[1L], mf[[1L]])
      }
    }
  }
  g <- ggplot(df, mapping, ..., environment = environment)
  GGbio(g, data = data)
}


setMethod("autoplot", "NanoStringRccSet",
function(object, ...,
         type = c("bindingDensity-mean",
                  "bindingDensity-sd",
                  "fov-mean",
                  "fov-sd",
                  "heatmap",
                  "mean-sd-features",
                  "mean-sd-samples",
                  "positiveControl"),
         log2scale = TRUE,
         elt = "exprs",
         tooltip_digits = 6L)
{
  args <- list(...)
  type <- match.arg(type)
  switch(type,
         "bindingDensity-mean" =,
         "bindingDensity-sd" = {
           stats <- summary(object, log2scale = log2scale, elt = elt)
           df <- pData(protocolData(object))[, "BindingDensity", drop = FALSE]
           if (log2scale) {
             y <- if (type == "bindingDensity-mean") "MeanLog2" else "SDLog2"
           } else {
             y <- if (type == "bindingDensity-mean") "Mean" else "SD"
           }
           mapping <- aes_string(x = "BindingDensity", y = y,
                                 tooltip = "ToolTip")
           df <- cbind(df, stats[, y, drop = FALSE])
           df[["ToolTip"]] <-
             sprintf("%s<br>%s = %s<br>%s = %s", rownames(df),
                     colnames(df)[1L], signif(df[,1L], tooltip_digits),
                     colnames(df)[2L], signif(df[,2L], tooltip_digits))
           p <- ggplot(df, mapping, ...) + geom_point_interactive(...)
         },
         "fov-mean" =,
         "fov-sd" = {
           stats <- summary(object, log2scale = log2scale, elt = elt)
           df <- pData(protocolData(object))[, c("FovCounted", "FovCount")]
           df <- data.frame(FovCounted = df[["FovCounted"]] / df[["FovCount"]],
                            row.names = rownames(df))
           if (log2scale) {
             y <- if (type == "fov-mean") "MeanLog2" else "SDLog2"
           } else {
             y <- if (type == "fov-mean") "Mean" else "SD"
           }
           mapping <- aes_string(x = "FovCounted", y = y, tooltip = "ToolTip")
           df <- cbind(df, stats[, y, drop = FALSE])
           df[["ToolTip"]] <-
             sprintf("%s<br>%s = %s<br>%s = %s", rownames(df),
                     colnames(df)[1L], signif(df[,1L], tooltip_digits),
                     colnames(df)[2L], signif(df[,2L], tooltip_digits))
           p <- ggplot(df, mapping, ...) + geom_point_interactive(...) +
             scale_x_continuous(name = "FOV Counted", labels = percent)
         },
         "heatmap" = {
           object <- endogenousSubset(object)
           scores <- assayDataElement2(object, elt)
           rownames(scores) <- featureData(object)[["GeneName"]]
           colnames(scores) <- protocolData(object)[["SampleID"]]
           if (log2scale)
             scores <- log2t(scores)
           scores <- t(pmin(pmax(scale(scores), -3), 3))
           p <- pheatmap(scores,
                         color =
                           colorRampPalette(c("darkblue",
                                              rev(brewer.pal(n = 7L,
                                                             name = "RdYlBu")),
                                              "darkred"))(100),
                         show_rownames = (nrow(scores) <= 64L),
                         show_colnames = (ncol(scores) <= 64L),
                         silent = TRUE,
                         ...)
         },
         "mean-sd-features" =,
         "mean-sd-samples" = {
           MARGIN <- 1L + (type == "mean-sd-samples")
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
         },
         "positiveControl" = {
           mapping <-
             aes_string(x = "ControlConc", y = elt, tooltip = "SampleName")
           p <- ggplot(positiveControlSubset(object), mapping) +
             geom_point_interactive(...) +
             scale_x_continuous(trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
})
