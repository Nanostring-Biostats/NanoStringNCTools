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
  vars <- setdiff(vars, c("FeatureName", "SampleName"))

  # Determine the types of variables
  hasFeatureVars <- any(vars %in% fvarLabels(data))
  hasSampleVars  <- any(vars %in% svarLabels(data))
  hasLog2Summaries <- any(vars %in% rownames(.summaryMetadata[["log2"]]))
  hasSummaries <- any(vars %in% rownames(.summaryMetadata[["moments"]]))
  hasQuantiles <- any(vars %in% rownames(.summaryMetadata[["quantiles"]]))
  if (hasQuantiles && !hasLog2Summaries)
    hasSummaries <- TRUE
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
    vars <- setdiff(vars, assayDataElts)
  }
  if (!all(vars %in% colnames(df))) {
    stop("\"mapping\" contains undefined variables")
  }
  df <- df[, vars, drop = FALSE]

  # Get assay data elements if needed
  if (length(assayDataElts) == 0L) {
    if (identical(rownames(df), featureNames(data))) {
      df <- cbind(data.frame(FeatureName = rownames(df),
                             stringsAsFactors = FALSE), df)
    } else {
      df <- cbind(data.frame(SampleName = rownames(df),
                             stringsAsFactors = FALSE), df)
    }
    rownames(df) <- NULL
  } else {
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
function(data, mapping = aes(), ..., extradata = NULL,
         tooltip_digits = 4L, environment = parent.frame())
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
          sprintf("%s<br>%s = %s", df[[tooltip]], names(mf)[1L],
                  signif(mf[[1L]], digits = tooltip_digits))
      }
    }
  }
  g <- ggplot(df, mapping, ..., environment = environment)
  GGbio(g, data = data)
}


geom_beeswarm_interactive <-
function(mapping = NULL, data = NULL,
         priority = c("ascending", "descending", "density", "random", "none"),
         cex = 1, groupOnX = NULL, dodge.width = 0, stat = "identity",
         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  position <- position_beeswarm(priority = priority, cex = cex,
                                groupOnX = groupOnX, dodge.width = dodge.width)
  layer(data = data, mapping = mapping, stat = stat,
        geom = GeomInteractivePoint, position = position,
        show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}


format_percent <- function(x) {
  sprintf("%s%%", format(100 * x, digits = 2L))
}
