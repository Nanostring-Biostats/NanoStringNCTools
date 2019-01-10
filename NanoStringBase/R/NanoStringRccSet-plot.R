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
  if (length(assayDataElts) > 0L) {
    vars <- setdiff(vars, c(assayDataElts, "FeatureName", "SampleName"))
  }
  if (!all(vars %in% colnames(df))) {
    stop("\"mapping\" contains undefined variables")
  }
  df <- df[, vars, drop = FALSE]
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
         type = c("MeanLog2-SDLog2-Features",
                  "MeanLog2-SDLog2-Samples",
                  "Mean-SD-Features",
                  "Mean-SD-Samples",
                  "PositiveControl-LogLog"),
         elt = "exprs",
         tooltip_digits = 6L)
{
  args <- list(...)
  type <- match.arg(type)
  switch(type,
         "MeanLog2-SDLog2-Features" =,
         "MeanLog2-SDLog2-Samples" =,
         "Mean-SD-Features" =,
         "Mean-SD-Samples" = {
           MARGIN <- 1L + (type %in% c("MeanLog2-SDLog2-Samples",
                                       "Mean-SD-Samples"))
           log2scale <- type %in% c("MeanLog2-SDLog2-Features",
                                    "MeanLog2-SDLog2-Samples")
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
         "PositiveControl-LogLog" = {
           mapping <-
             aes_string(x = "ControlConc", y = elt, tooltip = "SampleName")
           p <- ggplot(positiveControlSubset(object), mapping) +
             geom_point_interactive(...) +
             scale_x_continuous(trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
})
