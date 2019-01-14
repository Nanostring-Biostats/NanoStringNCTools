setGeneric("exprsplot", signature = "object",
           function(object, ...) standardGeneric("exprsplot"))

setMethod("exprsplot", "NanoStringRccSet",
function(object, ..., gene, group = NULL, log2scale = TRUE, elt = "exprs",
         points_overlay = TRUE, tooltip_digits = 4L)
{
  stopifnot(length(gene) == 1L &&
              (is.character(gene) ||
                 (is.numeric(gene) && (gene >= 1L & gene <= nrow(object)))))
  if (!is.null(group) && !all(group %in% svarLabels(object)))
    stop("all variables in \"group\" must be in 'svarLabels'")

  y <- assayDataElement2(object, elt)
  if (is.numeric(gene))
    gene <- featureData(object)[["GeneName"]][gene]
  y <- y[featureData(object)[["GeneName"]] == gene, ]
  if (is.null(group)) {
    x <- factor("")
    group <- ""
  } else if (length(group) == 1L)
    x <- factor(sData(object)[[group]])
  else {
    x <- do.call(interaction,
                 c(sData(object)[group],
                   list(drop = TRUE, sep = ":")))
    group <- paste(group, collapse = ":")
  }
  tooltip <- sprintf("%s<br>%s = %s<br>%s = %s", sampleNames(object),
                     group, x, gene, signif(y, tooltip_digits))
  df <- data.frame(group = x, gene = y, tooltip = tooltip,
                   stringsAsFactors = FALSE)
  p <- ggplot(df, aes(x = group, y = gene), ...) +
    geom_boxplot(outlier.shape = NA, ...) +
    scale_x_discrete(name = group)
  if (log2scale)
    p <- p + scale_y_continuous(name = gene, trans = "log2")
  else
    p <- p + scale_y_continuous(name = gene)
  if (points_overlay)
    p <- p + geom_beeswarm_interactive(aes(tooltip = tooltip), ...)
  p
})
