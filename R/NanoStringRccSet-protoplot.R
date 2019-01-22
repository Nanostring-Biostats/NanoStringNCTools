setGeneric("protoplot", signature = "object",
           function(object, ...) standardGeneric("protoplot"))

setMethod("protoplot", "NanoStringRccSet",
function(object, ...,
         type = c("bindingDensity-mean",
                  "bindingDensity-sd",
                  "heatmap-genes",
                  "heatmap-signatures",
                  "lane-bindingDensity",
                  "lane-fov",
                  "mean-sd-features",
                  "mean-sd-samples",
                  "positiveControl"),
         log2scale = TRUE,
         elt = "exprs",
         group = NULL,
         tooltip_digits = 4L)
{
  args <- list(...)
  blue <- "#4E79A7"
  orange <- "#F28E2B"
  red <- "#E15759"
  type <- match.arg(type)
  switch(type,
         "bindingDensity-mean" =,
         "bindingDensity-sd" = {
           if (log2scale) {
             y <- if (type == "bindingDensity-mean") "MeanLog2" else "SDLog2"
           } else {
             y <- if (type == "bindingDensity-mean") "Mean" else "SD"
           }
           mapping <- aes_string(x = "BindingDensity", y = y,
                                 tooltip = "SampleName")
           p <- ggplot(object, mapping, ...) + geom_point_interactive(...) +
             scale_x_continuous(name = "Binding Density")
         },
         "heatmap-genes" = {
           object <- endogenousSubset(object)
           scores <- assayDataElement2(object, elt)
           rownames(scores) <- featureData(object)[["GeneName"]]
           colnames(scores) <- protocolData(object)[["SampleID"]]
           p <- protoheatmap(scores, log2scale = log2scale, group = group,
                             object = object, ...)
         },
         "heatmap-signatures" = {
           scores <- signatureScores(object, elt)
           colnames(scores) <- protocolData(object)[["SampleID"]]
           p <- protoheatmap(scores, log2scale = log2scale, group = group,
                             object = object, ...)
         },
         "lane-bindingDensity" = {
           mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                 tooltip = "SampleName")
           p <- ggplot(object, mapping, ...) +
             geom_point_interactive(color = blue, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, 1.8, 2.25), linetype = 2L,
                        color = red)
         },
         "lane-fov" = {
           df <- pData(protocolData(object))
           df <- data.frame(FOVCounted = df[["FovCounted"]] / df[["FovCount"]],
                            row.names = rownames(df))
           mapping <- aes_string(x = "LaneID", y = "FOVCounted",
                                 tooltip = "SampleName")
           p <- ggplot(object, mapping, extradata = df, ...) +
             geom_point_interactive(color = blue, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_hline(yintercept = 0.75, linetype = 2L, color = red)
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
           posCtrl <- positiveControlSubset(object)
           posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
           x <- log2(featureData(posCtrl)[["ControlConc"]])
           extradata <-
             data.frame(RSquared =
                          esApply(posCtrl, 2L,
                                  function(y) cor(x, log2t(y, 0.5))))
           extradata[["Low R-Squared"]] <- extradata[["RSquared"]] < 0.95
           extradata[["CustomTooltip"]] <-
             sprintf("%s\nR-Squared = %.4f", rownames(extradata),
                     extradata[["RSquared"]])

           if (any(extradata[["Low R-Squared"]])) {
             mapping <-
               aes_string(x = "ControlConc", y = elt, group = "SampleName",
                          tooltip = "CustomTooltip",
                          color = as.name("Low R-Squared"))
             p <- ggplot(posCtrl, mapping, extradata = extradata, ...) +
               geom_line_interactive() +
               geom_point_interactive() +
               scale_color_manual(values = c(blue, red))
           } else {
             mapping <-
               aes_string(x = "ControlConc", y = elt, group = "SampleName",
                          tooltip = "CustomTooltip")
             p <- ggplot(posCtrl, mapping, extradata = extradata, ...) +
               geom_line_interactive(color = blue) +
               geom_point_interactive(color = blue)
           }
           p <- p +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
})


protoheatmap <- function(scores, log2scale, group, object, ...)
{
  if (log2scale)
    scores <- log2t(scores)
  scores <- t(pmin(pmax(scale(scores), -3), 3))

  if (is.null(group)) {
    annotation_row <- NA
  } else {
    annotation_row <- sData(object)[group]
    rownames(annotation_row) <- rownames(scores)
  }

  pheatmap(scores,
           color =
             colorRampPalette(c("darkblue",
                                rev(brewer.pal(n = 7L, name = "RdYlBu")),
                                "darkred"))(100),
           annotation_row = annotation_row,
           show_rownames = (nrow(scores) <= 64L),
           show_colnames = (ncol(scores) <= 64L),
           silent = TRUE,
           ...)
}
