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
         tooltip_digits = 4L)
{
  args <- list(...)
  pointColor <- "#4E79A7"
  refLineColor <- "#E15759"
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
         "heatmap-signatures" = {
           scores <- signatureScores(object, elt)
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
         "lane-bindingDensity" = {
           mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                 tooltip = "SampleName")
           p <- ggplot(object, mapping, ...) +
             geom_point_interactive(color = pointColor, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, 1.8, 2.25), linetype = 2L,
                        color = refLineColor)
         },
         "lane-fov" = {
           df <- pData(protocolData(object))
           df <- data.frame(FOVCounted = df[["FovCounted"]] / df[["FovCount"]],
                            row.names = rownames(df))
           mapping <- aes_string(x = "LaneID", y = "FOVCounted",
                                 tooltip = "SampleName")
           p <- ggplot(object, mapping, extradata = df, ...) +
             geom_point_interactive(color = pointColor, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_hline(yintercept = 0.75, linetype = 2L, color = refLineColor)
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
             aes_string(x = "ControlConc", y = elt, group = "SampleName")
           p <- ggplot(positiveControlSubset(object), mapping, ...) +
             geom_line_interactive(aes_string(tooltip = "SampleName"),
                                   color = pointColor) +
             geom_point_interactive(aes_string(tooltip = "SampleName"),
                                    color = pointColor, ...) +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
})
