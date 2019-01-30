setGeneric("protoplot", signature = "object",
           function(object, ...) standardGeneric("protoplot"))

setMethod("protoplot", "NanoStringRccSet",
function(object, ...,
         type = c("boxplot-feature",
                  "boxplot-signature",
                  "bindingDensity-mean",
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
         index = 1L,
         tooltipHeading = NULL,
         tooltipDigits = 4L)
{
  args <- list(...)
  colors <- c(blue = "#4E79A7", orange = "#F28E2B", red = "#E15759",
              darkgray = "#79706E")
  type <- match.arg(type)
  switch(type,
         "boxplot-feature" =,
         "boxplot-signature" = {
           if (type == "boxplot-feature") {
             scores <- assayDataElement2(object, elt)
             name <- fData(object)[index, "GeneName"]
           } else {
             scores <- signatureScores(object, elt)
             name <- rownames(scores)[index]
           }
           if (is.null(group)) {
             x <- rep("", ncol(scores))
           } else {
             x <- sData(object)[[group]]
           }
           y <- scores[index, ]
           if (is.null(group))
             tooltip <- sprintf("%s<br>%s = %s", colnames(scores),
                                name, signif(y, tooltipDigits))
           else
             tooltip <- sprintf("%s<br>%s = %s<br>%s = %s", colnames(scores),
                                group, x, name, signif(y, tooltipDigits))
           df <- data.frame(group = x, signature = y, tooltip = tooltip,
                            stringsAsFactors = FALSE)
           p <- ggplot(df, aes_string(x = "group", y = "signature")) +
             stat_boxplot(geom = "errorbar", width = 0.5) +
             geom_boxplot_interactive(aes_string(tooltip = "group"),
                                      outlier.shape = NA) +
             scale_x_discrete(name = group) + scale_y_continuous(name = name) +
             geom_beeswarm_interactive(aes_string(tooltip = "tooltip"),
                                       size = 2, color = colors[["blue"]])
         },
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
           p <- protoheatmap(scores, log2scale = log2scale, group = group,
                             object = object,  tooltipHeading = tooltipHeading,
                             ...)
         },
         "heatmap-signatures" = {
           scores <- signatureScores(object, elt)
           p <- protoheatmap(scores, log2scale = log2scale, group = group,
                             object = object, tooltipHeading = tooltipHeading,
                             ...)
         },
         "lane-bindingDensity" = {
           maxBD <- 2.25
           outlier <- protocolData(object)[["BindingDensity"]] > maxBD
           if (any(outlier)) {
             extradata <- data.frame(Outlier = outlier,
                                     row.names = sampleNames(object))
             mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                   tooltip = "SampleName", color = "Outlier")
             p <- ggplot(object, mapping, extradata = extradata, ...) +
               geom_point_interactive() +
               scale_color_manual(values = unname(colors[c("blue", "red")]))
           } else {
             mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                   tooltip = "SampleName")
             p <- ggplot(object, mapping, ...) +
               geom_point_interactive(color = colors[["blue"]])
           }
           p <- p +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, maxBD), linetype = 2L,
                        color = colors[["darkgray"]])
         },
         "lane-fov" = {
           extradata <- pData(protocolData(object))
           extradata <-
             data.frame(FOVCounted =
                          extradata[["FovCounted"]] / extradata[["FovCount"]],
                        row.names = rownames(extradata))
           outlier <- extradata[["FOVCounted"]] < 0.75
           if (any(outlier)) {
             df[["Outlier"]] <- outlier
             mapping <- aes_string(x = "LaneID", y = "FOVCounted",
                                   tooltip = "SampleName", color = "Outlier")
             p <- ggplot(object, mapping, extradata = extradata, ...) +
               geom_point_interactive() +
               scale_color_manual(values = unname(colors[c("blue", "red")]))
           } else {
             mapping <- aes_string(x = "LaneID", y = "FOVCounted",
                                   tooltip = "SampleName")
             p <- ggplot(object, mapping, extradata = extradata, ...) +
               geom_point_interactive(color = colors[["blue"]])
           }
           p <- p +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_hline(yintercept = 0.75, linetype = 2L,
                        color = colors[["darkgray"]])
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
                     colnames(df)[1L], signif(df[,1L], tooltipDigits),
                     colnames(df)[2L], signif(df[,2L], tooltipDigits))
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
               scale_color_manual(values = unname(colors[c("blue", "red")]))
           } else {
             mapping <-
               aes_string(x = "ControlConc", y = elt, group = "SampleName",
                          tooltip = "CustomTooltip")
             p <- ggplot(posCtrl, mapping, extradata = extradata, ...) +
               geom_line_interactive(color = colors[["blue"]]) +
               geom_point_interactive(color = colors[["blue"]])
           }
           p <- p +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
})


protoheatmap <-
function(scores, log2scale, group, object,
         tooltipHeading = NULL,
         scaleCutoff = 3,
         groupPalette =
           c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
             "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
         ...)
{
  if (log2scale)
    scores <- log2t(scores)

  if (!is.null(tooltipHeading))
    colnames(scores) <- sData(object)[[tooltipHeading]]

  scaleCutoff <- abs(scaleCutoff)
  scores <- t(pmin(pmax(scale(scores), - scaleCutoff), scaleCutoff))

  if (is.null(group)) {
    annotation_row <- NA
    annotation_colors <- NA
  } else {
    annotation_row <- sData(object)[group]
    annotation_row[] <-
      lapply(annotation_row, function(x) {
        if (is.factor(x)) {
          levels(x) <- trimws(levels(x))
          x[x == ""] <- NA
          x <- x[drop = TRUE]
        } else {
          x <- trimws(x)
          x[x == ""] <- NA
        }
        x <- addNA(x, ifany = TRUE)
        levels(x)[is.na(levels(x))] <- "<N/A>"
        x
      })
    rownames(annotation_row) <- rownames(scores)

    annotation_colors <- cumsum(sapply(annotation_row, nlevels))
    annotation_colors <- Map(`:`, c(1L, head(annotation_colors, -1L) + 1L),
                             annotation_colors)
    annotation_colors <- structure(lapply(annotation_colors, function(x) {
      x <- x %% length(groupPalette)
      x[x == 0L] <- length(groupPalette)
      groupPalette[x]
    }), names = colnames(annotation_row))
    for (j in seq_len(ncol(annotation_row))) {
      names(annotation_colors[[j]]) <- levels(annotation_row[[j]])
    }
  }

  pheatmap(scores,
           color =
             colorRampPalette(c("darkblue",
                                rev(brewer.pal(n = 7L, name = "RdYlBu")),
                                "darkred"))(100),
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           show_rownames = (nrow(scores) <= 64L),
           show_colnames = (ncol(scores) <= 64L),
           silent = TRUE,
           ...)
}
