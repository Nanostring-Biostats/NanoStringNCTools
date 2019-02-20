autoplot.NanoStringRccSet <-
function(object,
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
         index = 1L,
         geom = list(),
         tooltipHeading = NULL,
         tooltipDigits = 4L,
         heatmapGroup = NULL,
         ...)
{
  processGeom <- function(lst, nm, default_aes = aes()) {
    if (nm %in% names(lst)) {
      if ("color" %in% names(lst[[nm]])) {
        lst[[nm]][["colour"]] <- lst[[nm]][["color"]]
        lst[[nm]][["color"]] <- NULL
      }
    } else {
      lst[[nm]] <- aes()
    }
    if (length(default_aes) > 0L) {
      args <- setdiff(names(default_aes), c("tooltip", "onclick", "data_id"))
      lst[[nm]] <- lst[[nm]][names(lst[[nm]]) %in% args]
      lst[[nm]] <- c(lst[[nm]], default_aes[setdiff(args, names(lst[[nm]]))])
      oldClass(lst[[nm]]) <- "uneval"
    }
    lst
  }
  geom <- as.list(geom)
  geom <- processGeom(geom, "base")
  geom <- processGeom(geom, "point", GeomInteractivePoint$default_aes)
  geom <- processGeom(geom, "line", GeomInteractiveLine$default_aes)
  geom <- processGeom(geom, "boxplot", GeomInteractiveBoxplot$default_aes)

  ggpoint <- function(mapping, ...) {
    for (arg in names(geom[["point"]])) {
      if (is.name(geom[["point"]][[arg]])) {
        mapping[[arg]] <- geom[["point"]][[arg]]
        geom[["point"]][[arg]] <- NULL
      }
    }
    ggplot(object, mapping, elt = elt, ...) +
      do.call(geom_point_interactive, geom[["point"]])
  }

  ggline <- function(object, mapping, ...) {
    for (i in c("line", "point")) {
      for (arg in names(geom[[i]])) {
        if (is.name(geom[[i]][[arg]])) {
          mapping[[arg]] <- geom[[i]][[arg]]
          geom[[i]][[arg]] <- NULL
        }
      }
    }
    ggplot(object, mapping, elt = elt, ...) +
      do.call(geom_line_interactive, geom[["line"]]) +
      do.call(geom_point_interactive, geom[["point"]])
  }

  type <- match.arg(type)
  switch(type,
         "boxplot-feature" =,
         "boxplot-signature" = {
           if (type == "boxplot-feature") {
             scores <- assayDataElement2(object, elt)
             ytitle <- fData(object)[index, "GeneName"]
           } else {
             scores <- signatureScores(object, elt)
             ytitle <- rownames(scores)[index]
           }
           y <- scores[index, ]
           if (!is.name(geom[["base"]][["x"]])) {
             x <- rep.int("", ncol(scores))
             xtitle <- ""
           } else {
             x <- eval(geom[["base"]][["x"]], sData(object))
             xtitle <- as.character(geom[["base"]][["x"]])
           }
           if (!is.name(geom[["point"]][["color"]])) {
             color <- NULL
             colortitle <- ""
           } else {
             color <- eval(geom[["point"]][["color"]], sData(object))
             colortitle <- as.character(geom[["point"]][["color"]])
           }
           tooltip <- colnames(scores)
           if (is.name(geom[["base"]][["x"]])) {
             tooltip <- sprintf("%s<br>%s = %s", tooltip, xtitle, x)
           }
           tooltip <- sprintf("%s<br>%s = %s", tooltip, ytitle,
                              signif(y, tooltipDigits))
           if (is.name(geom[["point"]][["color"]])) {
             tooltip <- sprintf("%s<br>%s = %s", tooltip, colortitle, color)
           }
           df <- data.frame(x = x, score = y, tooltip = tooltip,
                            stringsAsFactors = FALSE)
           df[["color"]] <- color
           p <- ggplot(df, aes_string(x = "x", y = "score")) +
             stat_boxplot(geom = "errorbar",
                          width = geom[["boxplot"]][["size"]],
                          color = geom[["boxplot"]][["colour"]]) +
             do.call(geom_boxplot_interactive,
                     c(list(aes_string(tooltip = "x")),
                       geom[["boxplot"]],
                       outlier.shape = NA)) +
             scale_x_discrete(name = xtitle) + scale_y_continuous(name = ytitle)
           if (is.null(color)) {
             p <- p +
               do.call(geom_beeswarm_interactive,
                       c(list(aes_string(tooltip = "tooltip")),
                         geom[["point"]]))
           } else {
             p <- p +
               do.call(geom_beeswarm_interactive,
                       c(list(aes_string(color = "color", tooltip = "tooltip")),
                         geom[["point"]])) +
               guides(color = guide_legend(title = colortitle))
           }
           p
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
           p <- ggpoint(mapping, ...) +
             scale_x_continuous(name = "Binding Density")
         },
         "heatmap-genes" = {
           object <- endogenousSubset(object)
           scores <- assayDataElement2(object, elt)
           rownames(scores) <- featureData(object)[["GeneName"]]
           p <-
             protoheatmap(scores, log2scale = log2scale, group = heatmapGroup,
                          object = object,  tooltipHeading = tooltipHeading,
                          ...)
         },
         "heatmap-signatures" = {
           scores <- signatureScores(object, elt)
           if (anyNA(scores)) {
             whichNA <- which(is.na(rowSums(scores)))
             warning(sprintf("dropped %d signatures due to missing values",
                             length(whichNA)))
             scores <- scores[- whichNA, , drop = FALSE]
           }
           p <-
             protoheatmap(scores, log2scale = log2scale, group = heatmapGroup,
                          object = object, tooltipHeading = tooltipHeading,
                          ...)
         },
         "lane-bindingDensity" = {
           maxBD <- 2.25
           extradata <-
             data.frame(Outlier =
                          protocolData(object)[["BindingDensity"]] > maxBD,
                        row.names = sampleNames(object))
           mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                 tooltip = "SampleName")
           if (any(extradata[["Outlier"]])) {
             geom[["point"]] <- unclass(geom[["point"]])
             if (!is.name(geom[["point"]][["colour"]])) {
               mapping[["colour"]] <- as.name("Outlier")
               geom[["point"]][["colour"]] <- NULL
             } else if (!is.name(geom[["point"]][["shape"]])) {
               mapping[["shape"]] <- as.name("Outlier")
               geom[["point"]][["shape"]] <- NULL
             }
             oldClass(geom[["point"]]) <- "uneval"
           }
           p <- ggpoint(mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, maxBD), linetype = 2L,
                        color = "darkgray")
         },
         "lane-fov" = {
           extradata <- pData(protocolData(object))
           extradata <-
             data.frame(FOVCounted =
                          extradata[["FovCounted"]] / extradata[["FovCount"]],
                        row.names = rownames(extradata))
           extradata[["Outlier"]] <- extradata[["FOVCounted"]] < 0.75
           mapping <- aes_string(x = "LaneID", y = "FOVCounted",
                                 tooltip = "SampleName")
           if (any(extradata[["Outlier"]])) {
             geom[["point"]] <- unclass(geom[["point"]])
             if (!is.name(geom[["point"]][["colour"]])) {
               mapping[["colour"]] <- as.name("Outlier")
               geom[["point"]][["colour"]] <- NULL
             } else if (!is.name(geom[["point"]][["shape"]])) {
               mapping[["shape"]] <- as.name("Outlier")
               geom[["point"]][["shape"]] <- NULL
             }
             oldClass(geom[["point"]]) <- "uneval"
           }
           p <- ggpoint(mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_hline(yintercept = 0.75, linetype = 2L,
                        color = "darkgray")
         },
         "mean-sd-features" = {
           if (log2scale)
             mapping <- aes_string(x = "MeanLog2", y = "SDLog2",
                                   tooltip = "FeatureName")
           else
             mapping <- aes_string(x = "Mean", y = "SD",
                                   tooltip = "FeatureName")
           p <- ggpoint(mapping, ...)
         },
         "mean-sd-samples" = {
           if (log2scale)
             mapping <- aes_string(x = "MeanLog2", y = "SDLog2",
                                   tooltip = "SampleName")
           else
             mapping <- aes_string(x = "Mean", y = "SD",
                                   tooltip = "SampleName")
           p <- ggpoint(mapping, ...)
         },
         "positiveControl" = {
           posCtrl <- positiveControlSubset(object)
           posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
           x <- log2(featureData(posCtrl)[["ControlConc"]])
           extradata <-
             data.frame(RSquared =
                          assayDataApply(posCtrl, 2L,
                                         function(y) cor(x, log2t(y, 0.5)),
                                         elt = elt))
           extradata[["Low R-Squared"]] <- extradata[["RSquared"]] < 0.95
           extradata[["CustomTooltip"]] <-
             sprintf("%s\nR-Squared = %.4f", rownames(extradata),
                     extradata[["RSquared"]])

           mapping <-
             aes_string(x = "ControlConc", y = elt, group = "SampleName",
                        tooltip = "CustomTooltip")
           if (any(extradata[["Low R-Squared"]])) {
             mapping[["colour"]] <- as.name("Low R-Squared")
             for (i in c("line", "point")) {
               geom[[i]] <- unclass(geom[[i]])
               geom[[i]][["colour"]] <- NULL
               oldClass(geom[[i]]) <- "uneval"
             }
           }
           p <- ggline(posCtrl, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2")
         })
  p
}


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
