autoplot.NanoStringRccSet <-
function(object,
         type = c("boxplot-feature",
                  "boxplot-signature",
                  "bindingDensity-mean",
                  "bindingDensity-sd",
                  "ercc-linearity",
                  "ercc-lod",
                  "heatmap-genes",
                  "heatmap-signatures",
                  "housekeep-geom",
                  "lane-bindingDensity",
                  "lane-fov",
                  "lod",
                  "mean-sd-features",
                  "mean-sd-samples"),
         log2scale = TRUE,
         elt = "exprs",
         index = 1L,
         geomParams = list(),
         tooltipDigits = 4L,
         heatmapGroup = NULL,
         blacklist = NULL,
         ...)
{
  geomParams <- as.list(geomParams)
  geomParams <- update_geom_params("base", geomParams)
  geomParams <- update_geom_params("point", geomParams, GeomInteractivePoint$default_aes)
  geomParams <- update_geom_params("line", geomParams, GeomInteractiveLine$default_aes)
  geomParams <- update_geom_params("boxplot", geomParams, GeomInteractiveBoxplot$default_aes)

  ggpoint <- function(object, mapping, ...) {
    for (arg in names(geomParams[["point"]])) {
      if (is.name(geomParams[["point"]][[arg]])) {
        mapping[[arg]] <- geomParams[["point"]][[arg]]
        geomParams[["point"]][[arg]] <- NULL
      }
    }
    ggplot(object, mapping, elt = elt, ...) +
      do.call(geom_point_interactive, geomParams[["point"]])
  }

  ggline <- function(object, mapping, ...) {
    for (i in c("line", "point")) {
      for (arg in names(geomParams[[i]])) {
        if (is.name(geomParams[[i]][[arg]])) {
          mapping[[arg]] <- geomParams[[i]][[arg]]
          geomParams[[i]][[arg]] <- NULL
        }
      }
    }
    ggplot(object, mapping, elt = elt, ...) +
      do.call(geom_line_interactive, geomParams[["line"]]) +
      do.call(geom_point_interactive, geomParams[["point"]])
  }

  type <- match.arg(type)
  switch(type,
         "boxplot-feature" =,
         "boxplot-signature" = {
           if (type == "boxplot-feature") {
             scores <- assayDataElement2(object, elt)
             ytitle <- fData(object)[index, dimLabels(object)[1L]]
           } else {
             scores <- signatureScores(object, elt)
             if ( !is.null( blacklist ) )
             {
               scores <- scores[-which( rownames( scores ) %in% blacklist ), , drop = FALSE]
             }
             ytitle <- rownames(scores)[index]
           }
           colnames(scores) <- sData(object)[[dimLabels(object)[2L]]]
           y <- scores[index, ]
           if (!is.name(geomParams[["base"]][["x"]])) {
             x <- rep.int("", ncol(scores))
             xtitle <- ""
           } else {
             x <- as.character(eval(geomParams[["base"]][["x"]], sData(object)))
             xtitle <- as.character(geomParams[["base"]][["x"]])
           }
           if (!is.name(geomParams[["point"]][["colour"]])) {
             colour <- NULL
             colourtitle <- ""
           } else {
             colour <- eval(geomParams[["point"]][["colour"]], sData(object))
             colourtitle <- as.character(geomParams[["point"]][["colour"]])
             geomParams[["point"]] <- unclass(geomParams[["point"]])
             geomParams[["point"]][["colour"]] <- NULL
             oldClass(geomParams[["point"]]) <- "uneval"
           }
           tooltip <- colnames(scores)
           if (is.name(geomParams[["base"]][["x"]])) {
             tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, xtitle, x)
           }
           tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, ytitle,
                              signif(y, tooltipDigits))
           if (is.name(geomParams[["point"]][["colour"]])) {
             tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, colourtitle, colour)
           }
           df <- data.frame(x = x, score = y, tooltip = tooltip,
                            stringsAsFactors = FALSE)
           df[["colour"]] <- colour
           p <- ggplot(df, aes_string(x = "x", y = "score")) +
             stat_boxplot(geom = "errorbar",
                          width = geomParams[["boxplot"]][["size"]],
                          colour = geomParams[["boxplot"]][["colour"]]) +
             do.call(geom_boxplot_interactive,
                     c(list(aes_string(tooltip = "x")),
                       geomParams[["boxplot"]],
                       outlier.shape = NA)) +
             scale_x_discrete(name = xtitle) + scale_y_continuous(name = ytitle)
           if (is.null(colour)) {
             p <- p +
               do.call(geom_beeswarm_interactive,
                       c(list(aes_string(tooltip = "tooltip")),
                         geomParams[["point"]],
                         geomParams[["beeswarm"]]))
           } else {
             p <- p +
               do.call(geom_beeswarm_interactive,
                       c(list(aes_string(tooltip = "tooltip",
                                         colour = "colour")),
                         geomParams[["point"]],
                         geomParams[["beeswarm"]])) +
               guides(colour = guide_legend(title = colourtitle,
                                              ncol = 1L,
                                              title.position = "top")) +
               theme(legend.position = "right")
           }
           if (xtitle == "") {
             p <- p + theme(axis.text.x  = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.title.x = element_blank())
           } else {
             p <- p +
               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
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
           p <- ggpoint(object, mapping, ...) +
             scale_x_continuous(name = "Binding Density")
         },
         "ercc-linearity" = {
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
             sprintf("%s | R-Squared = %.4f", rownames(extradata),
                     extradata[["RSquared"]])

           mapping <-
             aes_string(x = "ControlConc", y = elt, group = "SampleName",
                        tooltip = "CustomTooltip")
           
           # Check if panel standard exists
           PSCol <- pscheck(object)
           # Add custom aesthetics
           if (any(extradata[["Low R-Squared"]]) | 
                 !is.null(PSCol)) {
             for (i in c("line", "point")) {
             # Separate for editing
               geomParams[[i]] <- unclass(geomParams[[i]])
               if (any(extradata[["Low R-Squared"]]) & 
                     !is.name(geomParams[[i]][["colour"]])) {
                 # Set color based on low R-squared
                 mapping[["colour"]] <- as.name("Low R-Squared")
                 # Remove default color
                 geomParams[[i]][["colour"]] <- NULL
               }
               if(i == "point") {
                 if (!is.null(PSCol) & 
                       !is.name(geomParams[[i]][["shape"]])) { 
                   # Get panel standard labels
                   PSLabels <- getpslabels(object, PSCol)
                   # Assign shape based on if panel standard
                   mapping[["shape"]] <- 
                     rep(PSLabels, length(featureData(posCtrl)[["ControlConc"]]))
                   # Remove default shape
                   geomParams[[i]][["shape"]] <- NULL
                 }
               }
               # Reset class if setting new color or shape
               oldClass(geomParams[[i]]) <- "uneval"
             }
           }

           p <- ggline(posCtrl, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2")
           # Add legend if panel standard provided
           if (!is.null(PSCol)) {
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top")) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none")
           }
         },
         "ercc-lod" = {
           negCtrl <- munge(negativeControlSubset(object),
                            mapping = aes_(exprs = as.name(elt)))
           negCtrl[["x"]] <- ""
           if (log2scale) {
             cutoff <- log2t(negCtrl[["exprs"]])
             cutoff <- 2^(mean(cutoff, na.rm = TRUE) +
                            qnorm(0.975) * sd(cutoff, na.rm = TRUE))
           } else {
             cutoff <- negCtrl[["exprs"]]
             cutoff <- mean(cutoff, na.rm = TRUE) +
               qnorm(0.975) * sd(cutoff, na.rm = TRUE)
           }
           posCtrl <- positiveControlSubset(object)
           posCtrl <- subset(posCtrl,
                             featureData(posCtrl)[["ControlConc"]] == 0.5)
           posCtrl <- munge(posCtrl, mapping = aes_(exprs = as.name(elt)))
           posCtrl[["tooltip"]] <-
             sprintf("%s | POS_E(0.5)&nbsp;=&nbsp;%s", posCtrl[["SampleName"]],
                     signif(posCtrl[["exprs"]], tooltipDigits))
           posCtrl[["x"]] <- ""
           posCtrl[["Outlier"]] <- posCtrl[["exprs"]] < cutoff
           mapping <- aes_string(x = "x", y = elt, tooltip = "tooltip")

           # Check if panel standard exists
           PSCol <- pscheck(object)
           # Discriminate outliers and/or panel standards if designated
           if (any(posCtrl[["Outlier"]]) | 
                 !is.null(PSCol)) {
             # Separate for editing
             geomParams[["point"]] <- unclass(geomParams[["point"]])
             if (any(posCtrl[["Outlier"]]) & 
                   !is.name(geomParams[["point"]][["colour"]])) {
               # Set color based on outlier
               mapping[["colour"]] <- as.name("Outlier")
               # Remove default point color
               geomParams[["point"]][["colour"]] <- NULL
             }
             if (!is.null(PSCol) & 
                   !is.name(geomParams[["point"]][["shape"]])) { 
               # Get panel standard labels
               PSLabels <- getpslabels(object, PSCol)
               # Assign shape based on if panel standard
               mapping[["shape"]] <- PSLabels
               # Remove default shape
               geomParams[["point"]][["shape"]] <- NULL
             }
             # Reset class if setting new color or shape
             oldClass(geomParams[["point"]]) <- "uneval"
           }

           p <- ggplot(negCtrl, aes_string(x = "x", y = "exprs")) +
             stat_boxplot(geom = "errorbar",
                          width = geomParams[["boxplot"]][["size"]],
                          colour = geomParams[["boxplot"]][["colour"]]) +
             do.call(geom_boxplot_interactive,
                     c(geomParams[["boxplot"]],
                       outlier.shape = NA)) +
             do.call(geom_beeswarm_interactive,
                     c(list(posCtrl, mapping = mapping),
                       geomParams[["point"]],
                       geomParams[["beeswarm"]])) +
             geom_hline(yintercept = cutoff, linetype = 2L,
                        colour = "darkgray") +
             scale_x_discrete(name = "") +
             scale_y_continuous(name = "Counts (log2)", trans = "log2") +
             theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank())
           # Add legend if panel standard provided
           if (!is.null(PSCol)) {
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top")) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none")
           }
         },
         "heatmap-genes" = {
           scores <-
             t(munge(endogenousSubset(object), ~ GeneMatrix,
                     elt = elt)[["GeneMatrix"]])
           p <- protoheatmap(scores, log2scale = log2scale,
                             group = heatmapGroup, object = object, ...)
         },
         "heatmap-signatures" = {
           scores <-
             t(munge(endogenousSubset(object), ~ SignatureMatrix,
                     elt = elt)[["SignatureMatrix"]])
           p <- protoheatmap(scores, log2scale = log2scale,
                             group = heatmapGroup, object = object, blacklist = blacklist, ...)
         },
         "housekeep-geom" = { 
           # Extract housekeeping geometric mean data
           hkSet <- as.data.frame(object[["hkStats"]])
           hkSet[["tooltip"]] <- 
             sprintf("%s | Geometric&nbsp;Mean&nbsp;=&nbsp;%s", rownames(hkSet),
                     signif(hkSet[["GeomMean"]], tooltipDigits))
           hkSet[["x"]] <- rownames(hkSet)
           
           # Default to all values passing
           hkSet[["Quality"]] <- "Passing >= 100" 
           # Indicate which are borderline or failing quality
           hkSet$Quality[hkSet$GeomMean < 100] <- "Borderline < 100"
           hkSet$Quality[hkSet$GeomMean < 32] <- "Failed < 32"
           
           # Reorder by geometric mean value
           hkSet <- transform(hkSet, x=reorder(x, GeomMean) ) 
           mapping <- aes_string(x = "x", y = "GeomMean",
                                 tooltip = "tooltip")
           # function call to check if panel standard exists
           PSCol <- pscheck(object)

           # Discriminate lower quality values by color and panel standards by shape 
           if (any(hkSet[["Quality"]] != "Passing >= 100") | 
                 !is.null(PSCol)) {
             # Separate for editing
             geomParams[["point"]] <- unclass(geomParams[["point"]])
             if (any(hkSet[["Quality"]] != "Passing >= 100") & 
                   !is.name(geomParams[["point"]][["colour"]])) {
               # Set color based on quality
               mapping[["colour"]] <- as.name("Quality")
               # Remove default point color
               geomParams[["point"]][["colour"]] <- NULL
             }
             # Set shapes based on panel standards
             if (!is.null(PSCol) & 
                   !is.name(geomParams[["point"]][["shape"]])) { 
               # Get sample and panel standard designations
               PSLabels <- getpslabels(object, PSCol)
               # Set shape based on label
               mapping[["shape"]] <- PSLabels
               # Remove default point shape
               geomParams[["point"]][["shape"]] <- NULL
             }
             # Reset class if setting new color or shape
             oldClass(geomParams[["point"]]) <- "uneval"
           }
           
           p <- ggpoint(hkSet, mapping, ...) +
             # Add lines indicating low quality or failing housekeepers
             geom_hline(yintercept = 32, linetype = 2L,
                          colour = "darkgray") +
             geom_hline(yintercept = 100, linetype = 2L,
                          colour = "darkgray") +
             guides(colour = guide_legend(title = "Housekeeper Quality",
                                            ncol = 1L,
                                            title.position = "top")) +
             scale_x_discrete(name="Sample") +
             scale_y_continuous(name="Geometric Mean") +
             theme(legend.position = "right") +
             theme(text=element_text(family="TT Arial"), 
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
           
           # If there are panel standards add a shape legend
           if (!is.null(PSCol)) {
             p <- p + 
                    guides(shape = guide_legend(title = "Sample Type",
                               ncol = 1L,
                               title.position = "top")) +
                    theme(legend.position = "right") +
                    scale_shape_manual(values = c(2, 16), guide = "none")
           }
         },
         "lane-bindingDensity" = {
           maxBD <- 2.25
           extradata <-
             data.frame(Outlier =
                          protocolData(object)[["BindingDensity"]] > maxBD,
                        row.names = sampleNames(object))
           mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                 tooltip = "SampleName")
           
           # Check if panel standard provided
           PSCol <- pscheck(object)
           # Discriminate outliers and/or panel standards if designated
           if (any(extradata[["Outlier"]]) | 
                 !is.null(PSCol)) {
             # Separate for editing
             geomParams[["point"]] <- unclass(geomParams[["point"]])
             if (any(extradata[["Outlier"]]) & 
                   !is.name(geomParams[["point"]][["colour"]])) {
               # Set color based on outlier
               mapping[["colour"]] <- as.name("Outlier")
               # Remove default point color
               geomParams[["point"]][["colour"]] <- NULL
             }
             if (!is.null(PSCol) & 
                   !is.name(geomParams[["point"]][["shape"]])) { 
               # Get panel standard labels
               PSLabels <- getpslabels(object, PSCol)
               # Assign shape based on if panel standard
               mapping[["shape"]] <- PSLabels
               # Remove default shape
               geomParams[["point"]][["shape"]] <- NULL
             }
             # Reset class if setting new color or shape
             oldClass(geomParams[["point"]]) <- "uneval"
           }

           p <- ggpoint(object, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, maxBD), linetype = 2L,
                        colour = "darkgray")
           # Add legend if panel standard provided
           if (!is.null(PSCol)) {
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top")) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none")
           }
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
           # Check if panel standard provided
           PSCol <- pscheck(object)
           # Discriminate outliers and/or panel standards if designated
           if (any(extradata[["Outlier"]]) | 
                 !is.null(PSCol)) {
             # Separate for editing
             geomParams[["point"]] <- unclass(geomParams[["point"]])
             if (any(extradata[["Outlier"]]) & 
                   !is.name(geomParams[["point"]][["colour"]])) {
               # Set color based on outlier
               mapping[["colour"]] <- as.name("Outlier")
               # Remove default point color
               geomParams[["point"]][["colour"]] <- NULL
             }
             if (!is.null(PSCol) & 
                   !is.name(geomParams[["point"]][["shape"]])) { 
               # Get panel standard labels
               PSLabels <- getpslabels(object, PSCol)
               # Assign shape based on if panel standard
               mapping[["shape"]] <- PSLabels
               # Remove default shape
               geomParams[["point"]][["shape"]] <- NULL
             }
             # Reset class if setting new color or shape
             oldClass(geomParams[["point"]]) <- "uneval"
           }

           p <- ggpoint(object, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_hline(yintercept = 0.75, linetype = 2L,
                        colour = "darkgray")
           # Add legend if panel standard provided
           if (!is.null(PSCol)) {
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top")) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none")
           }
         },
         "mean-sd-features" = {
           if (log2scale)
             mapping <- aes_string(x = "MeanLog2", y = "SDLog2",
                                   tooltip = dimLabels(object)[1L])
           else
             mapping <- aes_string(x = "Mean", y = "SD",
                                   tooltip = dimLabels(object)[1L])
           p <- ggpoint(object, mapping, ...)
         },
         "mean-sd-samples" = {
           if (log2scale)
             mapping <- aes_string(x = "MeanLog2", y = "SDLog2",
                                   tooltip = dimLabels(object)[2L])
           else
             mapping <- aes_string(x = "Mean", y = "SD",
                                   tooltip = dimLabels(object)[2L])
           p <- ggpoint(object, mapping, ...)
         })
  p
}


protoheatmap <-
function(scores, log2scale, group, object,
         labelsize = 10L,
         scaleCutoff = 3,
         groupPalette =
           c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
             "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
         blacklist = NULL,
         ...)
{
  if (anyNA(rownames(scores))) {
    scores <- scores[!is.na(rownames(scores)), , drop = FALSE]
  }
  if ( !is.null( blacklist ) )
  {
    scores <- scores[-which( rownames( scores ) %in% blacklist ), , drop = FALSE]
  }
  if (anyNA(colnames(scores))) {
    ok <- which(!is.na(colnames(scores)))
    object <- object[, ok]
    scores <- scores[, ok, drop = FALSE]
  }

  if (log2scale) {
    scores <- log2t(scores)
  }

  scaleCutoff <- abs(scaleCutoff)
  scores <- pmin(pmax(t(scale(t(scores))), - scaleCutoff), scaleCutoff)

  if (is.null(group)) {
    annotation_col <- NA
    annotation_colors <- NA
  } else {
    annotation_col <- sData(object)[group]
    annotation_col[] <-
      lapply(annotation_col, function(x) {
        if (is.factor(x)) {
          levels(x) <- trimws(levels(x))
          x[x == ""] <- NA
          x <- x[drop = TRUE]
        } else {
          x <- trimws(x)
          x[x == ""] <- NA
        }
        x <- addNA(x, ifany = TRUE)
        levels(x)[is.na(levels(x))] <- "N/A"
        x
      })
    rownames(annotation_col) <- make.unique(colnames(scores), sep = "_")

    annotation_colors <- cumsum(sapply(annotation_col, nlevels))
    annotation_colors <- Map(`:`, c(1L, head(annotation_colors, -1L) + 1L),
                             annotation_colors)
    annotation_colors <- structure(lapply(annotation_colors, function(x) {
      x <- x %% length(groupPalette)
      x[x == 0L] <- length(groupPalette)
      groupPalette[x]
    }), names = colnames(annotation_col))
    for (j in seq_len(ncol(annotation_col))) {
      names(annotation_colors[[j]]) <- levels(annotation_col[[j]])
    }
  }

  pheatmap(scores,
           color =
             colorRampPalette(c("darkblue",
                                rev(brewer.pal(n = 7L, name = "RdYlBu")),
                                "darkred"))(100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = (nrow(scores) <= 60L),
           show_colnames = (ncol(scores) <= 36L),
           silent = TRUE, legend.position = "bottom",
           fontsize = labelsize, cellheight = labelsize + 2,
           fontfamily = "HersheySans")
}

# Check if panel standards were provided
pscheck <-
function(currObj) {
  # Check if panel standard column attached to rcc set dimLabels
  if(length(dimLabels(currObj)) > 2) {
    # Get panel standard column label
    panelStdCol <- dimLabels(currObj)[3]
  } else {
    # Panel standard not required for assay
    panelStdCol <- NULL
  }
  return(panelStdCol)
}

# Get panel standard and sample designations
getpslabels <-
function(currObj, PSColumn) {
  panelStandardLabels <- pData(currObj)[, PSColumn]
  # Label for legend
  panelStandardLabels[panelStandardLabels == 0] <- "Sample"
  panelStandardLabels[panelStandardLabels == 1] <- "Panel Standard"
  return(panelStandardLabels)
}
