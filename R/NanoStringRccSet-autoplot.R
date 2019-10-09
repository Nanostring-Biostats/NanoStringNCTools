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
         tooltipID = "SampleName",
         ...)
{
  geomParams <- as.list(geomParams)
  geomParams <- update_geom_params("base", geomParams)
  geomParams <- update_geom_params("point", geomParams, GeomInteractivePoint$default_aes)
  geomParams <- update_geom_params("line", geomParams, GeomInteractiveLine$default_aes)
  geomParams <- update_geom_params("boxplot", geomParams, GeomInteractiveBoxplot$default_aes)

  fontFamily <- ifelse( is.null( theme_get()$text$family ) , "Arial" , theme_get()$text$family )

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
               scores <- scores[!rownames(scores) %in% blacklist, , drop = FALSE]
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
             if("palette" %in% names(geomParams)) {
               plot_palDF <- geomParams[["palette"]][["dataframe"]]
               pal_ind <- as.character(plot_palDF[["Variable"]]) == colourtitle
               plot_pal <- unlist(plot_palDF[pal_ind, "MainColor"])
               names(plot_pal) <- plot_palDF[pal_ind, "Level"]
             } else {
               plot_pal <- tableau_color_pal(palette = "Tableau 20")(20)
             }
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
                          lwd = 0.5, # add this to boxplot geomParams later, sets error bar line width = box line width
                          width = geomParams[["boxplot"]][["size"]],
                          colour = geomParams[["boxplot"]][["colour"]]) +
             do.call(geom_boxplot_interactive,
                     c(list(aes_string(tooltip = "x")),
                       geomParams[["boxplot"]],
                       outlier.shape = NA)) +
             scale_x_discrete(name = xtitle) +
             scale_y_continuous(name = "Score",
                                labels = function(x) {sprintf("%.1f", x)})
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
               scale_colour_manual(values = plot_pal) +
               guides(colour = guide_legend(title = colourtitle,
                                            ncol = 1L,
                                            title.position = "top"))
             if ( !is.null( geomParams[["showLegend"]] ) )
             {
               if ( geomParams[["showLegend"]][["legend"]] == "off" )
               {
                 p <- p + theme( legend.position = "none" )
               }
             }
             else
             {
               p <- p + theme( legend.position = "right" )
             }
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
                                         function(y) cor(x, log2t(y, 0.5))^2,
                                         elt = elt))
           extradata[["Passing Correlation Value"]] <- extradata[["RSquared"]] >= 0.95
           extradata[["CustomTooltip"]] <-
             sprintf("%s | R-Squared = %.4f", object[[tooltipID]],
                     extradata[["RSquared"]])

           mapping <-
             aes_string(x = "ControlConc", y = elt, group = "SampleName",
                        tooltip = "CustomTooltip")
           
           # Check if panel standard exists
           PSCol <- pscheck(object)
           # Add custom aesthetics
           for (i in c("line", "point")) {
             # Separate for editing
             geomParams[[i]] <- unclass(geomParams[[i]])
             # Set color based on low R-squared
             mapping[["colour"]] <- as.name("Passing Correlation Value")
             # Remove default color
             geomParams[[i]][["colour"]] <- NULL
             if(i == "point") {
               if (!is.null(PSCol)) { 
                 # Get panel standard labels
                 PSLabels <- getpslabels(object, PSCol)
               } else {
                 PSLabels <- rep("Sample", nrows(extradata))
               }
               # Assign shape based on if panel standard
               mapping[["shape"]] <- 
               rep(PSLabels, each = length(featureData(posCtrl)[["ControlConc"]]))
               # Remove default shape
               geomParams[[i]][["shape"]] <- NULL
               # Reset class if setting new color or shape
               oldClass(geomParams[[i]]) <- "uneval"
             }
           }

           p <- ggline(posCtrl, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Concentration (fM)", trans = "log2") +
             scale_y_continuous(trans = "log2") + 
             scale_colour_manual(values = c("#7ab800", "#E15759"),
                                   limits = c(TRUE, FALSE),
                                   drop = FALSE) +
              guides(colour = guide_legend(ncol = 1L,
                                             title.position = "top", 
                                             order=1))
           # Add legend if panel standard provided
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                             ncol = 1L,
                                             title.position = "top",
                                             order = 0,
                                             override.aes = list(color=c("#7ab800", "#7ab800")))) + 
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), 
                                    guide = "none", 
                                    limits= c("Panel Standard", "Sample"),
                                    drop = FALSE)
         },
         "ercc-lod" = {
           negCtrl <- munge(negativeControlSubset(object),
                            mapping = aes_(exprs = as.name(elt)))
           negCtrl[["x"]] <- negCtrl[["SampleName"]]
           cutoff <- negCtrl[["exprs"]]
           cutoff <- tapply(cutoff, negCtrl[["SampleName"]] ,function( x ) mean( x , na.rm = TRUE ) ) +
               2 * tapply( cutoff , negCtrl[["SampleName"]] , function( x ) sd( x , na.rm = TRUE ) )

           posCtrl <- positiveControlSubset(object)
           posCtrl <- subset(posCtrl,
                             featureData(posCtrl)[["ControlConc"]] == 0.5)
           posCtrl <- munge(posCtrl, mapping = aes_(exprs = as.name(elt)))
           posCtrl[["tooltip"]] <-
             sprintf("%s | POS_E(0.5)&nbsp;=&nbsp;%s", object[[tooltipID]],
                     signif(posCtrl[["exprs"]], tooltipDigits))
           posCtrl[["x"]] <- posCtrl[["SampleName"]]
           posCtrl[["Limit of Detection"]] <- "Passing"
           posCtrl[["Limit of Detection"]][posCtrl$exprs < cutoff[posCtrl[["SampleName"]]]] <- "Failed"
           mapping <- aes_string(x = "x", y = elt, tooltip = "tooltip")

           # Check if panel standard exists
           PSCol <- pscheck(object)
           # Discriminate outliers and panel standards if designated
           # Separate for editing
           geomParams[["point"]] <- unclass(geomParams[["point"]])
           # Set color based on outlier
           mapping[["colour"]] <- as.name("Limit of Detection")
           # Remove default point color
           geomParams[["point"]][["colour"]] <- NULL
           if (!is.null(PSCol)) { 
             # Get panel standard labels
             PSLabels <- getpslabels(object, PSCol)
           } else {
             # Label all as samples
             PSLabels <- rep("Sample", nrows(posCtrl))
           }
           # Assign shape based on if panel standard
           mapping[["shape"]] <- PSLabels
           # Remove default shape
           geomParams[["point"]][["shape"]] <- NULL
           # Reset class if setting new color or shape
           oldClass(geomParams[["point"]]) <- "uneval"

           # Rename with sample labels
           if ( !( tooltipID %in% "SampleName" ) )
           {
             negCtrl[["x"]] <- rep( pData( object )[[tooltipID]] , each = nrow( negativeControlSubset(object) ) )
             posCtrl[["x"]] <- pData( object )[[tooltipID]]
           }

           # Create data for critical value threshold lines
           indThreshold <- data.frame( x = seq_along( order(posCtrl[["x"]]) ) - 0.5 ,
                                       xend = seq_along( order(posCtrl[["x"]]) ) + 0.5 ,
                                       y = cutoff[order(posCtrl[["x"]])] ,
                                       yend = cutoff[order(posCtrl[["x"]])] )
           # Set x position for cutoff line text
           p <- ggplot(negCtrl, aes_string(x = "x", y = "exprs")) +
             stat_boxplot(geom = "errorbar",
                          width = geomParams[["boxplot"]][["size"]],
                          colour = geomParams[["boxplot"]][["colour"]]) +
             do.call(geom_boxplot_interactive,
                     c(geomParams[["boxplot"]],
                       outlier.shape = NA)) +
             geom_beeswarm_interactive( posCtrl ,
                       mapping = mapping ,
                       size = geomParams[["point"]][["size"]] ,
                       fill = geomParams[["point"]][["fill"]] ,
                       alpha = geomParams[["point"]][["alpha"]] ,
                       stroke = geomParams[["point"]][["stroke"]] ,
                       groupOnX = FALSE ) +
             geom_segment( aes( x = x , xend = xend , y = y , yend = y ) , indThreshold , color = "red" ) +
             scale_y_continuous( name = "Counts" , trans = "log2" ) +
             theme( axis.ticks.x = element_blank() ,
                    axis.title.x = element_blank() ) +
             scale_colour_manual(values = c("#7ab800", "#E15759"),
                                 limits = c("Passing", "Failed"),
                                 drop = FALSE) +
             guides(colour = guide_legend(ncol = 1L,
                                          title.position = "top", 
                                          order=1))
           if ( nrow( posCtrl ) <= 60L )
           {
             p <- p + theme( text = element_text( family = fontFamily ) , 
                             axis.text.x.bottom = element_text( angle = 90 , hjust = 1 , vjust = 0.5 , size = min( 8 , 180 / nrow( posCtrl ) ) ) )
           }
           else
           {
             p <- p + theme( axis.text.x.bottom = element_blank() ,
                             axis.ticks.x = element_blank() )
           }
           # Add panel standard legend
           p <- p + 
             guides(shape = guide_legend(title = "Sample Type",
                                         ncol = 1L,
                                         title.position = "top")) +
             theme(legend.position = "right") +
             scale_shape_manual(values = c(2, 16), guide = "none", 
                                limits= c("Panel Standard", "Sample"),
                                drop = FALSE)
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
             sprintf("%s | Geometric&nbsp;Mean&nbsp;=&nbsp;%s", object[[tooltipID]],
                     signif(hkSet[["GeomMean"]], tooltipDigits))
           hkSet[["x"]] <- rownames(hkSet)
           # Get plot x limit for cut-off text
           cutX <- length(hkSet[["x"]])

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
           # Separate for editing
           geomParams[["point"]] <- unclass(geomParams[["point"]])
           # Set color based on quality
           mapping[["colour"]] <- as.name("Quality")
           # Remove default point color
           geomParams[["point"]][["colour"]] <- NULL
           # Set shapes based on panel standards
           if (!is.null(PSCol)) { 
               # Get sample and panel standard designations
               PSLabels <- getpslabels(object, PSCol)
           } else {
               # Set all to samples if no panel standard
               PSLabels <- rep("Sample", nrow(hkSet))
           }
           # Set shape based on label
           mapping[["shape"]] <- PSLabels
           # Remove default point shape
           geomParams[["point"]][["shape"]] <- NULL
           # Reset class if setting new color or shape
           oldClass(geomParams[["point"]]) <- "uneval"
           
           p <- ggpoint(hkSet, mapping, ...) +
             # Add lines indicating low quality or failing housekeepers
             geom_hline(yintercept = 32, linetype = 2L,
                          colour = "darkgray") +
             geom_hline(yintercept = 100, linetype = 2L,
                          colour = "darkgray") +
             geom_text(aes(cutX, h, label = label, hjust = "right", vjust = 1.25),
                         data =
                           data.frame(h = c(32), label = c("Minimum Threshold = 32 counts"),
                         stringsAsFactors = FALSE),
                         color = "#79706E", 
                         size = 3, 
                         family = fontFamily, 
                         inherit.aes = FALSE) +
             geom_text(aes(cutX, h, label = label, hjust = "right", vjust = -0.25),
                       data =
                         data.frame(h = c(100), label = c("Borderline Threshold = 100 counts"),
                                    stringsAsFactors = FALSE),
                       color = "#79706E", 
                       size = 3, 
                       family = fontFamily, 
                       inherit.aes = FALSE) +
             guides(colour = guide_legend(title = "Housekeeper Quality",
                                            ncol = 1L,
                                            title.position = "top",
                                            order = 1)) +
             scale_x_discrete(name="Sample") +
             scale_y_continuous(name="Geometric Mean") +
             theme(legend.position = "right") +
             scale_colour_manual(values = c("#7ab800", "#BAB0AC", "#E15759"),
                                   limits = c("Passing >= 100", 
                                                "Borderline < 100", 
                                                "Failed < 32"),
                                   drop = FALSE)
           if( length(hkSet[["x"]]) <= 60L ) {
             p <- p + theme(text = element_text(family=fontFamily), 
                            axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5))
           } else {
             p <- p + theme( axis.text.x.bottom = element_blank(),
                             axis.ticks.x = element_blank())
           }
          # Add shape legend
           p <- p + 
                  guides(shape = guide_legend(title = "Sample Type",
                           ncol = 1L,
                           title.position = "top",
                           order = 0,
                           override.aes = list(color=c("#7ab800", "#7ab800")))) +
                  theme(legend.position = "right") +
                  scale_shape_manual(values = c(2, 16), guide = "none", 
                                         limits= c("Panel Standard", "Sample"),
                                         drop = FALSE)
         },
         "lane-bindingDensity" = {
           instrument <- substr( protocolData( object )[["ScannerID"]] , 5 , 5 )
           minBD <- 0.1
           maxBD <- 2.25
           if ( any( instrument %in% "P" ) )
           {
             SPRINT <- TRUE
             if ( all( instrument ) %in% "P" )
             {
               maxBD <- 1.8
             }
           }
           else
           {
             SPRINT <- FALSE
           }
           if ( length( unique( instrument ) ) > 1 )
           {
             warning( "More than one instrument type in RCC set.  Using SPRINT threshold of 1.8 instead of 2.25.\n" )
           }
             extradata <-
               data.frame("PassingBD" = unlist( apply( data.frame( bd = protocolData(object)[["BindingDensity"]] , i = instrument , min = minBD ) , 1 ,
                              function( x )
                                {
                                  maxBD <- switch( x[2] , 
                                                   A = 2.25 ,
                                                   B = 2.25 ,
                                                   C = 2.25 ,
                                                   D = 2.25 ,
                                                   G = 2.25 ,
                                                   H = 2.25 ,
                                                   P = 1.8 ,
                                                   default = 2.25 )
                                  return( x[1] > x[3] & x[1] < maxBD )
                              } ) ) ,
                          row.names = sampleNames( object ) )
           extradata[["CustomTooltip"]] <- object[[tooltipID]]
           mapping <- aes_string(x = "LaneID", y = "BindingDensity",
                                 tooltip = "CustomTooltip")
           
           # Check if panel standard provided
           PSCol <- pscheck(object)
           # Discriminate outliers and panel standards if designated
           # Separate for editing
           geomParams[["point"]] <- unclass(geomParams[["point"]])
           # Set color based on outlier
           mapping[["colour"]] <- as.name("PassingBD")
           # Remove default point color
           geomParams[["point"]][["colour"]] <- NULL
           if (!is.null(PSCol)) {
             # Get panel standard labels
             PSLabels <- getpslabels(object, PSCol)
           } else {
             # Set all to samples if no panel standard
             PSLabels <- rep("Sample", nrow(hkSet))
           }
           # Assign shape based on if panel standard
           mapping[["shape"]] <- PSLabels
           # Remove default shape
           geomParams[["point"]][["shape"]] <- NULL
           # Reset class if setting new color or shape
           oldClass(geomParams[["point"]]) <- "uneval"
           
           # Set x position for cutoff line text
           cutX = 11
           p <- ggpoint(object, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "Binding Density",
                                limits = c(0, NA_real_)) +
             geom_hline(yintercept = c(0.1, maxBD), linetype = 2L,
                        colour = "darkgray") +
             geom_text(aes(cutX, h, label = label, hjust = 0.55, vjust = 1.25),
                       data =
                         data.frame(h = c(0.1, maxBD), label = c("Minimum Binding Density", "Maximum Binding Density"),
                                    stringsAsFactors = FALSE),
                       color = "#79706E", 
                       size = 3, 
                       family = fontFamily, 
                       inherit.aes = FALSE) + 
             guides(colour = guide_legend(title = "Passing Binding Density",
                                          ncol = 1L,
                                          title.position = "top",
                                          order = 1)) +
             scale_colour_manual(values = c("#7ab800", "#E15759"),
                                 limits = c(TRUE, FALSE),
                                 drop = FALSE)
           if ( SPRINT )
           {
             p <- p + geom_hline(yintercept = 1.8, linetype = 2L,
                        colour = "darkgray") +
               geom_text(aes(cutX, h, label = label, hjust = 0.55, vjust = 1.25),
                         data =
                           data.frame(h = 1.8, label = "SPRINT Binding Density",
                                      stringsAsFactors = FALSE),
                         color = "#79706E", 
                         size = 3, 
                         family = fontFamily, 
                         inherit.aes = FALSE)
           }
           # Add legend for panel standard
           p <- p + 
              guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top")) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none", 
                                  limits= c("Panel Standard", "Sample"),
                                  drop = FALSE)
         },
         "lane-fov" = {
           extradata <- pData(protocolData(object))
           extradata <-
             data.frame(FOVCounted =
                          extradata[["FovCounted"]] / extradata[["FovCount"]],
                        row.names = rownames(extradata))
           extradata[["Imaging Quality"]] <- "Passing >= 75%"
           extradata$"Imaging Quality"[extradata$"FOVCounted" < 0.75] <- "Failed < 75%"
           extradata[["CustomTooltip"]] <- object[[tooltipID]]
           mapping <- aes_string( x = "LaneID" , y = "FOVCounted" ,
                                  tooltip = "CustomTooltip" )
           # Check if panel standard provided
           PSCol <- pscheck(object)
           # Discriminate outliers and/or panel standards if designated
           geomParams[["point"]] <- unclass(geomParams[["point"]])
           # Set color based on quality
           mapping[["colour"]] <- as.name("Imaging Quality")
           # Remove default point color
           geomParams[["point"]][["colour"]] <- NULL
           if (!is.null(PSCol)) { 
             # Get panel standard labels
             PSLabels <- getpslabels(object, PSCol)
           } else {
             # Set all to samples if no panel standard
             PSLabels <- rep("Sample", nrow(extraData))
           }
           # Assign shape based on if panel standard
           mapping[["shape"]] <- PSLabels
           # Remove default shape
           geomParams[["point"]][["shape"]] <- NULL
           # Reset class if setting new color or shape
           oldClass(geomParams[["point"]]) <- "uneval"

           # Set x position for cutoff line text
           cutX = 11
           p <- ggpoint(object, mapping, extradata = extradata, ...) +
             scale_x_continuous(name = "Lane", breaks = 1:12,
                                limits = c(1L, 12L)) +
             scale_y_continuous(name = "FOV Counted", labels = format_percent,
                                limits = c(0, 1)) +
             geom_text(aes(cutX, h, label = label, hjust = 0.1, vjust = 1.25),
                       data =
                         data.frame(h = c(0.75), label = c("75% Passing"),
                                    stringsAsFactors = FALSE),
                       color = "#79706E", 
                       size = 3, 
                       family = fontFamily, 
                       inherit.aes = FALSE) +
             geom_hline(yintercept = 0.75, linetype = 2L,
                        colour = "darkgray") + 
             guides(colour = guide_legend(title = "Imaging Quality",
                                          ncol = 1L,
                                          title.position = "top",
                                          order = 1)) +
             scale_colour_manual(values = c("#7ab800", "#E15759"),
                                 limits = c("Passing >= 75%", "Failed < 75%"),
                                 drop = FALSE)
           # Add legend for panel standard
             p <- p + 
               guides(shape = guide_legend(title = "Sample Type",
                                           ncol = 1L,
                                           title.position = "top",
                                           order = 0,
                                           override.aes = list(color=c("#7ab800", "#7ab800")))) +
               theme(legend.position = "right") +
               scale_shape_manual(values = c(2, 16), guide = "none", 
                                    limits= c("Panel Standard", "Sample"),
                                    drop = FALSE)
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
           labelsize = 9L,
           scaleCutoff = 3,
           groupPalette =
             c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
               "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
           blacklist = NULL,
           annotation_colors = NULL,
           ...)
  {
  if (anyNA(rownames(scores))) {
    scores <- scores[!is.na(rownames(scores)), , drop = FALSE]
  }
  if ( !is.null( blacklist ) )
  {
    scores <- scores[!rownames(scores) %in% blacklist, , drop = FALSE]
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
    if(is.null(annotation_colors)) {
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
  }
  
  pheatmap(scores,
           color =
             colorRampPalette(c("darkblue",
                                rev(brewer.pal(n = 7L, name = "RdYlBu")),
                                "darkred"))(100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_rows = (nrow(scores) > 2),
           cluster_cols = (ncol(scores) > 2),
           show_rownames = (nrow(scores) <= 60L),
           show_colnames = (ncol(scores) <= 36L),
           silent = TRUE, legend.position = "bottom",
           fontsize = labelsize, 
           angle_col = 90, 
           cellheight = ifelse(nrow(scores) <= 60L, labelsize + 2, NA),
           cellwidth = ifelse(ncol(scores) <= 36L, labelsize + 2, NA),
           fontfamily = ifelse( is.null( theme_get()$text$family ) , "HersheySans" , theme_get()$text$family ) ,
           fontfamily_col = ifelse( is.null( theme_get()$text$family ) , "HersheySans" , theme_get()$text$family ) )
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
