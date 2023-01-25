autoplot.NanoStringRccSet <- function(object, type = c("boxplot-feature", "boxplot-signature", 
    "bindingDensity-mean", "bindingDensity-sd", "ercc-linearity", "ercc-lod", "heatmap-genes", 
    "heatmap-signatures", "housekeep-geom", "lane-bindingDensity", "lane-fov", "mean-sd-features", 
    "mean-sd-samples"), log2scale = TRUE, elt = "exprs", index = 1L, geomParams = list(), 
    tooltipDigits = 4L, heatmapGroup = NULL, blacklist = NULL, tooltipID = NULL, qcCutoffs = list(Housekeeper = c(failingCutoff = 32, 
        passingCutoff = 100), Imaging = c(fovCutoff = 0.75), BindingDensity = c(minimumBD = 0.1, 
        maximumBD = 2.25, maximumBDSprint = 1.8), ERCCLinearity = c(correlationValue = 0.95), 
        ERCCLoD = c(standardDeviations = 2)), scalingFactor = 1L, show_rownames_gene_limit = 60L, 
    show_colnames_gene_limit = 36L, show_rownames_sig_limit = 60L, show_colnames_sig_limit = 36L, 
    subSet = NULL, ...) {
    if (is.null(tooltipID)) {
        tooltipID <- dimLabels(object)[2L]
    }
    if ((length(geomParams) > 0)) {
        for (i in seq_along(geomParams)) {
            for (j in seq_along(geomParams[[i]])) {
                if (is(geomParams[[i]][[j]], "name")) {
                    charColName <- as.character(geomParams[[i]][[j]])
                    if (substr(charColName, nchar(charColName), nchar(charColName)) == "_") {
                        newLabel = substr(charColName, 1, (nchar(charColName) - 1))
                        pData(object)[newLabel] <- pData(object)[charColName]
                        if ("palette" %in% names(geomParams)) {
                            levels(geomParams[["palette"]][["dataframe"]][["Variable"]])[levels(geomParams[["palette"]][["dataframe"]][["Variable"]]) == 
                                geomParams[[i]][[j]]] <- newLabel
                        }
                        geomParams[[i]][[j]] <- as.name(newLabel)
                    }
                }
            }
        }
    }
    geomParams <- as.list(geomParams)
    geomParams <- update_geom_params("base", geomParams)
    geomParams <- update_geom_params("point", geomParams, GeomInteractivePoint$default_aes[!names(GeomInteractivePoint$default_aes) %in% 
        c("hover_css", "selected_css", "tooltip_fill", "hover_nearest")])
    geomParams <- update_geom_params("line", geomParams, GeomInteractiveLine$default_aes[!names(GeomInteractiveLine$default_aes) %in% 
        c("hover_css", "selected_css", "tooltip_fill", "hover_nearest")])
    geomParams <- update_geom_params("boxplot", geomParams, GeomInteractiveBoxplot$default_aes[!names(GeomInteractiveBoxplot$default_aes) %in% 
        c("hover_css", "selected_css", "tooltip_fill", "hover_nearest", "outlier.data_id",
          "outlier.tooltip", "outlier.onclick", "outlier.hover_css", 
          "outlier.selected_css", "outlier.tooltip_fill", "outlier.hover_nearest")])
    fontFamily <- ifelse(is.null(theme_get()$text$family), "Arial", theme_get()$text$family)
    ggpoint <- function(object, mapping, ...) {
        for (arg in names(geomParams[["point"]])) {
            if (is.name(geomParams[["point"]][[arg]])) {
                mapping[[arg]] <- geomParams[["point"]][[arg]]
                geomParams[["point"]][[arg]] <- NULL
            }
        }
        ggplot(object, mapping, elt = elt, ...) + do.call(geom_point_interactive, geomParams[["point"]])
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
        ggplot(object, mapping, elt = elt, ...) + do.call(geom_line_interactive, geomParams[["line"]]) + 
            do.call(geom_point_interactive, geomParams[["point"]])
    }
    type <- match.arg(type)
    switch(type, `boxplot-feature` = , `boxplot-signature` = {
        if (type == "boxplot-feature") {
            scores <- assayDataElement2(object, elt)
            ytitle <- fData(object)[index, dimLabels(object)[1L]]
            data_label <- "Expression"
        } else {
            scores <- signatureScores(object, elt)
            ytitle <- rownames(scores)[index]
            data_label <- "Score"
            if (!is.null(blacklist)) {
                scores <- scores[!rownames(scores) %in% blacklist, , drop = FALSE]
            }
        }
        colnames(scores) <- sData(object)[[dimLabels(object)[2L]]]
        if (!is.null(subSet)) {
            scores <- scores[, subSet]
        }
        if (type == "boxplot-feature") {
            y <- scores[index, ]
        } else {
            if (ytitle %in% rownames(scores)) {
                y <- scores[ytitle, ]
            } else {
                return(NULL)
            }
        }
        if (!is.name(geomParams[["base"]][["x"]])) {
            x <- rep.int("", ncol(scores))
            xtitle <- ""
        } else {
            if (is.null(subSet)) {
                x <- as.character(eval(geomParams[["base"]][["x"]], sData(object)))
            } else {
                x <- as.character(eval(geomParams[["base"]][["x"]], sData(object)))[subSet]
            }
            xtitle <- as.character(geomParams[["base"]][["x"]])
        }
        if (!is.name(geomParams[["point"]][["colour"]])) {
            colour <- NULL
            colourtitle <- ""
        } else {
            if (is.null(subSet)) {
                colour <- eval(geomParams[["point"]][["colour"]], sData(object))
            } else {
                colour <- eval(geomParams[["point"]][["colour"]], sData(object))[subSet]
            }
            colourtitle <- as.character(geomParams[["point"]][["colour"]])
            geomParams[["point"]] <- unclass(geomParams[["point"]])
            geomParams[["point"]][["colour"]] <- NULL
            oldClass(geomParams[["point"]]) <- "uneval"
            if ("palette" %in% names(geomParams)) {
                plot_palDF <- geomParams[["palette"]][["dataframe"]]
                pal_ind <- as.character(plot_palDF[["Variable"]]) == colourtitle
                plot_pal <- unlist(plot_palDF[pal_ind, "MainColor"])
                names(plot_pal) <- plot_palDF[pal_ind, "Level"]
            } else {
                plot_pal <- ggthemes::tableau_color_pal(palette = "Tableau 20")(20)
            }
        }
        tooltip <- colnames(scores)
        if (is.name(geomParams[["base"]][["x"]])) {
            tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, xtitle, x)
        }
        tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, ytitle, signif(y, tooltipDigits))
        if (is.name(geomParams[["point"]][["colour"]])) {
            tooltip <- sprintf("%s | %s&nbsp;=&nbsp;%s", tooltip, colourtitle, colour)
        }
        x <- gsub(pattern = "\n", replacement = "", x)
        tooltip <- gsub(pattern = "\n", replacement = "", tooltip)
        df <- data.frame(x = x, score = y, tooltip = tooltip, stringsAsFactors = FALSE)
        df[["colour"]] <- colour
        errorbar_width <- scalingFactor
        geomParams[["boxplot"]][["size"]] <- rel(0.2) * scalingFactor
        p <- ggplot(df, aes_string(x = "x", y = "score")) + stat_boxplot(geom = "errorbar", 
            lwd = rel(0.3) * scalingFactor, width = errorbar_width, colour = geomParams[["boxplot"]][["colour"]]) + 
            do.call(geom_boxplot_interactive, c(list(aes_string(tooltip = "x")), geomParams[["boxplot"]], 
                outlier.shape = NA)) + scale_x_discrete(name = xtitle, labels = lapply(levels(factor(df[["x"]])), 
            strwrpr)) + scale_y_continuous(name = data_label, labels = function(x) {
            sprintf("%.1f", x)
        })
        if (is.null(colour)) {
            p <- p + do.call(geom_beeswarm_interactive, c(list(aes_string(tooltip = "tooltip")), 
                geomParams[["point"]], geomParams[["beeswarm"]]))
        } else {
            p <- p + do.call(geom_beeswarm_interactive, c(list(aes_string(tooltip = "tooltip", 
                colour = "factor(colour)")), geomParams[["point"]], geomParams[["beeswarm"]])) + 
                scale_colour_manual(values = plot_pal) + guides(colour = guide_legend(title = colourtitle, 
                ncol = 1L, title.position = "top"))
            if (!is.null(geomParams[["showLegend"]])) {
                if (geomParams[["showLegend"]][["legend"]] == "off") {
                    p <- p + theme(legend.position = "none")
                }
            } else {
                p <- p + theme(legend.position = "right")
            }
        }
        if (xtitle == "") {
            p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.title.x = element_blank())
        } else {
            p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        }
        p
    }, `bindingDensity-mean` = , `bindingDensity-sd` = {
        if (log2scale) {
            y <- if (type == "bindingDensity-mean") "MeanLog2" else "SDLog2"
        } else {
            y <- if (type == "bindingDensity-mean") "Mean" else "SD"
        }
        mapping <- aes_string(x = "BindingDensity", y = y, tooltip = "SampleName")
        p <- ggpoint(object, mapping, ...) + scale_x_continuous(name = "Binding Density")
    }, `ercc-linearity` = {
        posCtrl <- positiveControlSubset(object)
        posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
        x <- log2(featureData(posCtrl)[["ControlConc"]])
        extradata <- data.frame(RSquared = assayDataApply(posCtrl, 2L, function(y) cor(x, 
            log2t(y, 0.5))^2, elt = elt))
        extradata[["Passing Correlation Value"]] <- extradata[["RSquared"]] >= qcCutoffs[["ERCCLinearity"]][["correlationValue"]]
        extradata[["CustomTooltip"]] <- sprintf("%s | R-Squared = %.4f", sData(object)[[tooltipID]], 
            extradata[["RSquared"]])
        mapping <- aes_string(x = "ControlConc", y = elt, group = "SampleName", tooltip = "CustomTooltip")
        PSCol <- pscheck(object)
        RSCol <- rscheck(object)
        for (i in c("line", "point")) {
            geomParams[[i]] <- unclass(geomParams[[i]])
            mapping[["colour"]] <- as.name("Passing Correlation Value")
            geomParams[[i]][["colour"]] <- NULL
            if (i == "point") {
                if (!is.null(PSCol) && !is.null(RSCol)) {
                    PSLabels <- getpslabels(object, PSCol, RSCol)
                } else if (!is.null(PSCol) && is.null(RSCol)) {
                    PSLabels <- getpslabels(object, PSCol, RSColumn = NULL)
                } else {
                    PSLabels <- rep("Sample", nrow(extradata))
                }
                mapping[["shape"]] <- rep(PSLabels, each = length(featureData(posCtrl)[["ControlConc"]]))
                geomParams[[i]][["shape"]] <- NULL
                oldClass(geomParams[[i]]) <- "uneval"
            }
        }
        geomParams[["line"]][["size"]] <- 0.5 * scalingFactor
        geomParams[["point"]][["size"]] <- 3 * scalingFactor
        geomParams[["point"]][["stroke"]] <- 1 * scalingFactor
        p <- ggline(posCtrl, mapping, extradata = extradata, ...) + scale_x_continuous(name = "Concentration (fM)", 
            trans = "log2") + scale_y_continuous(trans = "log2") + scale_colour_manual(values = c("#7ab800", 
            "#E15759"), limits = c(TRUE, FALSE), drop = FALSE) + guides(colour = guide_legend(ncol = 1L, 
            title.position = "top", order = 1))
        p <- p + guides(shape = guide_legend(title = "Sample Type", ncol = 1L, title.position = "top", 
            order = 0, override.aes = list(color = c("#7ab800", "#7ab800")))) + theme(legend.position = "right") + 
            scale_shape_manual(values = c(2, 16), guide = "none", limits = c("Panel Standard", 
                "Sample"), drop = FALSE)
        p <- p + theme(axis.ticks.x = element_line(size = scalingFactor * 0.2), axis.ticks.length = unit(3 * 
            scalingFactor, "pt"), axis.text = element_text(size = scalingFactor * 7), axis.title = element_text(size = scalingFactor * 
            10, face = "bold"), legend.title = element_text(size = scalingFactor * 8, face = "bold"), 
            legend.key.size = unit(20 * scalingFactor, "pt"), legend.text = element_text(size = scalingFactor * 
                6, margin = margin(t = scalingFactor * 5, b = scalingFactor * 5, unit = "pt")), 
            panel.border = element_rect(fill = NA, color = "black", size = scalingFactor * 
                0.25))
    }, `ercc-lod` = {
        negCtrl <- munge(negativeControlSubset(object), mapping = aes_(exprs = as.name(elt)))
        negCtrl[["x"]] <- negCtrl[["SampleName"]]
        cutoff <- negCtrl[["exprs"]]
        cutoff <- tapply(cutoff, negCtrl[["SampleName"]], function(x) mean(x, na.rm = TRUE)) + 
            qcCutoffs[["ERCCLoD"]][["standardDeviations"]] * tapply(cutoff, negCtrl[["SampleName"]], 
                function(x) sd(x, na.rm = TRUE))
        posCtrl <- positiveControlSubset(object)
        posCtrl <- subset(posCtrl, featureData(posCtrl)[["ControlConc"]] == 0.5)
        posCtrl <- munge(posCtrl, mapping = aes_(exprs = as.name(elt)))
        posCtrl[["tooltip"]] <- sprintf("%s | POS_E(0.5)&nbsp;=&nbsp;%s", sData(object)[[tooltipID]], 
            signif(posCtrl[["exprs"]], tooltipDigits))
        posCtrl[["x"]] <- posCtrl[["SampleName"]]
        posCtrl[["Limit of Detection"]] <- "Passing"
        posCtrl[["Limit of Detection"]][posCtrl$exprs < cutoff[posCtrl[["SampleName"]]]] <- "Failed"
        mapping <- aes_string(x = "x", y = elt, tooltip = "tooltip")
        PSCol <- pscheck(object)
        RSCol <- rscheck(object)
        geomParams[["point"]] <- unclass(geomParams[["point"]])
        mapping[["colour"]] <- as.name("Limit of Detection")
        geomParams[["point"]][["colour"]] <- NULL
        if (!is.null(PSCol) && !is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSCol)
        } else if (!is.null(PSCol) && is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSColumn = NULL)
        } else {
            PSLabels <- rep("Sample", nrow(posCtrl))
        }
        mapping[["shape"]] <- PSLabels
        geomParams[["point"]][["shape"]] <- NULL
        geomParams[["point"]][["size"]] <- 3 * scalingFactor
        geomParams[["point"]][["stroke"]] <- 1 * scalingFactor
        geomParams[["boxplot"]][["size"]] <- 0.5 * scalingFactor
        oldClass(geomParams[["point"]]) <- "uneval"
        geomParams[["line"]][["size"]] <- 0.5 * scalingFactor
        if (dimLabels(object)[2L] != "SampleName") {
            negCtrl[["x"]] <- as.character(rep(sData(object)[[dimLabels(object)[2L]]], 
                each = nrow(negativeControlSubset(object))))
            posCtrl[["x"]] <- as.character(sData(object)[[dimLabels(object)[2L]]])
        }
        indThreshold <- data.frame(x = seq_along(order(posCtrl[["x"]])) - 0.5, xend = seq_along(order(posCtrl[["x"]])) + 
            0.5, y = cutoff[order(posCtrl[["x"]])], yend = cutoff[order(posCtrl[["x"]])], 
            passingThreshold = rep("Passing Threshold", length(posCtrl[["x"]])))
        p <- ggplot(negCtrl, aes_string(x = "x", y = "exprs")) + stat_boxplot(geom = "errorbar", 
            colour = geomParams[["boxplot"]][["colour"]], size = 0.5 * scalingFactor) + 
            do.call(geom_boxplot_interactive, c(geomParams[["boxplot"]], outlier.shape = NA)) + 
            geom_beeswarm_interactive(posCtrl, mapping = mapping, size = geomParams[["point"]][["size"]], 
                fill = geomParams[["point"]][["fill"]], alpha = geomParams[["point"]][["alpha"]], 
                stroke = geomParams[["point"]][["stroke"]], groupOnX = FALSE) + geom_segment(aes(x = x, 
            xend = xend, y = y, yend = y, linetype = passingThreshold), indThreshold, color = "red", 
            size = 0.5 * scalingFactor) + scale_y_continuous(name = "Counts", trans = "log2") + 
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank()) + scale_colour_manual(values = c("#7ab800", 
            "#E15759"), limits = c("Passing", "Failed"), drop = FALSE) + scale_linetype_manual("Passing Threshold", 
            values = c(`Passing Threshold` = 1), name = NULL) + guides(colour = guide_legend(ncol = 1L, 
            title.position = "top", order = 1))
        if (nrow(posCtrl) <= 60L) {
            p <- p + theme(text = element_text(family = fontFamily), axis.text.x.bottom = element_text(angle = 90, 
                hjust = 1, vjust = 0.5, size = min(8, 420/nrow(posCtrl)) * scalingFactor))
        } else {
            p <- p + theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank())
        }
        p <- p + guides(shape = guide_legend(title = "Sample Type", ncol = 1L, title.position = "top")) + 
            theme(legend.position = "right") + scale_shape_manual(values = c(2, 16), guide = "none", 
            limits = c("Panel Standard", "Sample"), drop = FALSE)
        p <- p + theme(axis.ticks.y = element_line(size = scalingFactor * 0.2), axis.ticks.length = unit(3 * 
            scalingFactor, "pt"), axis.text = element_text(size = scalingFactor * 7), axis.title = element_text(size = scalingFactor * 
            10, face = "bold"), legend.title = element_text(size = scalingFactor * 8, face = "bold"), 
            legend.key.size = unit(20 * scalingFactor, "pt"), legend.text = element_text(size = scalingFactor * 
                6, margin = margin(t = scalingFactor * 5, b = scalingFactor * 5, unit = "pt")), 
            panel.border = element_rect(fill = NA, color = "black", size = scalingFactor * 
                0.25))
    }, `heatmap-genes` = {
        scores <- t(munge(endogenousSubset(object), ~GeneMatrix, elt = elt)[["GeneMatrix"]])
        p <- protoheatmap(scores, log2scale = log2scale, group = heatmapGroup, object = object, 
            show_rownames_limit = show_rownames_gene_limit, show_colnames_limit = show_colnames_gene_limit, 
            ...)
    }, `heatmap-signatures` = {
        scores <- t(munge(endogenousSubset(object), ~SignatureMatrix, elt = elt)[["SignatureMatrix"]])
        p <- protoheatmap(scores, log2scale = log2scale, group = heatmapGroup, object = object, 
            blacklist = blacklist, show_rownames_limit = show_rownames_sig_limit, show_colnames_limit = show_colnames_sig_limit, 
            ...)
    }, `housekeep-geom` = {
        hkGenes <- housekeepingSubset(object)
        hkStats <- summary(hkGenes, 2L, elt = "exprs")
        hkSet <- as.data.frame(hkStats)
        hkSet[["tooltip"]] <- sprintf("%s | Geometric&nbsp;Mean&nbsp;=&nbsp;%s", sData(object)[[tooltipID]], 
            signif(hkSet[["GeomMean"]], tooltipDigits))
        hkSet[["x"]] <- sData(hkGenes)[[dimLabels(hkGenes)[2L]]]
        cutX <- length(hkSet[["x"]])
        failingCutoff <- qcCutoffs[["Housekeeper"]][["failingCutoff"]]
        passingCutoff <- qcCutoffs[["Housekeeper"]][["passingCutoff"]]
        qcBorderlineText <- sprintf("Borderline < %s", passingCutoff)
        qcFailedText <- sprintf("Failed < %s", failingCutoff)
        qcPassedText <- sprintf("Passing >= %s", passingCutoff)
        hkSet[["Quality"]] <- qcPassedText
        hkSet$Quality[hkSet$GeomMean < passingCutoff] <- qcBorderlineText
        hkSet$Quality[hkSet$GeomMean < failingCutoff] <- qcFailedText
        hkSet <- transform(hkSet, x = reorder(x, GeomMean))
        mapping <- aes_string(x = "x", y = "GeomMean", tooltip = "tooltip")
        PSCol <- pscheck(object)
        RSCol <- rscheck(object)
        geomParams[["point"]] <- unclass(geomParams[["point"]])
        mapping[["colour"]] <- as.name("Quality")
        geomParams[["point"]][["colour"]] <- NULL
        if (!is.null(PSCol) && !is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSCol)
        } else if (!is.null(PSCol) && is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSColumn = NULL)
        } else {
            PSLabels <- rep("Sample", nrow(hkSet))
        }
        mapping[["shape"]] <- PSLabels
        geomParams[["point"]][["shape"]] <- NULL
        geomParams[["point"]][["size"]] <- 3 * scalingFactor
        geomParams[["point"]][["stroke"]] <- 1 * scalingFactor
        oldClass(geomParams[["point"]]) <- "uneval"
        p <- ggpoint(hkSet, mapping, ...) + geom_hline(yintercept = failingCutoff, linetype = 2L, 
            colour = "darkgray", size = 0.5 * scalingFactor) + geom_hline(yintercept = passingCutoff, 
            linetype = 2L, colour = "darkgray", size = 0.5 * scalingFactor) + geom_text(aes(cutX, 
            h, label = label, hjust = "right", vjust = 1.25), data = data.frame(h = c(failingCutoff), 
            label = c(sprintf("Minimum Threshold = %s counts", failingCutoff)), stringsAsFactors = FALSE), 
            color = "#79706E", size = 3 * scalingFactor, family = fontFamily, inherit.aes = FALSE) + 
            geom_text(aes(cutX, h, label = label, hjust = "right", vjust = -0.25), data = data.frame(h = c(passingCutoff), 
                label = c(sprintf("Borderline Threshold = %s counts", passingCutoff)), 
                stringsAsFactors = FALSE), color = "#79706E", size = 3 * scalingFactor, 
                family = fontFamily, inherit.aes = FALSE) + guides(colour = guide_legend(title = "Housekeeper Quality", 
            ncol = 1L, title.position = "top", order = 1, title.hjust = 0)) + scale_x_discrete(name = "Sample") + 
            scale_y_continuous(name = "Geometric Mean") + theme(legend.position = "right") + 
            scale_colour_manual(values = c("#7ab800", "#BAB0AC", "#E15759"), limits = c(qcPassedText, 
                qcBorderlineText, qcFailedText), drop = FALSE)
        if (length(hkSet[["x"]]) <= 60L) {
            p <- p + theme(text = element_text(family = fontFamily), axis.text.x.bottom = element_text(angle = 90, 
                hjust = 1, vjust = 0.5, size = 6 * scalingFactor))
        } else {
            p <- p + theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank())
        }
        p <- p + guides(shape = guide_legend(title = "Sample Type", ncol = 1L, title.position = "top", 
            order = 0, override.aes = list(color = c("#7ab800", "#7ab800")))) + theme(legend.position = "right") + 
            scale_shape_manual(values = c(2, 16), guide = "none", limits = c("Panel Standard", 
                "Sample"), drop = FALSE)
        p <- p + theme(axis.ticks.x = element_line(size = scalingFactor * 0.2), axis.ticks.length = unit(3 * 
            scalingFactor, "pt"), axis.text = element_text(size = scalingFactor * 6), axis.title = element_text(size = scalingFactor * 
            10, face = "bold"), legend.title = element_text(size = scalingFactor * 6, face = "bold"), 
            legend.text = element_text(size = scalingFactor * 6, margin = margin(t = scalingFactor * 
                5, b = scalingFactor * 5, unit = "pt")), panel.border = element_rect(fill = NA, 
                color = "black", size = scalingFactor * 0.25))
    }, `lane-bindingDensity` = {
        instrument <- substr(protocolData(object)[["ScannerID"]], 5, 5)
        minBD <- qcCutoffs[["BindingDensity"]][["minimumBD"]]
        maxBD <- qcCutoffs[["BindingDensity"]][["maximumBD"]]
        maxBDSprint <- qcCutoffs[["BindingDensity"]][["maximumBDSprint"]]
        if (any(instrument %in% "P")) {
            SPRINT <- TRUE
            if (all(instrument) %in% "P") {
                maxBD <- maxBDSprint
            }
        } else {
            SPRINT <- FALSE
        }
        if (length(unique(instrument)) > 1) {
            warning(sprintf("More than one instrument type in RCC set.  Using SPRINT threshold of %s instead of %s.", 
                maxBDSprint, maxBD))
        }
        extradata <- data.frame(PassingBD = unlist(apply(data.frame(bd = protocolData(object)[["BindingDensity"]], 
            i = instrument, min = minBD), 1, function(x) {
            maxBD <- switch(x[2], A = qcCutoffs[["BindingDensity"]][["maximumBD"]], B = qcCutoffs[["BindingDensity"]][["maximumBD"]], 
                C = qcCutoffs[["BindingDensity"]][["maximumBD"]], D = qcCutoffs[["BindingDensity"]][["maximumBD"]],
                E = qcCutoffs[["BindingDensity"]][["maximumBD"]],
                G = qcCutoffs[["BindingDensity"]][["maximumBD"]], H = qcCutoffs[["BindingDensity"]][["maximumBD"]], 
                P = qcCutoffs[["BindingDensity"]][["maximumBDSprint"]], default = qcCutoffs[["BindingDensity"]][["maximumBD"]])
            return(x[1] >= x[3] & x[1] <= maxBD)
        })), row.names = sampleNames(object))
        extradata[["CustomTooltip"]] <- sData(object)[[tooltipID]]
        mapping <- aes_string(x = "LaneID", y = "BindingDensity", tooltip = "CustomTooltip")
        PSCol <- pscheck(object)
        RSCol <- rscheck(object)
        geomParams[["point"]] <- unclass(geomParams[["point"]])
        mapping[["colour"]] <- as.name("PassingBD")
        geomParams[["point"]][["colour"]] <- NULL
        if (!is.null(PSCol) && !is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSCol)
        } else if (!is.null(PSCol) && is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSColumn = NULL)
        } else {
            hkGenes <- housekeepingSubset(object)
            hkStats <- summary(hkGenes, 2L, elt = "exprs")
            hkSet <- as.data.frame(hkStats)
            PSLabels <- rep("Sample", nrow(hkSet))
        }
        mapping[["shape"]] <- PSLabels
        geomParams[["point"]][["shape"]] <- NULL
        geomParams[["point"]][["size"]] <- 3 * scalingFactor
        geomParams[["point"]][["stroke"]] <- 1 * scalingFactor
        oldClass(geomParams[["point"]]) <- "uneval"
        cutX = 11
        p <- ggpoint(object, mapping, extradata = extradata, ...) + scale_x_continuous(name = "Lane", 
            breaks = seq_len(12), limits = c(1L, 12L)) + scale_y_continuous(name = "Binding Density", 
            limits = c(0, NA_real_)) + geom_hline(yintercept = c(minBD, maxBD), linetype = 2L, 
            colour = "darkgray", size = scalingFactor * 0.5) + geom_text(aes(cutX, h, label = label, 
            hjust = 0.55, vjust = 1.25), data = data.frame(h = c(minBD, maxBD), label = c("Minimum Binding Density", 
            "Maximum Binding Density"), stringsAsFactors = FALSE), color = "#79706E", size = scalingFactor * 
            3, family = fontFamily, inherit.aes = FALSE) + guides(colour = guide_legend(title = "Passing Binding Density", 
            ncol = 1L, title.position = "top", order = 1)) + scale_colour_manual(values = c("#7ab800", 
            "#E15759"), limits = c(TRUE, FALSE), drop = FALSE)
        if (SPRINT) {
            p <- p + geom_hline(yintercept = maxBDSprint, linetype = 2L, colour = "darkgray") + 
                geom_text(aes(cutX, h, label = label, hjust = 0.55, vjust = 1.25), data = data.frame(h = maxBDSprint, 
                  label = "SPRINT Binding Density", stringsAsFactors = FALSE), color = "#79706E", 
                  size = scalingFactor * 3, family = fontFamily, inherit.aes = FALSE)
        }
        p <- p + guides(shape = guide_legend(title = "Sample Type", ncol = 1L, title.position = "top")) + 
            scale_shape_manual(values = c(2, 16), guide = "none", limits = c("Panel Standard", 
                "Sample"), drop = FALSE)
        p <- p + theme(axis.ticks.x = element_line(size = scalingFactor * 0.2), axis.ticks.length = unit(3 * 
            scalingFactor, "pt"), axis.text = element_text(size = scalingFactor * 7), axis.title = element_text(size = scalingFactor * 
            10, face = "bold"), legend.title = element_text(size = scalingFactor * 8, face = "bold"), 
            legend.text = element_text(size = scalingFactor * 6, margin = margin(t = scalingFactor * 
                5, b = scalingFactor * 5, unit = "pt")), panel.border = element_rect(fill = NA, 
                color = "black", size = scalingFactor * 0.25))
    }, `lane-fov` = {
        extradata <- pData(protocolData(object))
        extradata <- data.frame(FOVCounted = extradata[["FovCounted"]]/extradata[["FovCount"]], 
            row.names = rownames(extradata))
        extradata[["Imaging Quality"]] <- sprintf("Passing >= %s%%", qcCutoffs[["Imaging"]][["fovCutoff"]] * 
            100)
        extradata$"Imaging Quality"[extradata$FOVCounted < qcCutoffs[["Imaging"]][["fovCutoff"]]] <- sprintf("Failed < %s%%", 
            qcCutoffs[["Imaging"]][["fovCutoff"]] * 100)
        extradata[["CustomTooltip"]] <- sData(object)[[tooltipID]]
        mapping <- aes_string(x = "LaneID", y = "FOVCounted", tooltip = "CustomTooltip")
        PSCol <- pscheck(object)
        RSCol <- rscheck(object)
        geomParams[["point"]] <- unclass(geomParams[["point"]])
        mapping[["colour"]] <- as.name("Imaging Quality")
        geomParams[["point"]][["colour"]] <- NULL
        if (!is.null(PSCol) && !is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSCol)
        } else if (!is.null(PSCol) && is.null(RSCol)) {
            PSLabels <- getpslabels(object, PSCol, RSColumn = NULL)
        } else {
            PSLabels <- rep("Sample", nrow(extradata))
        }
        mapping[["shape"]] <- PSLabels
        geomParams[["point"]][["size"]] <- 3 * scalingFactor
        geomParams[["point"]][["stroke"]] <- 1 * scalingFactor
        geomParams[["point"]][["shape"]] <- NULL
        oldClass(geomParams[["point"]]) <- "uneval"
        cutX = 11
        p <- ggpoint(object, mapping, extradata = extradata, ...) + scale_x_continuous(name = "Lane", 
            breaks = seq_len(12), limits = c(1L, 12L)) + scale_y_continuous(name = "FOV Counted", 
            labels = format_percent, limits = c(0, 1)) + geom_text(aes(cutX, h, label = label, 
            hjust = 0.1, vjust = 1.25), data = data.frame(h = c(qcCutoffs[["Imaging"]][["fovCutoff"]]), 
            label = c(sprintf("%s%% Passing", qcCutoffs[["Imaging"]][["fovCutoff"]] * 100)), 
            stringsAsFactors = FALSE), color = "#79706E", size = 3 * scalingFactor, family = fontFamily, 
            inherit.aes = FALSE) + geom_hline(yintercept = qcCutoffs[["Imaging"]][["fovCutoff"]], 
            linetype = 2L, colour = "darkgray", size = 0.5 * scalingFactor) + guides(colour = guide_legend(title = "Imaging Quality", 
            ncol = 1L, title.position = "top", order = 1)) + scale_colour_manual(values = c("#7ab800", 
            "#E15759"), limits = c(sprintf("Passing >= %s%%", qcCutoffs[["Imaging"]][["fovCutoff"]] * 
            100), sprintf("Failed < %s%%", qcCutoffs[["Imaging"]][["fovCutoff"]] * 100)), 
            drop = FALSE)
        p <- p + guides(shape = guide_legend(title = "Sample Type", ncol = 1L, title.position = "top", 
            order = 0, override.aes = list(color = c("#7ab800", "#7ab800")))) + theme(legend.position = "right") + 
            scale_shape_manual(values = c(2, 16), guide = "none", limits = c("Panel Standard", 
                "Sample"), drop = FALSE)
        p <- p + theme(axis.ticks.x = element_line(size = scalingFactor * 0.2), axis.ticks.length = unit(scalingFactor * 
            3, "pt"), axis.text = element_text(size = scalingFactor * 9), axis.title = element_text(size = scalingFactor * 
            10, face = "bold"), legend.title = element_text(size = scalingFactor * 10, 
            face = "bold"), legend.text = element_text(size = scalingFactor * 6, margin = margin(t = scalingFactor * 
            5, b = scalingFactor * 5, unit = "pt")), panel.border = element_rect(fill = NA, 
            color = "black", size = scalingFactor * 0.25))
    }, `mean-sd-features` = {
        if (log2scale) mapping <- aes_string(x = "MeanLog2", y = "SDLog2", tooltip = dimLabels(object)[1L]) else mapping <- aes_string(x = "Mean", 
            y = "SD", tooltip = dimLabels(object)[1L])
        p <- ggpoint(object, mapping, ...)
    }, `mean-sd-samples` = {
        if (log2scale) mapping <- aes_string(x = "MeanLog2", y = "SDLog2", tooltip = dimLabels(object)[2L]) else mapping <- aes_string(x = "Mean", 
            y = "SD", tooltip = dimLabels(object)[2L])
        p <- ggpoint(object, mapping, ...)
    })
    p
}
protoheatmap <- function(scores, log2scale, group, object, labelsize = 9L, scaleCutoff = 3, 
    groupPalette = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", 
        "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "maroon", "#4069FF", "#00FFFF", "#FF0707", 
        "#FF69B4", "#A66293", "#33CDF4", "#FFD861", "#ABABAB", "#318026"), blacklist = NULL, 
    annotation_colors = NULL, show_rownames_limit = 60L, show_colnames_limit = 36L, ...) {
    if (anyNA(rownames(scores))) {
        scores <- scores[!is.na(rownames(scores)), , drop = FALSE]
    }
    if (!is.null(blacklist)) {
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
    scores <- pmin(pmax(t(scale(t(scores))), -scaleCutoff), scaleCutoff)
    if (is.null(group)) {
        annotation_col <- NA
        annotation_colors <- NA
    }
    else {
        for (i in seq_along(group)) {
            if (substr(group[i], nchar(group[i]), nchar(group[i])) == "_") {
                newLabel = substr(group[i], 1, (nchar(group[i]) - 1))
                pData(object)[newLabel] <- sData(object)[group[i]]
                annotation_colors[[newLabel]] <- annotation_colors[[group[i]]]
                group[i] <- newLabel
            }
        }
        annotation_col <- sData(object)[group]
        annotation_col[] <- lapply(annotation_col, function(x) {
            if (is.factor(x)) {
                levels(x) <- trimws(levels(x))
                x[x == ""] <- NA
                x <- x[drop = TRUE]
            }
            else {
                x <- trimws(x)
                x[x == ""] <- NA
            }
            x <- addNA(x, ifany = TRUE)
            levels(x)[is.na(levels(x))] <- "N/A"
            x
        })
        rownames(annotation_col) <- make.unique(colnames(scores), sep = "_")
        if (is.null(annotation_colors)) {
            annotation_colors <- cumsum(vapply(annotation_col, nlevels, FUN.VALUE=numeric(1)))
            annotation_colors <- Map(`:`, c(1L, head(annotation_colors, -1L) + 1L), annotation_colors)
            annotation_colors <- structure(lapply(annotation_colors, function(x) {
                x <- x%%length(groupPalette)
                x[x == 0L] <- length(groupPalette)
                groupPalette[x]
            }), names = colnames(annotation_col))
            for (j in seq_len(ncol(annotation_col))) {
                names(annotation_colors[[j]]) <- levels(annotation_col[[j]])
            }
        }
        if (!is.null(annotation_colors)) {
            for (i in seq_along(annotation_colors)) {
                if (annotation_colors[i] == "Subtype") {
                  subtype_colors <- c(TIS = "maroon", LumA = "#4069FF", LumB = "#00FFFF", 
                    Basal = "#FF0707", `HER2-E` = "#FF69B4")
                  annotation_colors[["Subtype"]] <- subtype_colors[which(names(subtype_colors) %in% 
                    names(annotation_colors[[i]]))]
                  break
                }
            }
        }
    }
    pheatmap(scores, color = colorRampPalette(c("darkblue", rev(brewer.pal(n = 7L, name = "RdYlBu")), 
        "darkred"))(100), annotation_col = annotation_col, annotation_colors = annotation_colors, 
        cluster_rows = (nrow(scores) > 2), cluster_cols = (ncol(scores) > 2), show_rownames = (nrow(scores) <= 
            show_rownames_limit), show_colnames = (ncol(scores) <= show_colnames_limit), 
        silent = TRUE, legend.position = "bottom", fontsize = labelsize, angle_col = 90, 
        cellheight = ifelse(nrow(scores) <= 60L, labelsize + 2, NA), cellwidth = ifelse(ncol(scores) <= 
            36L, labelsize + 2, NA), fontfamily = ifelse(is.null(theme_get()$text$family), 
            "HersheySans", theme_get()$text$family), fontfamily_col = ifelse(is.null(theme_get()$text$family), 
            "HersheySans", theme_get()$text$family))
}
pscheck <- function(currObj) {
    if (length(dimLabels(currObj)) > 2) {
        panelStdCol <- dimLabels(currObj)[3]
    }
    else {
        panelStdCol <- NULL
    }
    return(panelStdCol)
}
rscheck <- function(currObj) {
    if (length(dimLabels(currObj)) > 3) {
        refSampleCol <- dimLabels(currObj)[4]
    }
    else {
        refSampleCol <- NULL
    }
    return(refSampleCol)
}
getpslabels <- function(currObj, PSColumn, RSColumn) {
    panelStandardLabels <- pData(currObj)[, PSColumn]
    panelStandardLabels[panelStandardLabels == 0] <- "Sample"
    panelStandardLabels[panelStandardLabels == 1] <- "Panel Standard"
    if (!is.null(RSColumn)) {
        referenceSampleLabels <- pData(currObj)[, RSColumn]
        for (i in seq_len(length(referenceSampleLabels))) {
            if (referenceSampleLabels[i] == 1) {
                panelStandardLabels[i] <- "Reference Sample"
            }
        }
    }
    return(panelStandardLabels)
}
