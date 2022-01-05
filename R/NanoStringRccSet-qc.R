setGeneric("setQCFlags", signature = "object", function(object, ...) standardGeneric("setQCFlags"))
setMethod("setQCFlags", "NanoStringRccSet", function(object, qcCutoffs = list(Housekeeper = c(failingCutoff = 32, 
    passingCutoff = 100), Imaging = c(fovCutoff = 0.75), BindingDensity = c(minimumBD = 0.1, 
    maximumBD = 2.25, maximumBDSprint = 1.8), ERCCLinearity = c(correlationValue = 0.95), 
    ERCCLoD = c(standardDeviations = 2)), hkGenes = NULL, ReferenceSampleColumn = NULL, 
    ...) {
    fovPercentLB = qcCutoffs[["Imaging"]][["fovCutoff"]]
    bindDenRange = c(qcCutoffs[["BindingDensity"]][["minimumBD"]], qcCutoffs[["BindingDensity"]][["maximumBD"]])
    posCtrlRsqLB = qcCutoffs[["ERCCLinearity"]][["correlationValue"]]
    negCtrlSDUB = qcCutoffs[["ERCCLoD"]][["standardDeviations"]]
    minHKGeoMean = qcCutoffs[["Housekeeper"]][["failingCutoff"]]
    blHKGeoMean = qcCutoffs[["Housekeeper"]][["passingCutoff"]]
    maxBindDen = qcCutoffs[["BindingDensity"]][["maximumBD"]]
    maxBindDenSprint = qcCutoffs[["BindingDensity"]][["maximumBDSprint"]]
    stopifnot(isSinglePercent(fovPercentLB))
    stopifnot(is.numeric(bindDenRange) && length(bindDenRange) == 2L && !anyNA(bindDenRange) && 
        bindDenRange[1L] >= 0 && bindDenRange[1L] < bindDenRange[2L])
    stopifnot(isSinglePercent(posCtrlRsqLB))
    negCtrl <- negativeControlSubset(object)
    posCtrl <- positiveControlSubset(object)
    posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
    controlConc <- featureData(posCtrl)[["ControlConc"]]
    prData <- protocolData(object)
    x <- log2(controlConc)
    if (is.null(hkGenes)) {
        subHKGenes <- housekeepingSubset(object)
    }
    else {
        preHKGenes <- housekeepingSubset(object)
        subHKGenes <- subset(preHKGenes, featureData(preHKGenes)[["GeneName"]] %in% hkGenes)
    }
    hkStats <- summary(subHKGenes, 2L, elt = "exprs")
    if (!is.null(ReferenceSampleColumn)) {
        pam50HKGenes <- c("MRPL19", "SF3A1", "PUM1", "ACTB", "PSMC4", "RPLP0", "GUSB", 
            "TFRC")
        subHKGenes <- subHKGenes[which(featureData(subHKGenes)[["GeneName"]] %in% pam50HKGenes), 
            which(as.logical(pData(subHKGenes)[[ReferenceSampleColumn]]))]
        thkStats <- summary(subHKGenes, 2L, elt = "exprs")
        hkStats[rownames(thkStats), ] <- thkStats
    }
    Binding <- unlist(apply(data.frame(prData[["BindingDensity"]], substr(protocolData(object)[["ScannerID"]], 
        5, 5), bindDenRange[1L]), 1, function(x) {
        maxBD <- switch(x[2], A = maxBindDen, B = maxBindDen, C = maxBindDen, D = maxBindDen, E = maxBindDen, 
            D = maxBindDen,G = maxBindDen, H = maxBindDen, P = maxBindDenSprint, default = maxBindDen)
        return(x[1] < x[3] | x[1] > maxBD)
    }))
    negCtrld <- munge(negCtrl, mapping = aes_(exprs = as.name("exprs")))
    cutoff <- negCtrld[["exprs"]]
    cutoff <- tapply(cutoff, negCtrld[["SampleName"]], function(x) mean(x, na.rm = TRUE)) + 
        negCtrlSDUB * tapply(cutoff, negCtrld[["SampleName"]], function(x) sd(x, na.rm = TRUE))
    prData[["QCFlags"]] <- cbind(Imaging = prData[["FovCounted"]]/prData[["FovCount"]] < 
        fovPercentLB, Binding = Binding, Linearity = assayDataApply(posCtrl, 2L, function(y) {
        cxy <- cor(x, log2t(y, 0.5))^2 < posCtrlRsqLB
        cxy <- ifelse(is.na(cxy), TRUE, cxy)
        return(cxy)
    }), LoD = apply(exprs(posCtrl[controlConc == 0.5, ]), 2L, max) < cutoff, Housekeeping = hkStats[, 
        "GeomMean"] < minHKGeoMean)
    prData[["QCBorderlineFlags"]] <- cbind(Imaging = rep(FALSE, nrow(prData@data)), Binding = rep(FALSE, 
        nrow(prData@data)), Linearity = rep(FALSE, nrow(prData@data)), LoD = rep(FALSE, 
        nrow(prData@data)), Housekeeping = hkStats[, "GeomMean"] > minHKGeoMean & hkStats[, 
        "GeomMean"] < blHKGeoMean)
    QCResults <- apply(prData[["QCFlags"]], 1L, function(x) {
        y <- sum(x) == 0L
        y <- ifelse(is.na(y), TRUE, y)
        return(y)
    })
    if (sum(QCResults) < 5) {
        stop("Unable to run 360 Report: less than five samples passed QC")
    }
    protocolData(object) <- prData
    preproc(object)["QCFlags"] <- list(match.call(call = sys.call(sys.parent(1L))))
    object
})
isSinglePercent <- function(x) {
    is.numeric(x) && length(x) == 1L && !anyNA(x) && x >= 0 && x <= 1
}
