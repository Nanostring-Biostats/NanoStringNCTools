setGeneric("normalize", signature = "object", function(object, ...) standardGeneric("normalize"))
setMethod("normalize", "NanoStringRccSet", function(object, type = c("nSolver", "PositiveControl-Log2Log2", 
    "Housekeeping-Log2"), fromElt = "exprs", toElt = "exprs_norm", ...) {
    type <- match.arg(type)
    switch(type, `Housekeeping-Log2` = {
        hkGenes <- housekeepingSubset(object)
        stats <- log2(summary(hkGenes, 2L, elt = fromElt)[, "SizeFactor"])
        mat <- log2t(assayDataElement2(object, fromElt))
        mat <- .safe.as.integer(2^sCenter(mat, stats))
    }, nSolver = {
        assayDataElement(object, toElt) <- 1 + assayDataElement(object, fromElt)
        posCtrl <- positiveControlSubset(object)
        pcG <- summary(posCtrl, 2L, elt = toElt)[, "GeomMean"]
        assayDataElement(object, toElt) <- sweep(assayDataElement(object, toElt), 2L, mean(pcG)/pcG, 
            FUN = "*")
        housekeepingSet <- housekeepingSubset(object)
        hkG <- summary(housekeepingSet, 2L, elt = toElt)[, "GeomMean"]
        mat <- log2(sweep(assayDataElement(object, toElt), 2L, mean(hkG)/hkG, FUN = "*"))
    }, `PositiveControl-Log2Log2` = {
        posCtrl <- positiveControlSubset(object)
        posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
        y <- summary(posCtrl, 1L, elt = fromElt)[, "GeomMean"]
        coefs <- assayDataApply(posCtrl, 2L, function(x) coef(lm(log2t(y) ~ log2t(x))))
        mat <- assayDataElement2(object, fromElt)
        mat <- sweep(log2t(mat), 2L, coefs[2L, ], FUN = "*")
        mat <- .safe.as.integer(2^sweep(mat, 2L, coefs[1L, ], FUN = "+"))
    })
    assayDataElement(object, toElt) <- mat
    preproc(object)[[toElt]] <- match.call()
    object
})
