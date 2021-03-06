setMethod("[", "NanoStringRccSet", function(x, i, j, ..., drop = FALSE) {
    x <- callNextMethod()
    weights <- weights(signatures(x))
    grps <- signatures(x)@groups
    fncs <- signatures(x)@func
    if (length(weights) > 0L) {
        valid <- c("(Intercept)", featureData(x)[["GeneName"]])
        keep <- unlist(lapply(weights, function(y) all(names(y) %in% valid)))
        if (!all(keep)) {
            signatures(x) <- new("SignatureSet", weights = weights[keep], groups = grps[keep], 
                func = fncs[keep])
        }
    }
    x
})
setMethod("subset", "NanoStringRccSet", function(x, subset, select, ...) {
    if (!missing(subset)) {
        x <- x[eval(substitute(subset), fData(x), parent.frame(2L)), ]
    }
    if (!missing(select)) {
        x <- x[, eval(substitute(select), sData(x), parent.frame(2L))]
    }
    x
})
setGeneric("endogenousSubset", signature = "x", function(x, subset, select) standardGeneric("endogenousSubset"))
setMethod("endogenousSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(featureData(x)[["CodeClass"]] == "Endogenous"), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("housekeepingSubset", signature = "x", function(x, subset, select) standardGeneric("housekeepingSubset"))
setMethod("housekeepingSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(featureData(x)[["CodeClass"]] == "Housekeeping"), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("negativeControlSubset", signature = "x", function(x, subset, select) standardGeneric("negativeControlSubset"))
setMethod("negativeControlSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(featureData(x)[["CodeClass"]] == "Negative"), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("positiveControlSubset", signature = "x", function(x, subset, select) standardGeneric("positiveControlSubset"))
setMethod("positiveControlSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(featureData(x)[["CodeClass"]] == "Positive"), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("controlSubset", signature = "x", function(x, subset, select) standardGeneric("controlSubset"))
setMethod("controlSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(featureData(x)[["IsControl"]]), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("nonControlSubset", signature = "x", function(x, subset, select) standardGeneric("nonControlSubset"))
setMethod("nonControlSubset", "NanoStringRccSet", function(x, subset, select) {
    call <- match.call()
    call$x <- x[which(!featureData(x)[["IsControl"]]), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
setGeneric("signatureSubset", signature = "x", function(x, subset, select) standardGeneric("signatureSubset"))
setMethod("signatureSubset", "NanoStringRccSet", function(x, subset, select) {
    genes <- unique(names(unlist(unname(weights(signatures(x))))))
    call <- match.call()
    call$x <- x[which(featureData(x)[["GeneName"]] %in% genes), ]
    call[[1L]] <- as.name("subset")
    eval(call, parent.frame())
})
