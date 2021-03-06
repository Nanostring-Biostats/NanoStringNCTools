setGeneric("signatures", signature = "object", function(object) standardGeneric("signatures"))
setMethod("signatures", "NanoStringRccSet", function(object) object@signatures)
setGeneric("signatures<-", signature = c("object", "value"), function(object, value) standardGeneric("signatures<-"))
setReplaceMethod("signatures", c("NanoStringRccSet", "SignatureSet"), function(object, 
    value) {
    object@signatures <- value
    object
})
.sigCalc <- function(X, sigWeights) {
    t(vapply(sigWeights, function(wts) {
        if ("(Intercept)" %in% names(wts)) {
            X <- cbind(`(Intercept)` = 1, X)
        }
        if (all(names(wts) %in% colnames(X))) {
            X <- X[, names(wts), drop = FALSE]
            (X %*% wts)[, 1L]
        }
        else {
            structure(rep.int(NA_real_, nrow(X)), names = rownames(X))
        }
    }, FUN.VALUE=rep(numeric(1), nrow(X))))
}
setGeneric("signatureScores", signature = "object", function(object, ...) standardGeneric("signatureScores"))
setMethod("signatureScores", "NanoStringRccSet", function(object, elt = "exprs") {
    if (length(signatures(object)) == 0L) {
        return(matrix(numeric(), nrow = 0L, ncol = ncol(object), dimnames = list(NULL, 
            colnames(object))))
    }
    exprs <- t(assayDataElement2(object, elt))
    colnames(exprs) <- featureData(object)[["GeneName"]]
    sigFuncList <- signatureFuncs(object)
    linWeights <- weights(signatures(object))[names(sigFuncList)[which(sigFuncList %in% 
        "default")]]
    nonLinFuncs <- sigFuncList[which(!(sigFuncList %in% "default"))]
    scores <- .sigCalc(exprs, linWeights)
    while (length(idx <- which(rowSums(is.na(scores)) > 0L))) {
        subscores <- .sigCalc(cbind(exprs, t(scores[-idx, , drop = FALSE])), weights(signatures(object))[idx])
        if (all(is.na(subscores))) 
            break
        else scores[idx, ] <- subscores
    }
    nonLinScores <- t(vapply(nonLinFuncs, function(x, elt) eval(parse(text = paste(x, "( object , fromElt = \"", 
        elt, "\" )", sep = ""))), FUN.VALUE=rep(numeric(1), length(sampleNames(object))), elt))
    if (ncol(nonLinScores) > 0) {
        scores <- rbind(scores, nonLinScores)
    }
    return(scores[names(weights(signatures(object))), ])
})
setGeneric("signatureGroups", signature = "object", function(object, ...) standardGeneric("signatureGroups"))
setMethod("signatureGroups", "NanoStringRccSet", function(object) {
    groups(object@signatures)
})
setGeneric("setSignatureGroups<-", signature = c("object", "value"), function(object, value) standardGeneric("setSignatureGroups<-"))
setReplaceMethod("setSignatureGroups", c("NanoStringRccSet", "factor"), function(object, 
    value) {
    groups(object@signatures) <- value
    return(object)
})
setReplaceMethod("setSignatureGroups", c("NanoStringRccSet", "character"), function(object, 
    value) {
    groups(object@signatures) <- value
    return(object)
})
setGeneric("signatureFuncs", signature = "object", function(object, ...) standardGeneric("signatureFuncs"))
setMethod("signatureFuncs", "NanoStringRccSet", function(object) {
    getSigFuncs(object@signatures)
})
setGeneric("setSignatureFuncs<-", signature = c("object", "value"), function(object, value) standardGeneric("setSignatureFuncs<-"))
setReplaceMethod("setSignatureFuncs", c("NanoStringRccSet", "character"), function(object, 
    value) {
    setSigFuncs(object@signatures) <- value
    return(object)
})
