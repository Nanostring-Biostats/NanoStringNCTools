# signatureWeights Accessor and Replacer
setGeneric("signatureWeights", signature = "object",
           function(object) standardGeneric("signatureWeights"))
setMethod("signatureWeights", "NanoStringRccSet",
          function(object) object@signatureWeights)

setGeneric("signatureWeights<-", signature = c("object", "value"),
           function(object, value) standardGeneric("signatureWeights<-"))
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "NumericList"),
                 function(object, value) {
                   object@signatureWeights <- value
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "ANY"),
                 function(object, value) {
                   object@signatureWeights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "list"),
                 function(object, value) {
                   value <- lapply(value, function(x) {
                     if(is.matrix(x) && ncol(x) == 1L)
                       structure(x[, 1L, drop = TRUE],
                                 names = rownames(x))
                     else
                       x
                   })
                   object@signatureWeights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("signatureWeights", c("NanoStringRccSet", "NULL"),
                 function(object, value) {
                   object@signatureWeights <- NumericList()
                   object
                 })

# signatureNames Accessor
setGeneric("signatureNames", signature = "object",
           function(object) standardGeneric("signatureNames"))
setMethod("signatureNames", "NanoStringRccSet",
          function(object) {
            names(signatureWeights(object))
          })

# signatureLength Accessor
setGeneric("signatureLength", signature = "object",
           function(object) standardGeneric("signatureLength"))
setMethod("signatureLength", "NanoStringRccSet",
          function(object) {
            sapply(signatureWeights(object),
                   function(x) sum(names(x) != "(Intercept)"))
          })

# signatureScores Accessor and Replacer
.sigCalc <- function(X, sigWeights)
{
  t(sapply(sigWeights,
           function(wts) {
             if ("(Intercept)" %in% names(wts))
               X <- rbind("(Intercept)" = 1, X)
             if (all(names(wts) %in% rownames(X))) {
               X <- X[names(wts), , drop = FALSE]
               colSums(wts * X)
             } else {
               rep.int(NA_real_, ncol(X))
             }
           }))
}
setGeneric("signatureScores", signature = "object",
           function(object, ...) standardGeneric("signatureScores"))
setMethod("signatureScores", "NanoStringRccSet",
          function(object, elt = "exprs") {
            if (length(signatureWeights(object)) == 0L) {
              return(matrix(numeric(), nrow = 0L, ncol = ncol(object),
                            dimnames = list(NULL, colnames(object))))
            }
            exprs <- assayDataElement2(object, elt)
            rownames(exprs) <- featureData(object)[["GeneName"]]
            scores <- .sigCalc(exprs, signatureWeights(object))
            while (length(idx <- which(rowSums(is.na(scores)) > 0L))) {
              subscores <- .sigCalc(rbind(exprs, scores[-idx, , drop = FALSE]),
                                    signatureWeights(object)[idx])
              if (all(is.na(subscores)))
                break
              else
                scores[idx, ] <- subscores
            }
            scores
          })