# signatures Accessor and Replacer
setGeneric("signatures", signature = "object",
           function(object) standardGeneric("signatures"))
setMethod("signatures", "NanoStringRccSet",
          function(object) object@signatures)

setGeneric("signatures<-", signature = c("object", "value"),
           function(object, value) standardGeneric("signatures<-"))
setReplaceMethod("signatures", c("NanoStringRccSet", "SignatureSet"),
                 function(object, value) {
                   object@signatures <- value
                   object
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
            if (length(signatures(object)) == 0L) {
              return(matrix(numeric(), nrow = 0L, ncol = ncol(object),
                            dimnames = list(NULL, colnames(object))))
            }
            exprs <- assayDataElement2(object, elt)
            rownames(exprs) <- featureData(object)[["GeneName"]]
            scores <- .sigCalc(exprs, weights(signatures(object)))
            while (length(idx <- which(rowSums(is.na(scores)) > 0L))) {
              subscores <- .sigCalc(rbind(exprs, scores[-idx, , drop = FALSE]),
                                    weights(signatures(object))[idx])
              if (all(is.na(subscores)))
                break
              else
                scores[idx, ] <- subscores
            }
            scores
          })
