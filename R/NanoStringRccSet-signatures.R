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
             if ("(Intercept)" %in% names(wts)) {
               X <- cbind("(Intercept)" = 1, X)
             }
             if (all(names(wts) %in% colnames(X))) {
               X <- X[, names(wts), drop = FALSE]
               (X %*% wts)[, 1L]
             } else {
               structure(rep.int(NA_real_, nrow(X)), names = rownames(X))
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
            exprs <- t(assayDataElement2(object, elt))
            colnames(exprs) <- featureData(object)[["GeneName"]]
            scores <- .sigCalc(exprs, weights(signatures(object)))
            while (length(idx <- which(rowSums(is.na(scores)) > 0L))) {
              subscores <-
                .sigCalc(cbind(exprs, t(scores[-idx, , drop = FALSE])),
                         weights(signatures(object))[idx])
              if (all(is.na(subscores)))
                break
              else
                scores[idx, ] <- subscores
            }
            scores
          })

setGeneric( "signatureGroups" , signature = "object" ,
           function (object , ... ) standardGeneric( "signatureGroups" ) )
setMethod("signatureGroups", "NanoStringRccSet",
          function( object ) {
            groups( object@signatures )
          } )
