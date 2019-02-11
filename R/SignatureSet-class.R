# Class definition
setClass("SignatureSet",
         contains = "VersionedBiobase",
         slots = c(weights = "NumericList"),
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(SignatureSet = "1.0.0")),
           weights = NumericList()))

# Initialization method
setMethod("initialize", "SignatureSet",
function(.Object, weights = NumericList(), ...)
{
  callNextMethod(.Object,
                 weights = as(weights, "NumericList"),
                 ...)
})

# Show method
setMethod("show", signature = "SignatureSet",
function(object) {
  callNextMethod(object)
  cat("weights: ")
  if (length(weights(object)) == 0L)
    cat("none\n")
  else
    cat("use 'weights(object)'")
})

# Constructors
SignatureSet <-
function(weights = NumericList(), ...)
{
  new2("SignatureSet", weights = weights, ...)
}

# Weights Accessor and Replacer
setMethod("weights", "SignatureSet",
          function(object, ...) object@weights)

setGeneric("weights<-", signature = c("object", "value"),
           function(object, value) standardGeneric("weights<-"))
setReplaceMethod("weights", c("SignatureSet", "NumericList"),
                 function(object, value) {
                   object@weights <- value
                   object
                 })
setReplaceMethod("weights", c("SignatureSet", "ANY"),
                 function(object, value) {
                   object@weights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("weights", c("SignatureSet", "list"),
                 function(object, value) {
                   value <- lapply(value, function(x) {
                     if(is.matrix(x) && ncol(x) == 1L)
                       structure(x[, 1L, drop = TRUE],
                                 names = rownames(x))
                     else
                       x
                   })
                   object@weights <- as(value, "NumericList")
                   object
                 })
setReplaceMethod("weights", c("SignatureSet", "NULL"),
                 function(object, value) {
                   object@weights <- NumericList()
                   object
                 })

# Additional methods
setMethod("length", "SignatureSet",
          function(x) {
            length(weights(x))
          })
setMethod("names", "SignatureSet",
          function(x) {
            names(weights(x))
          })
setMethod("lengths", "SignatureSet",
          function(x, use.names = TRUE) {
            sapply(weights(x),
                   function(x) sum(names(x) != "(Intercept)"),
                   USE.NAMES = use.names)
          })
