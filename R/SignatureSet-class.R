# Class definition
setClass("SignatureSet",
         contains = "VersionedBiobase",
         slots = c(weights = "NumericList", groups = "factor", version = "character"),
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(SignatureSet = "1.0.0")),
           weights = NumericList(),
           groups = factor(),
           version = character()))

# Initialization method
setMethod("initialize", "SignatureSet",
function(.Object, weights = NumericList(), ...)
{
  callNextMethod(.Object,
                 weights = as(weights, "NumericList"),
                 groups = factor( names( weights ) ),
                 version = "0.0.1",
                 ...)
})

# Initialization method
setMethod("initialize", "SignatureSet",
          function(.Object, weights = NumericList(), groups = factor(), ...)
          {
            callNextMethod(.Object,
                           weights = as(weights, "NumericList"),
                           groups = factor( groups ),
                           version = "0.0.1",
                           ...)
          })

# Initialization method
setMethod("initialize", "SignatureSet",
          function(.Object, weights = NumericList(), groups = factor(), version = character(), ...)
          {
            callNextMethod(.Object,
                           weights = as(weights, "NumericList"),
                           groups = factor( groups ),
                           version = as( version , "character" ),
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

SignatureSet <-
  function(weights = NumericList(), groups = factor(), ...)
  {
    new2("SignatureSet", weights = weights , groups = groups, ...)
  }

SignatureSet <-
  function(weights = NumericList(), groups = factor(), version = character(), ...)
  {
    new2("SignatureSet", weights = weights , groups = groups, version = version, ...)
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
setReplaceMethod("weights", c("SignatureSet", "CompressedNumericList"),
                 function(object, value) {
                   object@weights <- as(value, "NumericList")
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

# Groups Accessor and Replacer
setGeneric("groups", signature = c ("object" ) ,
           function( object , ... ) standardGeneric( "groups" ) )
setMethod( "groups" , "SignatureSet" ,
          function( object , ... ) object@groups )

setGeneric("groups<-", signature = c("object", "value"),
           function(object, value) standardGeneric("groups<-"))
setReplaceMethod("groups", c("SignatureSet", "factor"),
                 function(object, value) {
                   object@groups <- ifelse( inherits( value , "factor" ) , as( value , "factor" ) , factor( value ) )
                   object
                 })
setReplaceMethod("groups", c("SignatureSet", "ANY"),
                 function(object, value) {
                   object@groups <- ifelse( inherits( value , "factor" ) , as( value , "factor" ) , factor( value ) )
                   object
                 })
setReplaceMethod("groups", c("SignatureSet", "NULL"),
                 function(object, value) {
                   object@groups <- factor()
                   object
                 })

# Version Accessor and Replacer
setGeneric( "version" , signature = c( "object" ) ,
           function( object , ... ) standardGeneric( "version" ) )
setMethod( "version" , "SignatureSet" ,
          function( object , ... ) object@version )

setGeneric("version<-", signature = c("object", "value"),
           function(object, value) standardGeneric("version<-"))
setReplaceMethod("version", c("SignatureSet", "character"),
                 function(object, value) {
                   object@version <- value
                   object
                 })
setReplaceMethod("version", c("SignatureSet", "ANY"),
                 function(object, value) {
                   object@version <- as(value, "character")
                   object
                 })
setReplaceMethod("version", c("SignatureSet", "NULL"),
                 function(object, value) {
                   object@version <- character()
                   object
                 })

# Additional methods
# setMethod("version", "SignatureSet",
#           function( object ) {
#             return( object@version )
#           })
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

## Set validity functions to check for continuity between weigthts and groups
## A valid "track" object has the same number of x, y values
validSignatureSet <- function( object )
{
  errorMessage <- NULL
  if ( length( object@weights ) != length( object@groups ) )
  {
    errorMessage <- paste( "Unequal weights and groups.  Lengths: " , length( object@weights ) , ", " , length( object@groups ) , sep="" )
  }
  if ( length( object@weights ) != length( union( names( object@weights ) , names( object@groups ) ) ) )
  {
    errorMessage <- switch( is.null( errorMessage ),
                            paste( "Unequal weights and groups.  Lengths: " , length( object@weights ) , ", " , length( object@groups ) , sep="" ) ,
                            paste( errorMessage , "Weight names do not match group member names" , sep = "\n" )
    )
  }
  ifelse( is.null( errorMessage ) , TRUE , errorMessage )
}

## assign the function as the validity method for the class
setValidity( "SignatureSet" , validSignatureSet )
