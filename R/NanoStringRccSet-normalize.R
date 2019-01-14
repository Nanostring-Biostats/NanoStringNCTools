setMethod("normalize", "NanoStringRccSet",
function(object, type = "PositiveControl-Log2Log2",
         fromElt = "exprs", toElt = "exprs_norm", ...)
{
  type <- match.arg(type)

  switch(type,
         "PositiveControl-Log2Log2" = {
           # Get the coefficient estimates for log2-log2 regression
           posCtrl <- positiveControlSubset(object)
           y <- summary(posCtrl, 1L)[, "GeomMean"]
           coefs <- esApply(posCtrl, 2L,
                            function(x) coef(lm(log2t(y) ~ log2t(x))))

           # Regress expression values to the mean           
           mat <- assayDataElement2(object, fromElt)
           mat <- sweep(log2t(mat), 2L, coefs[2L,], FUN = "*")
           mat <- .safe.as.integer(2 ^ sweep(mat, 2L, coefs[1L,], FUN = "+"))
         })

  assayDataElement(object, toElt) <- mat
  preproc(object)[[toElt]] <- match.call()
  object
})
