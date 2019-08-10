setGeneric("setQCFlags", signature = "object",
           function(object, ...) standardGeneric("setQCFlags"))
setMethod("setQCFlags", "NanoStringRccSet",
function(object,
         fovPercentLB = 0.75,
         bindDenRange = c(0.1, 2.25),
         posCtrlRsqLB = 0.95,
         negCtrlSDUB = 2,
         hkGenes = NULL,
         minHKGeoMean = 32,
         blHKGeoMean = 100,
         ...)
{
  stopifnot(isSinglePercent(fovPercentLB))
  stopifnot(is.numeric(bindDenRange) &&
              length(bindDenRange) == 2L &&
              !anyNA(bindDenRange) &&
              bindDenRange[1L] >= 0 &&
              bindDenRange[1L] < bindDenRange[2L])
  stopifnot(isSinglePercent(posCtrlRsqLB))

  # Extract Negative and Positive Controls
  negCtrl <- negativeControlSubset(object)
  posCtrl <- positiveControlSubset(object)
  posCtrl <- posCtrl[featureData(posCtrl)[["ControlConc"]] >= 0.5, ]
  controlConc <- featureData(posCtrl)[["ControlConc"]]

  # Update protocolData with QC Flags
  prData <- protocolData(object)
  x <- log2(controlConc)

  # The generic case uses the genes as labeled
  if (is.null(hkGenes)) { 
    subHKGenes <- housekeepingSubset(object)
  } else {
    # Create NanoStringRCCSet for only the listed housekeepers
    preHKGenes <- housekeepingSubset(object)
    # Only use those housekeepers which data was collected for
    subHKGenes <- subset(preHKGenes, featureData( preHKGenes )[["GeneName"]] %in% hkGenes)
  }
  #Get geometric mean
  hkStats <- summary( subHKGenes , 2L , elt = "exprs" )

  prData[["QCFlags"]] <-
    cbind(Imaging =
            prData[["FovCounted"]] / prData[["FovCount"]] < fovPercentLB,
          Binding =
            prData[["BindingDensity"]] < bindDenRange[1L] |
            prData[["BindingDensity"]] > bindDenRange[2L],
          Linearity =
            assayDataApply(posCtrl, 2L,
                           function(y) cor(x, log2t(y, 0.5))^2 < posCtrlRsqLB),
          LoD =
            apply(exprs(posCtrl[controlConc == 0.5, ]), 2L, max) <=
            assayDataApply(negCtrl, 2L,
                           function(x) mean(x) + negCtrlSDUB * sd(x)),
          Housekeeping = 
            hkStats[, "GeomMean"] < minHKGeoMean)

  prData[["QCBorderlineFlags"]] <-
    cbind(Imaging = 
            rep(FALSE, nrow(prData@data)),
          Binding = 
            rep(FALSE, nrow(prData@data)),
          Linearity = 
            rep(FALSE, nrow(prData@data)),
          LoD = 
            rep(FALSE, nrow(prData@data)),
          Housekeeping =
            hkStats[, "GeomMean"] > minHKGeoMean &
            hkStats[, "GeomMean"] < blHKGeoMean)
  
  QCResults <- apply(prData[["QCFlags"]], 1L, function(x) sum(x) == 0L)
  if (sum(QCResults) < 5) {
    stop( "Unable to run 360 Report: less than five samples passed QC" )
  }

  protocolData(object) <- prData

  # Add method call to preproc list
  preproc(object)["QCFlags"] <-
    list(match.call(call = sys.call(sys.parent(1L))))

  object
})


isSinglePercent <- function(x) {
  is.numeric(x) && length(x) == 1L && !anyNA(x) && x >= 0 && x <= 1
}
