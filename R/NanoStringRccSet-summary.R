skewness <- function(x, na.rm = FALSE) {
    if (anyNA(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        else return(NA_real_)
    }
    n <- length(x)
    if (n < 3L) 
        return(NA_real_)
    x <- x - mean(x)
    ((sqrt(n * (n - 1L))/(n - 2L)) * (sum(x^3)/n))/((sum(x^2)/n)^1.5)
}
kurtosis <- function(x, na.rm = FALSE) {
    if (anyNA(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        else return(NA_real_)
    }
    n <- length(x)
    if (n < 4L) 
        return(NA_real_)
    x <- x - mean(x)
    ((n + 1L) * (n - 1L) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - (3 * (n - 1L))/(n + 1L)))/((n - 
        2L) * (n - 3L))
}
geomMean <- function(x, na.rm = FALSE) {
    if (anyNA(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        else return(NA_real_)
    }
    if (min(x) < 0) {
        return(NA_real_)
    }
    exp(mean(logt(x, thresh = 0.5)))
}
.summaryMetadata <- list(log2 = data.frame(labelDescription = c("Geometric Mean", "Geometric Mean Size Factor", 
    "Mean Log2", "SD Log2"), row.names = c("GeomMean", "SizeFactor", "MeanLog2", "SDLog2"), 
    stringsAsFactors = FALSE), moments = data.frame(labelDescription = c("Mean", "Standard Deviation", 
    "Skewness", "Excess Kurtosis"), row.names = c("Mean", "SD", "Skewness", "Kurtosis"), 
    stringsAsFactors = FALSE), quantiles = data.frame(labelDescription = c("Minimum", "First Quartile", 
    "Median", "Third Quartile", "Maximum"), row.names = c("Min", "Q1", "Median", "Q3", 
    "Max"), stringsAsFactors = FALSE))
.marginalSummary <- function(x, log2scale = TRUE) {
    if (anyNA(x)) 
        x <- x[!is.na(x)]
    quartiles <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    names(quartiles) <- rownames(.summaryMetadata[["quantiles"]])
    if (log2scale) {
        log2X <- log2t(x, thresh = 0.5)
        stats <- structure(c(geomMean(x), NA_real_, mean(log2X), sd(log2X)), names = rownames(.summaryMetadata[["log2"]]))
    }
    else {
        stats <- structure(c(mean(x), sd(x), skewness(x), kurtosis(x)), names = rownames(.summaryMetadata[["moments"]]))
    }
    c(stats, quartiles)
}
setMethod("summary", "NanoStringRccSet", function(object, MARGIN = 2L, GROUP = NULL, log2scale = TRUE, 
    elt = "exprs", signatureScores = FALSE, ...) {
    stopifnot(MARGIN %in% c(1L, 2L))
    FUN <- function(x) {
        if (signatureScores) {
            applyFUN <- signatureScoresApply
        }
        else {
            applyFUN <- assayDataApply
        }
        stats <- t(applyFUN(x, MARGIN = MARGIN, FUN = .marginalSummary, log2scale = log2scale, 
            elt = elt))
        if (log2scale) {
            stats[, "SizeFactor"] <- 2^(stats[, "MeanLog2"] - mean(stats[, "MeanLog2"]))
        }
        stats
    }
    if (is.null(GROUP)) {
        FUN(object)
    }
    else {
        esBy(object, GROUP = GROUP, FUN = FUN, simplify = FALSE)
    }
})
