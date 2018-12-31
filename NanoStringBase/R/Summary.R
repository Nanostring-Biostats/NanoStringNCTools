# Skewness
skewness <-
function(x, na.rm = FALSE)
{
  # Handle missing values
  if (anyNA(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    else
      return(NA_real_)
  }

  # Handle small sample sizes
  n <- length(x)
  if (n < 3L)
    return(NA_real_)

  # Calculate skewness
  x <- x - mean(x)
  ((sqrt(n * (n - 1L)) / (n - 2L)) * (sum(x^3) / n)) / ((sum(x^2) / n)^1.5)
}


# Excess kurtosis
kurtosis <-
function(x, na.rm = FALSE)
{
  # Handle missing values
  if (anyNA(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    else
      return(NA_real_)
  }

  # Handle small sample sizes
  n <- length(x)
  if (n < 4L)
    return(NA_real_)

  # Calculate excess kurtosis
  x <- x - mean(x)
  ((n + 1L) * (n - 1L) *
      ((sum(x^4) / n) / (sum(x^2) / n)^2 - (3 * (n - 1L)) / (n + 1L))) /
    ((n - 2L) * (n - 3L))
}


# Geometric mean
geomMean <-
function(x, level = 0.95, na.rm = FALSE)
{
  stats <-
    structure(rep.int(NA_real_, 3L),
              names =
                sprintf(c("GM_LCL%g", "GeomMean","GM_UCL%g"), 100 * level))

  # Handle missing values
  if (anyNA(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    else
      return(stats)
  }

  # Handle negative values
  if (any(x < 0)) {
    return(stats)
  }

  # Calculate geometric mean
  if (length(x) < 2L) {
    stats[2L] <- exp(mean(logt(x, thresh = 0.5)))
  } else {
    test <- t.test(logt(x, thresh = 0.5), conf.level = level)
    stats[2L] <- exp(test[["estimate"]])
    stats[c(1L, 3L)] <- exp(test[["conf.int"]])
  }

  stats
}
