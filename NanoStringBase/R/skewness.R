skewness <-
function(x, na.rm = FALSE)
{
  if (anyNA(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    else
      return(NA_real_)
  }

  n <- length(x)
  x <- x - mean(x)
  ((sqrt(n * (n - 1L)) / (n - 2L)) * (sum(x^3) / n)) / ((sum(x^2) / n)^1.5)
}

kurtosis <-
function(x, na.rm = FALSE)
{
  if (anyNA(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    else
      return(NA_real_)
  }
  
  n <- length(x)
  x <- x - mean(x)
  ((n + 1L) * (n - 1L) *
      ((sum(x^4) / n) / (sum(x^2) / n)^2 - (3 * (n - 1L)) / (n + 1L))) /
    ((n - 2L) * (n - 3L))
}
