logt <-
function(x, thresh = 0.5)
{
  if (min(x, na.rm = TRUE) < thresh) {
    x[!is.na(x) & x >= 0 & x < thresh] <- thresh
  }
  log(x)
}

log2t <-
function(x, thresh = 0.5)
{
  if (min(x, na.rm = TRUE) < thresh) {
    x[!is.na(x) & x >= 0 & x < thresh] <- thresh
  }
  log2(x)
}
