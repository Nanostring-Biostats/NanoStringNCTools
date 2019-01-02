# Zero-corrected logarithm functions
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


# Convenience wrappers around sweep function
colThresh <- function(x, STATS) sweep(x, 2L, STATS, FUN = "pmax")
colCenter <- function(x, STATS) sweep(x, 2L, STATS, FUN = "-")
colScale  <- function(x, STATS) sweep(x, 2L, STATS, FUN = "/")
