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
.safe.as.integer <- function(x) {
  x <- round(x)
  storage.mode(x) <- "integer"
  x
}

fThresh <- function(x, STATS) sweep(x, 1L, STATS, FUN = "pmax")
fCenter <- function(x, STATS) sweep(x, 1L, STATS, FUN = "-")
fScale  <- function(x, STATS) sweep(x, 1L, STATS, FUN = "/")

sThresh <- function(x, STATS) sweep(x, 2L, STATS, FUN = "pmax")
sCenter <- function(x, STATS) sweep(x, 2L, STATS, FUN = "-")
sScale  <- function(x, STATS) sweep(x, 2L, STATS, FUN = "/")

fIntThresh <- function(x, STATS) .safe.as.integer(sweep(x, 1L, STATS, FUN = "pmax"))
fIntCenter <- function(x, STATS) .safe.as.integer(sweep(x, 1L, STATS, FUN = "-"))
fIntScale  <- function(x, STATS) .safe.as.integer(sweep(x, 1L, STATS, FUN = "/"))

sIntThresh <- function(x, STATS) .safe.as.integer(sweep(x, 2L, STATS, FUN = "pmax"))
sIntCenter <- function(x, STATS) .safe.as.integer(sweep(x, 2L, STATS, FUN = "-"))
sIntScale  <- function(x, STATS) .safe.as.integer(sweep(x, 2L, STATS, FUN = "/"))
