log2c0 <-
function(x)
{
  if (min(x, na.rm = TRUE) < 0.5) {
    x[!is.na(x) & x >= 0 & x < 0.5] <- 0.5
  }
  log2(x)
}
