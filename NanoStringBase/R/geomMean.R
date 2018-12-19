geomMean <-
function(x, trim = 0, na.rm = FALSE, ...)
{
  exp(mean(log(x), trim = trim, na.rm = na.rm, ...))
}
