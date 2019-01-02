\name{colThresh}
\alias{colThresh}
\alias{colCenter}
\alias{colScale}

\title{Convenience Functions for Column Sweep Operations}

\description{
  Convenience functions for matrix thresholding, centering, and scaling based
  upon column statistics.
}

\usage{
colThresh(x, STATS)
colCenter(x, STATS)
colScale(x, STATS)
}

\arguments{
  \item{x}{a numeric array.}
  \item{STATS}{the summary statistic for thresholding, centering, or scaling.}
}

\value{
  An array with the same shape as \code{x} that has been modified by
  thresholding, centering, or scaling.
}

\details{
 These functions are convenience wrappers for the following code:
  \describe{
    \item{\code{colThresh}:}{\code{sweep(x, 2L, STATS, FUN = "pmax")}}
    \item{\code{colCenter}:}{\code{sweep(x, 2L, STATS, FUN = "-")}}
    \item{\code{colScale}:}{\code{sweep(x, 2L, STATS, FUN = "/")}}
  }
}

\author{Patrick Aboyoun}

\seealso{\code{\link{sweep}}}

\examples{
# Find reasonable column minimums
thresh <- apply(stack.x, 2L, quantile, 0.05)

# Threshold column values
identical(colThresh(stack.x, thresh),
          sweep(stack.x, 2L, thresh, FUN = "pmax"))

# Substract column values
identical(colCenter(stack.x, thresh),
          sweep(stack.x, 2L, thresh))

# Scale to common mean
identical(colScale(stack.x, colMeans(stack.x) / mean(colMeans(stack.x))),
          sweep(stack.x, 2L, colMeans(stack.x) / mean(colMeans(stack.x)),
                FUN = "/"))
}

\keyword{array}
\keyword{iteration}
