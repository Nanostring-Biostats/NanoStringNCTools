\name{log2c0}
\alias{log2c0}

\title{Safe Logarithm Base 2}

\description{
  Safe log2 calculation where values within [0, 0.5] are thresholded to 0.5
  prior to the transformation.
}

\usage{
log2c0(x)
}

\arguments{
  \item{x}{a numeric or complex vector.}
}

\value{
  A vector of the same length as \code{x} containing the transformed values.
}

\details{
  For non-negative elements in \code{x}, calculates \code{log2(pmax(x, 0.5))}.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{log2}}}

\examples{
log2c0(0:8)
identical(log2c0(0:8), log2(c(0.5, 1:8)))
}

\keyword{math}
