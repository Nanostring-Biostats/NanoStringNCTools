\name{log2t}
\alias{log2t}

\title{Logarithm Base 2 With Thresholding}

\description{
  Safe log2 calculation where values within [0, thresh] are thresholded to
  thresh prior to the transformation.
}

\usage{
log2t(x, thresh = 0.5)
}

\arguments{
  \item{x}{a numeric or complex vector.}
  \item{thresh}{a positive number specifying the threshold.}
}

\value{
  A vector of the same length as \code{x} containing the transformed values.
}

\details{
  For non-negative elements in \code{x}, calculates \code{log2(pmax(x, thresh))}.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{log2}}}

\examples{
log2t(0:8)
identical(log2t(0:8), log2(c(0.5, 1:8)))
}

\keyword{math}
