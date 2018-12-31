\name{log2t}
\alias{logt}
\alias{log2t}

\title{Logarithm With Thresholding}

\description{
  Safe log and log2 calculations where values within [0, thresh) are
  thresholded to thresh prior to the transformation.
}

\usage{
logt(x, thresh = 0.5)
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
  For non-negative elements in \code{x}, calculates \code{log(pmax(x, thresh))}
  or \code{log2(pmax(x, thresh))}.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{log}}, \code{\link{log2}}}

\examples{
logt(0:8)
identical(logt(0:8), log(c(0.5, 1:8)))

log2t(0:8)
identical(log2t(0:8), log2(c(0.5, 1:8)))
}

\keyword{math}
