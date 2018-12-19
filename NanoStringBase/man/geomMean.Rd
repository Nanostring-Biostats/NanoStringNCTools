\name{geomMean}
\alias{geomMean}

\title{Geometric Mean}

\description{
  The (trimmed) geometric mean calculated using \code{exp(mean(log(x), \dots))}.
}

\usage{
geomMean(x, trim = 0, na.rm = FALSE, \dots)
}

\arguments{
  \item{x}{An \R object.}
  \item{trim}{the fraction of observations within [0, 0.5] to be
    trimmed from each end of the sorted values of \code{x} before
    the geometric mean is computed.}
  \item{na.rm}{a logical value indicating whether \code{NA}
    values should be removed from the geometric mean calculation.}
  \item{\dots}{additional arguments for \code{?base::\link[base]{mean}}.}
}

\value{
  See \code{?base::\link[base]{mean}} for the value returned by the
  default method.
}

\seealso{
  \code{base::\link[base]{mean}}.
}

\examples{
set.seed(12345)
x <- rlnorm(100)
xm <- geomMean(x)
c(xm, mean(x, trim = 0.10))
}

\keyword{univar}
