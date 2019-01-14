\name{sThresh}

\alias{fThresh}
\alias{fCenter}
\alias{fScale}
\alias{sThresh}
\alias{sCenter}
\alias{sScale}

\alias{fIntThresh}
\alias{fIntCenter}
\alias{fIntScale}
\alias{sIntThresh}
\alias{sIntCenter}
\alias{sIntScale}

\alias{fAbove}
\alias{fBelow}
\alias{fAtLeast}
\alias{fAtMost}
\alias{sAbove}
\alias{sBelow}
\alias{sAtLeast}
\alias{sAtMost}

\title{Convenience Functions for Assay Data Element Sweep Operations}

\description{
  Convenience functions for matrix thresholding, centering, and scaling based
  upon margin statistics.
}

\usage{
# Loop over features
fThresh(x, STATS)
fCenter(x, STATS)
fScale(x, STATS)

## Round results to integers
fIntThresh(x, STATS)
fIntCenter(x, STATS)
fIntScale(x, STATS)

## Comparisons
fAbove(x, STATS)
fBelow(x, STATS)
fAtLeast(x, STATS)
fAtMost(x, STATS)


# Loop over samples
sThresh(x, STATS)
sCenter(x, STATS)
sScale(x, STATS)

# Round results to integers
sIntThresh(x, STATS)
sIntCenter(x, STATS)
sIntScale(x, STATS)

## Comparisons
sAbove(x, STATS)
sBelow(x, STATS)
sAtLeast(x, STATS)
sAtMost(x, STATS)
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
    \item{\code{fThresh}:}{\code{sweep(x, 1L, STATS, FUN = "pmax")}}
    \item{\code{fCenter}:}{\code{sweep(x, 1L, STATS, FUN = "-")}}
    \item{\code{fScale}:}{\code{sweep(x, 1L, STATS, FUN = "/")}}
    \item{\code{fIntThresh}:}{\code{round(sweep(x, 1L, STATS, FUN = "pmax"))}}
    \item{\code{fIntCenter}:}{\code{round(sweep(x, 1L, STATS, FUN = "-"))}}
    \item{\code{fIntScale}:}{\code{round(sweep(x, 1L, STATS, FUN = "/"))}}
    \item{\code{fAbove}:}{\code{sweep(x, 1L, STATS, FUN = ">")}}
    \item{\code{fBelow}:}{\code{sweep(x, 1L, STATS, FUN = "<")}}
    \item{\code{fAtLeast}:}{\code{sweep(x, 1L, STATS, FUN = ">=")}}
    \item{\code{fAtMost}:}{\code{sweep(x, 1L, STATS, FUN = "<=")}}
    \item{\code{sThresh}:}{\code{sweep(x, 2L, STATS, FUN = "pmax")}}
    \item{\code{sCenter}:}{\code{sweep(x, 2L, STATS, FUN = "-")}}
    \item{\code{sScale}:}{\code{sweep(x, 2L, STATS, FUN = "/")}}
    \item{\code{sIntThresh}:}{\code{round(sweep(x, 2L, STATS, FUN = "pmax"))}}
    \item{\code{sIntCenter}:}{\code{round(sweep(x, 2L, STATS, FUN = "-"))}}
    \item{\code{sIntScale}:}{\code{round(sweep(x, 2L, STATS, FUN = "/"))}}
    \item{\code{sAbove}:}{\code{sweep(x, 2L, STATS, FUN = ">")}}
    \item{\code{sBelow}:}{\code{sweep(x, 2L, STATS, FUN = "<")}}
    \item{\code{sAtLeast}:}{\code{sweep(x, 2L, STATS, FUN = ">=")}}
    \item{\code{sAtMost}:}{\code{sweep(x, 2L, STATS, FUN = "<=")}}
  }
}

\author{Patrick Aboyoun}

\seealso{\code{\link{sweep}}}

\examples{
# Find reasonable column minimums
thresh <- apply(stack.x, 2L, quantile, 0.05)

# Threshold column values
identical(sThresh(stack.x, thresh),
          sweep(stack.x, 2L, thresh, FUN = "pmax"))

# Substract column values
identical(sCenter(stack.x, thresh),
          sweep(stack.x, 2L, thresh))

# Scale to common mean
identical(sScale(stack.x, colMeans(stack.x) / mean(colMeans(stack.x))),
          sweep(stack.x, 2L, colMeans(stack.x) / mean(colMeans(stack.x)),
                FUN = "/"))

# Scale to common mean, rounded to the nearest integer
sIntScale(stack.x, colMeans(stack.x) / mean(colMeans(stack.x)))
}

\keyword{array}
\keyword{iteration}
