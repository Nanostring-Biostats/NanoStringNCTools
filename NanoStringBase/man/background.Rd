\name{bgThresh}
\alias{bgThresh}
\alias{bgSubtr}
\alias{fgScale}

\title{Background and Foreground Correction Helper Functions}

\description{
  Helper functions for Background Thresholding and Subtraction,
  and Foreground Scaling.
}

\usage{
bgThresh(x, STATS)
bgSubtr(x, STATS)
fgScale(x, STATS)
}

\arguments{
  \item{x}{a numeric array.}
  \item{STATS}{the summary statistic for thresholding, subtraction, or scaling.}
}

\value{
  An array with the same shape as \code{x} that has been modified by
  thresholding, subtraction, or scaling.
}

\details{
 These functions are convenience wrappers for the following code:
  \describe{
    \item{\code{bgThresh}:}{\code{sweep(x, 2L, STATS, FUN = "pmax")}}
    \item{\code{bgSubtr}:}{\code{sweep(x, 2L, STATS, FUN = "-")}}
    \item{\code{fgScale}:}{\code{sweep(x, 2L, STATS / mean(STATS), FUN = "/")}}
  }
}

\author{Patrick Aboyoun}

\seealso{\code{\link{sweep}}}

\examples{
# Find reasonable column minimums
thresh <- apply(stack.x, 2L, quantile, 0.05)

# Threshold column values
identical(bgThresh(stack.x, thresh),
          sweep(stack.x, 2L, thresh, FUN = "pmax"))

# Substract column values
identical(bgSubtr(stack.x, thresh),
          sweep(stack.x, 2L, thresh))

# Scale to common mean
identical(fgScale(stack.x, colMeans(stack.x)),
          sweep(stack.x, 2L, colMeans(stack.x) / mean(colMeans(stack.x)),
                FUN = "/"))
}

\keyword{array}
\keyword{iteration}
