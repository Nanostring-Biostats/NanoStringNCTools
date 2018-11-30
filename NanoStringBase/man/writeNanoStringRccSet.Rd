\name{writeNanoStringRccSet}
\alias{writeNanoStringRccSet}
\concept{NanoStringRccSet}

\title{Write NanoString Reporter Code Count (RCC) files}

\description{
  Write NanoString Reporter Code Count (RCC) files from an instance of class
  \code{\linkS4class{NanoStringRccSet}}.
}

\usage{
writeNanoStringRccSet(x, dir = getwd())
}

\arguments{
  \item{x}{an instance of class \code{\linkS4class{NanoStringRccSet}.}}
  \item{dir}{An optional character string representing the path to the
  directory for the RCC files.}
}

\details{
  Writes a set of NanoString Reporter Code Count (RCC) files based upon \code{x}
  in \code{dir}.
}

\value{
  A character vector containing the paths for all the newly created RCC files.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{NanoStringRccSet}}, \code{\link{readNanoStringRccSet}}}

\examples{
}
\keyword{file}
\keyword{manip}
