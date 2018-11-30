\name{readNanoStringRccSet}
\alias{readNanoStringRccSet}
\concept{NanoStringRccSet}

\title{Read 'NanoStringRccSet'}

\description{
  Create an instance of class \code{\linkS4class{NanoStringRccSet}} by reading
  data from NanoString Reporter Code Count (RCC) files.
}

\usage{
readNanoStringRccSet(rccFiles, rlfFile = NULL)
}

\arguments{
  \item{rccFiles}{A character vector containing the paths to the RCC files.}
  \item{rlfFile}{An optional character string representing the path to the
  corresponding RLF file.}
}

\value{
  An instance of the \code{\linkS4class{NanoStringRccSet}} class.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{NanoStringRccSet}}, \code{\link{writeNanoStringRccSet}}}

\examples{
}
\keyword{file}
\keyword{manip}
