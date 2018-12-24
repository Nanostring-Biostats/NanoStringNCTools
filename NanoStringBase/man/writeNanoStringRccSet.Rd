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
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringBase")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
solidTumorNoRlfPheno <- readNanoStringRccSet(rccs)
writeNanoStringRccSet(solidTumorNoRlfPheno, tempdir())
for (i in seq_along(rccs)) {
  stopifnot(identical(readLines(rccs[i]),
                      readLines(file.path(tempdir(), basename(rccs[i])))))
}
}

\keyword{file}
\keyword{manip}
