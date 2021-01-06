\name{readRccFile}
\alias{readRccFile}

\title{Read RCC File}

\description{
  Read a NanoString Reporter Code Count (RCC) file.
}

\usage{
readRccFile(file)
}

\arguments{
  \item{file}{A character string containing the path to the RCC file.}
}

\value{
  An list object with five elements:
  \item{\code{"Header"}}{a \code{data.frame} object containing the header
    information.}
  \item{\code{"Sample_Attributes"}}{a \code{data.frame} object containing the
    attributes of the sample.}
  \item{\code{"Lane_Attributes"}}{a \code{data.frame} object containing the
    attributes of the lane.}
  \item{\code{"Code_Summary"}}{a \code{data.frame} object containing the
    reporter code counts.}
  \item{\code{"Messages"}}{A character vector containing messages, if any.}
}

\author{Patrick Aboyoun}

\seealso{\code{\link{readNanoStringRccSet}}}

\examples{
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rccData <- lapply(rccs, readRccFile)
}

\keyword{file}
\keyword{manip}
