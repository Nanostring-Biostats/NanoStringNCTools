\name{readRlfFile}
\alias{readRlfFile}

\title{Read RLF File}

\description{
  Read a NanoString Reporter Library File (RLF) file.
}

\usage{
readRlfFile(file)
}

\arguments{
  \item{file}{A character string containing the path to the RLF file.}
}

\value{
  An instance of the \code{\linkS4class{DataFrame}} class containing columns:
  \item{"BarcodeClass"}{barcode class}
  \item{"GeneName"}{gene name}
  \item{"Accession"}{accession number}
  \item{\ldots}{additional columns}
}

\author{Patrick Aboyoun}

\seealso{\code{\link{readNanoStringRccSet}}}

\examples{
}
\keyword{file}
\keyword{manip}
