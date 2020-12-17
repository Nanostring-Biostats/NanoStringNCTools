\name{normalize}
\alias{normalize}
\title{
Normalize RCCSet
}
\description{
This package performs normalization on NanoStringRccSet data using one of three methods.
}
\usage{
normalize(object, \ldots)
}
\arguments{
  \item{object}{
    \code{object} NanoStringRccSet object
}
  \item{\dots}{
    \code{object} additional arguments to pass on to normalize function
}
}
\details{
Normalization is performed in one of three ways with data pulled from one slot of assayData and inserted into another.  It is possible to overwrite the original slot of assayData if the fromElt and toElt are set to the same slot.  \code{nSolver} normalization uses positive controls to scale and housekeepers to standardize the data and mimics the normalization performed by default in the nSolver software.  The \code{Housekeeping-Log2} normalization calculates the log2 sizeFactor of the housekeeping genes and then takes 2^ log2 expression data centered by the log transformed sizeFactor.  \code{PositiveControl-Log2Log2} regresses the log2 positive control probes greater than 0.5 concentration on their geometric mean and then uses the intercept and slope to predict normalized values from the log2 transformed expression values.  The predictions are then rescaled by 2^.
Additional parameters with NanoStringRccSet method include:
\code{type} normalization method to use.  Options are \code{nSolver}, \code{Housekeeping-Log2}, and \code{PositiveControl-Log2Log2}
\code{fromElt} assayData slot from which to pull raw data
\code{toElt} assayData slot to which normalized data will be inserted
}
\value{
The function returns a new NanoStringRccSet with either an additional assayData slot of normalized data, or overwrites the original assayData depending on whether fromElt and toElt are identical.
}
\references{
NanoString nSolver User Manual
https://www.nanostring.com/download_file/view/1168
}
\author{Patrick Aboyoun}
\examples{
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")

solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)

solidTumor <- normalize(solidTumor, "nSolver" , fromElt = "exprs", toElt = "exprs_norm")
head( assayDataElement( solidTumor , elt = "exprs_norm" ) )
}
\keyword{ normalize }
