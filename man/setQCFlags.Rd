\name{setQCFlags}
\alias{setQCFlags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set flags for QC of the assayData in a NanoStringRccSet.
}
\description{
This function takes a list containing the quality control (QC) thresholds for data in a NanoStringRccSet and then returns a matrix of QC retults by sample to protocolData.
}
\usage{
setQCFlags(object , qcCutoffs , hkGenes, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A valid NanoStringRccSet object with all housekeeping genes, positive control probes, and negative control probes present
}
  \item{qcCutoffs}{
A list with members named \code{Housekeeper}, \code{Imaging}, \code{BindingDensity}, \code{ERCCLinearity}, and \code{ERCCLoD}
}
  \item{hkGenes}{
A vector of housekeeping gene names if alternative genes to those defined in the panel are to be used
}
  \item{\dots}{
Any extra arguments to the function
}
}
\details{
This function checks that the housekeeping genes, positive control, and negative control probes or genes are within acceptable boundaries.  Borderline thresholds and fail thresholds are defined and each sample receives a row in a matrix that contains flags indicating either borderline or failing performance.

\code{Housekeeper} is a vector with names members.  \code{failingCutoff} sets the lower bound of housekeeper gene expression such that samples with a value below this threshold are labeled as failures.  \code{passingCutoff} sets a lower bound of housekeeper gene expression such that samples with a value below this threshold are labeled as borderline.  Values greater than or equal to either threshold are labeled as either borderline or passing.  The default values are \code{failingCutoff} = 32 and \code{passingCutoff} = 100.

\code{Imaging} is a vector with a single named member \code{fovCutoff}.  This threshold determines the mimimum proportion of FOV to be counted.  The default value is 0.75.

\code{BindingDensity} is a named vector with members \code{minimumBD}, \code{maximumBD}, and \code{maximumBDSprint}.  \code{minimumBD} sets a minimum threshold for binding density across machine platforms.  \code{maximumBD} sets a maxmimum binding density for non-Sprint machines while \code{maximumBDSprint} does the same for Sprint machines.  The default values are \code{minimumBD} = 0.1, \code{maximumBD} = 2.25, and \code{maximumBDSprint} = 1.8.

\code{ERCCLinearity} is a named vector with a single member \code{correlationValue}.  This member sets a minimum threshold for the correlation between the observed counts of positive controls and their theoretical concentration.  The default value is 0.95.

\code{ERCCLoD} is a named vector with a single member \code{standardDeviations}.  This sets a minimum threshold for the 0.5uMol concentration to be above the geoMean of the negative controls in units of standard deviation of the negative controls.  The default value is 2.
}
\value{
This function returns a new NanoStringRccSet with a matrix of QC pass criteria added to the \code{protocolData} slot called \code{QCFlags}.
}
\author{
Nicole Ortogero, Zhi Yang, and Patrick Aboyoun
}
\examples{
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")

solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)

solidTumor <- setQCFlags( solidTumor )
head( protocolData( solidTumor )[["QCFlags"]] )
}
\keyword{ QC }
