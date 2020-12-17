\name{setQCFlags}
\alias{setQCFlags}
\title{
Set flags for QC of the assayData in a NanoStringRccSet.
}
\description{
This function takes a list containing the quality control (QC) thresholds for data in a NanoStringRccSet and then returns a matrix of QC retults by sample to protocolData.
}
\usage{
       setQCFlags(object, \ldots)
}
\arguments{
  \item{object}{
  A valid NanoStringRccSet object with all housekeeping genes, positive control probes, and negative control probes present
  }

  \item{\dots}{
  Additional arguments to pass
  }
}
\details{
This function checks that the housekeeping genes, positive control, and negative control probes or genes are within acceptable boundaries.
Additional parameters with NanoStringRccSet method include:

\code{qcCutoffs} An optional list with members named \code{Housekeeper}, \code{Imaging}, \code{BindingDensity}, \code{ERCCLinearity}, and \code{ERCCLoD}

\code{hkGenes} An optional vector of housekeeping gene names if alternative genes to those defined in the panel are to be used

\code{ReferenceSampleColumn} An optional character string indicating the \code{pData} column containing reference sample information

Borderline thresholds and fail thresholds are defined and each sample receives a row in a matrix that contains flags indicating either borderline or failing performance.

\code{Housekeeper} is a vector with names members.  \code{failingCutoff} sets the lower bound of housekeeper gene expression such that samples with a value below this threshold are labeled as failures.  \code{passingCutoff} sets a lower bound of housekeeper gene expression such that samples with a value below this threshold are labeled as borderline.  Values greater than or equal to either threshold are labeled as either borderline or passing.  The default values are \code{failingCutoff} = 32 and \code{passingCutoff} = 100.

\code{Imaging} is a vector with a single named member \code{fovCutoff}.  This threshold determines the mimimum proportion of FOV to be counted.  The default value is 0.75.

\code{BindingDensity} is a named vector with members \code{minimumBD}, \code{maximumBD}, and \code{maximumBDSprint}.  \code{minimumBD} sets a minimum threshold for binding density across machine platforms.  \code{maximumBD} sets a maxmimum binding density for non-Sprint machines while \code{maximumBDSprint} does the same for Sprint machines.  The default values are \code{minimumBD} = 0.1, \code{maximumBD} = 2.25, and \code{maximumBDSprint} = 1.8.

\code{ERCCLinearity} is a named vector with a single member \code{correlationValue}.  This member sets a minimum threshold for the correlation between the observed counts of positive controls and their theoretical concentration.  The default value is 0.95.

\code{ERCCLoD} is a named vector with a single member \code{standardDeviations}.  This sets a minimum threshold for the 0.5uMol concentration to be above the geoMean of the negative controls in units of standard deviation of the negative controls.  The default value is 2.
}
\value{
This function returns a new \code{\linkS4class{NanoStringRccSet}} with matrices of QC pass and QC borderline criteria added to the \code{protocolData} slots called \code{QCFlags} and \code{QCBorderlineFlags}, respectively.
}
\examples{
# Create NanoStringRccSet from data files
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")
solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)

#Set QC flags with default cutoffs
solidTumorDefaultQC <- setQCFlags(solidTumor)
head( protocolData( solidTumorDefaultQC )[["QCFlags"]] )
head( protocolData( solidTumorDefaultQC )[["QCBorderlineFlags"]] )

#Update cutoffs
newQCCutoffs <- list(
  Housekeeper = c("failingCutoff" = 32,"passingCutoff" = 100) ,
  Imaging = c("fovCutoff" = 0.75) ,
  BindingDensity = c("minimumBD" = 0.1, "maximumBD" = 2.25, "maximumBDSprint" = 1.8) ,
  ERCCLinearity = c("correlationValue" = 0.98) ,
  ERCCLoD = c("standardDeviations" = 2)
  )

#Set QC flags with new cutoffs
solidTumorNewQC <- setQCFlags(solidTumor, qcCutoffs=newQCCutoffs)

#Compare QC results with default and new cutoffs
head( protocolData( solidTumorDefaultQC )[["QCFlags"]] )
head( protocolData( solidTumorNewQC )[["QCFlags"]] )
}
\keyword{datasets}
