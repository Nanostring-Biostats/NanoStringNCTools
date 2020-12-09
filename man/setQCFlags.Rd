\name{setQCFlags}
\alias{setQCFlags}
\title{
Flag QC data
}
\description{
Flag NanoStringRccSet samples that do not meet QC threshold requirements.
}
\usage{
  \code{setQCFlags(object,
             qcCutoffs = list(
               Housekeeper = c("failingCutoff" = 32,"passingCutoff" = 100) ,
               Imaging = c("fovCutoff" = 0.75) ,
               BindingDensity = c("minimumBD" = 0.1, "maximumBD" = 2.25, "maximumBDSprint" = 1.8) ,
               ERCCLinearity = c("correlationValue" = 0.95) ,
               ERCCLoD = c("standardDeviations" = 2)
             ),
             hkGenes = NULL,
             ReferenceSampleColumn = NULL,
             \ldots)}
}
\arguments{
  \item{object}{A NanoStringRccSet object}
  \item{qcCutoffs}{A list of QC cutoffs}
  \item{hkGenes}{A vector of housekeeping gene names}
  \item{ReferenceSampleColumn}{A string indicating the \code{pData} column containing reference sample information}
}
\value{
  An instance of the \code{\linkS4class{NanoStringRccSet}} class with QC flags.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
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



#Update cutoffs
newQCCutoffs <- list(
               Housekeeper = c("failingCutoff" = 10,"passingCutoff" = 100) ,
               Imaging = c("fovCutoff" = 0.5) ,
               BindingDensity = c("minimumBD" = 0.1, "maximumBD" = 2.25, "maximumBDSprint" = 1.8) ,
               ERCCLinearity = c("correlationValue" = 0.96) ,
               ERCCLoD = c("standardDeviations" = 3)
             )

#Set QC flags with new cutoffs
solidTumorNewQC <- setQCFlags(solidTumor, qcCutoffs=newQCCutoffs)
}
\keyword{datasets}
