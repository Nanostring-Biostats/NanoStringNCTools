\name{readNanoStringRccSet}
\alias{readNanoStringRccSet}
\concept{NanoStringRccSet}

\title{Read 'NanoStringRccSet'}

\description{
  Create an instance of class \code{\linkS4class{NanoStringRccSet}} by reading
  data from NanoString Reporter Code Count (RCC) files.
}

\usage{
readNanoStringRccSet(rccFiles, rlfFile = NULL,
                     phenoDataFile = NULL,
                     phenoDataRccColName = "^RCC",
                     phenoDataColPrefix = "")
}

\arguments{
  \item{rccFiles}{A character vector containing the paths to the RCC files.}
  \item{rlfFile}{An optional character string representing the path to the
    corresponding RLF file.}
  \item{phenoDataFile}{An optional character string representing the path
    to the corresponding phenotypic csv data file.}
  \item{phenoDataRccColName}{The regular expression that specifies the RCC
    column in the \code{phenoDataFile}.}
  \item{phenoDataColPrefix}{An optional prefix to add to the phenoData
    column names to distinguish them from the names of assayData matrices,
    featureData columns, and protocolData columns.}
}

\value{
  An instance of the \code{\linkS4class{NanoStringRccSet}} class.
}

\author{Patrick Aboyoun}

\seealso{\code{\link{NanoStringRccSet}}, \code{\link{writeNanoStringRccSet}}}

\examples{
# Data file paths
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")

# Just RCC data
solidTumorNoRlfPheno <- readNanoStringRccSet(rccs)
varLabels(solidTumorNoRlfPheno)
fvarLabels(solidTumorNoRlfPheno)

# RCC and RLF data
solidTumorNoPheno <- readNanoStringRccSet(rccs, rlfFile = rlf)
setdiff(fvarLabels(solidTumorNoPheno), fvarLabels(solidTumorNoRlfPheno))

# All data
solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)
varLabels(solidTumor)
design(solidTumor) <- ~ BRAFGenotype + Treatment

# All data with phenoData prefix
solidTumorPhenoPrefix <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno,
                       phenoDataColPrefix = "PHENO_")
varLabels(solidTumorPhenoPrefix)
design(solidTumorPhenoPrefix) <- ~ PHENO_BRAFGenotype + PHENO_Treatment
}

\keyword{file}
\keyword{manip}
