\name{normalize}
\alias{normalize}
\title{
Normalize RCCSet
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\code{normalize(object, type, fromElt = "exprs", toElt = "exprs_norm", \ldots)}
}
\arguments{
  \item{object}{
    \code{object} NanoStringRCCSet object
}
  \item{type}{
    \code{object} normalization method to use
}
  \item{fromElt}{
    \code{object} assayData slot from which to pull raw data
}
  \item{toElt}{
    \code{object} assayData slot to which normalized data will be inserted
}
  \item{\dots}{
    \code{object} additional arguments to pass on to normalize function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Patrick Aboyoun}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")

solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)

solidTumor <- normalize(solidTumor, "nSolver" , fromElt = "exprs", toElt = "exprs_norm")
}
\keyword{ normalize }
