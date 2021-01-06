\name{geom_beeswarm_interactive}
\alias{geom_beeswarm_interactive}

\title{Geometry for Interactive Bee Swarm Points}

\description{
  The interactive version of \code{\link[ggbeeswarm]{geom_beeswarm}} from
  \pkg{ggbeeswarm}.
}

\usage{
geom_beeswarm_interactive(mapping = NULL, data = NULL,
                          priority = c("ascending", "descending", "density",
                                       "random", "none"),
                          cex = 1, groupOnX = NULL, dodge.width = 0,
                          stat = "identity", na.rm = FALSE, show.legend = NA,
                          inherit.aes = TRUE, \ldots)
}

\arguments{
  \item{mapping}{The aesthetic mapping.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{data}{The data to be displayed at this layer.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{priority}{Method used to perform point layout.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{cex}{Scaling for adjusting point spacing.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{groupOnX}{Indicator for jittering on x-axis.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{dodge.width}{Dodge amount for points from different aesthetic groups.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{stat}{The statistical transformation to use on the data for this layer.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{na.rm}{Indicator for removing missing values with a warning.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{show.legend}{Indicator for including this layer in the legend.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{inherit.aes}{Indicator for inheriting the aesthetics.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
  \item{\ldots}{Additional arguments.
    See \code{\link[ggbeeswarm]{geom_beeswarm}}.}
}

\value{
  The interactive geometry based on \code{\link[ggbeeswarm]{geom_beeswarm}}.
}

\author{Patrick Aboyoun}

\examples{
# Create NanoStringRccSet from data files
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringNCTools")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")
solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)

eg_data <- as.data.frame(assayDataElement(solidTumor, "exprs")[1:5, 1])
eg_data[["tooltip"]] <- names(eg_data)
geom_beeswarm_interactive(aes_string(tooltip = "tooltip"),
                          data=eg_data)
}

\seealso{\code{\link[ggbeeswarm]{geom_beeswarm}}}

\keyword{graphics}
