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

\seealso{\code{\link[ggbeeswarm]{geom_beeswarm}}}

\examples{
}

\keyword{graphics}
