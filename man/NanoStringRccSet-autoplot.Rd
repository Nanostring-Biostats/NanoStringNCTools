\name{NanoStringRccSet-autoplot}

\alias{autoplot}
\alias{NanoStringRccSet-autoplot}
\alias{autoplot.NanoStringRccSet}

\title{
Plot NanoStringRccSet Data
}
\description{
Generate common plots to visualize and QC NanoStringRccSet data.
}
\usage{
       autoplot(object,
           type = c("boxplot-feature",
                    "boxplot-signature",
                    "bindingDensity-mean",
                    "bindingDensity-sd",
                    "ercc-linearity",
                    "ercc-lod",
                    "heatmap-genes",
                    "heatmap-signatures",
                    "housekeep-geom",
                    "lane-bindingDensity",
                    "lane-fov",
                    "mean-sd-features",
                    "mean-sd-samples"),
           log2scale = TRUE,
           elt = "exprs",
           index = 1L,
           geomParams = list(),
           tooltipDigits = 4L,
           heatmapGroup = NULL,
           blacklist = NULL,
           tooltipID = dimLabels( object )[2L],
           qcCutoffs = list(
             Housekeeper = c("failingCutoff" = 32,"passingCutoff" = 100) ,
             Imaging = c("fovCutoff" = 0.75) ,
             BindingDensity = c("minimumBD" = 0.1, "maximumBD" = 2.25, "maximumBDSprint" = 1.8) ,
             ERCCLinearity = c("correlationValue" = 0.95) ,
             ERCCLoD = c("standardDeviations" = 2)
           ),
           scalingFactor=1L,
           show_rownames_gene_limit=60L,
           show_colnames_gene_limit=36L,
           show_rownames_sig_limit=60L,
           show_colnames_sig_limit=36L,
           subSet = NULL ,
           \ldots)
}
\arguments{
 \item{object}{A NanoStringRccSet object}
 \item{type}{Character string referencing the type of plot to generate}
 \item{log2scale}{An optional boolean indicating expression data is on log2 scale}
 \item{elt}{An optional character string of the expression matrix name}
 \item{index}{An optional integer giving the feature of interest row location}
 \item{geomParams}{An option \code{list} of parameters for geometry}
 \item{tooltipDigits}{An optional integer for number of tooltip decimal places to display}
 \item{heatmapGroup}{An optional character string referencing \code{pData} column to color samples by in heatmap}
 \item{blacklist}{An optional character vector of features not to plot}
 \item{tooltipID}{An optional character string referencing \code{pData} column to use for sample ID in the tooltip}
 \item{qcCutoffs}{An optional \code{list} of QC cutoffs}
 \item{scalingFactor}{An optional numeric value indicating a scaling factor to apply to plot drawing}
 \item{show_rownames_gene_limit}{An optional integer limit on number of features to display row-wise}
 \item{show_colnames_gene_limit}{An optional integer limit on number of features to display column-wise}
 \item{show_rownames_sig_limit}{An optional integer limit on number of signatures to display row-wise}
 \item{show_colnames_sig_limit}{An optional integer limit on number of signatures to display column-wise}
}
\value{
  A \code{ggplot} or \code{pheatmap} plot depending on the type of plot generated
}
\details{
   \describe{
      \item{\code{"boxplot-feature"}}{Generate feature boxplots}
      \item{\code{"boxplot-signature"}}{Generate signature boxplots}
      \item{\code{"bindingDensity-mean"}}{Plot binding density displayed as average expression}
      \item{\code{"bindingDensity-sd"}}{Plot binding density displayed as standard deviation of expression}
      \item{\code{"ercc-linearity"}}{Assess linearity of ERCCs}
      \item{\code{"ercc-lod"}}{Assess limit of detection based on ERCC expression}
      \item{\code{"heatmap-genes"}}{Generate a heatmap from feature expression}
      \item{\code{"heatmap-signatures"}}{Generate a heatmap from signature expression}
      \item{\code{"housekeep-geom"}}{Plot geometric mean of housekeeper genes}
      \item{\code{"lane-bindingDensity"}}{View binding density by lane}
      \item{\code{"lane-fov"}}{Assess image quality by lane}
      \item{\code{"mean-sd-features"}}{Plot mean versus standard deviation feature-wise}
      \item{\code{"mean-sd-samples"}}{Plot mean versus standard deviation sample-wise}
   }
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

# Assess experiment linearity
autoplot(solidTumor, "ercc-linearity")

# Plot a feature's expression across all samples
autoplot(solidTumor, "boxplot-feature", index=2)
}
\keyword{datasets}
