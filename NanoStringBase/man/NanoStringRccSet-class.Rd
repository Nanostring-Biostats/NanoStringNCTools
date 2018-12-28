\name{NanoStringRccSet-class}
\docType{class}

\alias{class:NanoStringRccSet}
\alias{NanoStringRccSet-class}

% constructors
\alias{NanoStringRccSet}
\alias{NanoStringRccSet,environment-method}
\alias{NanoStringRccSet,matrix-method}
\alias{NanoStringRccSet,missing-method}
\alias{NanoStringRccSet,NanoStringRccSet-method}

% accessing
\alias{sData}
\alias{sData,NanoStringRccSet-method}
\alias{svarLabels}
\alias{svarLabels,NanoStringRccSet-method}

% summarizing
\alias{summary,NanoStringRccSet-method}

% subsetting
\alias{subset,NanoStringRccSet-method}
\alias{endogenousSet}
\alias{endogenousSet,NanoStringRccSet-method}
\alias{housekeepingSet}
\alias{housekeepingSet,NanoStringRccSet-method}
\alias{negativeControlSet}
\alias{negativeControlSet,NanoStringRccSet-method}
\alias{positiveControlSet}
\alias{positiveControlSet,NanoStringRccSet-method}
\alias{controlSet}
\alias{controlSet,NanoStringRccSet-method}
\alias{nonControlSet}
\alias{nonControlSet,NanoStringRccSet-method}

% looping
\alias{esApply,NanoStringRccSet-method}
\alias{endogenousApply}
\alias{endogenousApply,NanoStringRccSet-method}
\alias{housekeepingApply}
\alias{housekeepingApply,NanoStringRccSet-method}
\alias{negativeControlApply}
\alias{negativeControlApply,NanoStringRccSet-method}
\alias{positiveControlApply}
\alias{positiveControlApply,NanoStringRccSet-method}
\alias{esBy}
\alias{esBy,NanoStringRccSet-method}

% transforming
\alias{transform,NanoStringRccSet-method}

% evaluating
\alias{with,NanoStringRccSet-method}

\title{Class to Contain NanoString Expression Level Assays}

\description{
  The \code{NanoStringRccSet} class extends the
  \code{\linkS4class{ExpressionSet}} class for NanoString Reporter Code Count
  (RCC) data.
}

\usage{
NanoStringRccSet(assayData,
                 phenoData = annotatedDataFrameFrom(assayData, byrow = FALSE),
                 featureData = annotatedDataFrameFrom(assayData, byrow = TRUE),
                 experimentData = MIAME(),
                 annotation = character(),
                 protocolData = annotatedDataFrameFrom(assayData, byrow = FALSE),
                 \ldots)
}

\arguments{
  \item{assayData}{A \code{matrix} or \code{environment} containing the RCCs.}
  \item{phenoData}{An \code{\link{AnnotatedDataFrame}} containing the phenotypic
    data.}
  \item{featureData}{An \code{\link{AnnotatedDataFrame}} containing columns
    \code{"BarcodeClass"}, \code{"GeneName"}, \code{"Accession"},
    \code{"IsControl"}, and \code{"ControlConc"}.}
  \item{experimentData}{An optional \code{\link{MIAME}} instance with meta-data
    about the experiment.}
  \item{annotation}{A character string for the \code{"GeneRLF"}.}
  \item{protocolData}{An \code{\link{AnnotatedDataFrame}} containing columns
    \code{"FileVersion"}, \code{"SoftwareVersion"}, \code{"SystemType"},
    \code{"SampleID"}, \code{"SampleOwner"}, \code{"SampleComments"},
    \code{"SampleDate"}, \code{"SystemAPF"}, \code{"AssayType"},
    \code{"LaneID"}, \code{"FovCount"}, \code{"FovCounted"}, \code{"ScannerID"},
    \code{"StagePosition"}, \code{"BindingDensity"}, \code{"CartridgeID"}, and
    \code{"CartridgeBarcode"}.}
  \item{\ldots}{Additional arguments for \code{\link{ExpressionSet}}.}
}

\section{Accessing}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} accessor
  methods, \code{NanoStringRccSet} objects have the following:
  \describe{
    \item{}{
      \code{sData(object)}: extracts the \code{data.frame} containing the
      sample data, \code{cbind(pData(object), pData(protocolData(object)))}
    }
    \item{}{
      \code{svarLabels(object)}: extracts the sample data column names,
      \code{c(varLabels(object), varLabels(protocolData(object)))}
    }
  }
}

\section{Summarizing}{
  \describe{
    \item{}{
      \code{summary(object, MARGIN, GROUP = NULL, elt = "exprs")}: Marginal
      summaries (Count, Non-Missing Count, Mean, Standard Deviation, Skewness,
      Excess Kurtosis, Minimum, First Quartile, Median, Third Quartile,
      Maximum, Median Absolute Deviation, Median Polish Row / Column Effect)
      of the \code{elt} \code{\link[Biobase]{assayData}} matrix along either
      the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension.
    }
  }
}

\section{Subsetting}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} subsetting
  methods, \code{NanoStringRccSet} objects have the following:
  \describe{
    \item{}{
      \code{subset(x, subset, select, \ldots)}: Subset the feature and sample
      dimensions using the \code{subset} and \code{select} arguments
      respectively. The \code{subset} argument will be evaluated with
      respect to the \code{\link[Biobase]{featureData}}, while the
      \code{select} argument will be evaluated with respect to the
      \code{\link[Biobase]{phenoData}} and \code{\link[Biobase]{protocolData}}.
    }
    \item{}{
      \code{endogenousSet(object)}: Extracts the endogenous barcode class
      feature subset of \code{object}.
    }
    \item{}{
      \code{housekeepingSet(object)}: Extracts the housekeeping barcode class
      feature subset of \code{object}.
    }
    \item{}{
      \code{negativeControlSet(object)}: Extracts the negative control barcode
      class feature subset of \code{object}.
    }
    \item{}{
      \code{positiveControlSet(object)}: Extracts the positive control barcode
      class feature subset of \code{object}.
    }
    \item{}{
      \code{controlSet(object)}: If the \code{"IsControl"} indicator column is
      present in the feature data, extracts the feature subset representing the
      controls of \code{object}.
    }
    \item{}{
      \code{nonControlSet(object)}: If the \code{"IsControl"} indicator column
      is present in the feature data, extracts the feature subset representing
      the non-controls of \code{object}.
    }
  }
}

\section{Looping}{
  In the looping methods below, the \code{X}, \code{MARGIN}, \code{FUN}, and
  \code{\ldots} arguments follow the definition in the
  \code{\link[Biobase]{esApply}} generic. The additional \code{elt} argument
  specifies which \code{assayData} matrix to loop over.
  \describe{
    \item{}{
      \code{esApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop over the
      feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension
      for all of \code{X}.
    }
    \item{}{
      \code{endogenousApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2})
      dimension for the endogenous barcode class feature subset of \code{X}.
    }
    \item{}{
      \code{housekeepingApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2})
      dimension for the housekeeping barcode class feature subset of \code{X}.
    }
    \item{}{
      \code{negativeControlApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2})
      dimension for the negative control barcode class feature subset of
      \code{X}.
    }
    \item{}{
      \code{positiveControlApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2})
      dimension for the positive control barcode class feature subset of
      \code{X}.
    }
    \item{}{
      \code{esBy(X, GROUP, FUN, \ldots, simplify = TRUE)}: Split
      \code{X} by \code{GROUP} column within \code{featureData},
      \code{phenoData}, or \code{protocolData} and apply \code{FUN}
      to each partition.
    }
  }
}

\section{Transforming}{
  \describe{
    \item{}{
      \code{transform(`_data`, \ldots)}: Similar to the
      \code{\link[base]{transform}} generic in the \pkg{base} package, creates
      or modifies one or more \code{\link[Biobase]{assayData}} matrices based
      upon \code{name = value} pairs in \code{\ldots}. The transforms specified
      in \code{\ldots} are appended to the preprocessing list in
      \code{\link{experimentData}}, which can be extracted using the
      \code{\link{preproc}} method.
    }
  }
}

\section{Evaluating}{
  \describe{
    \item{}{
      \code{with(data, expr, \ldots)}: Evaluate expression \code{expr} with
      respect to \code{\link[Biobase]{assayData}},
      \code{\link[Biobase]{featureData}}, \code{\link[Biobase]{phenoData}},
      and \code{\link[Biobase]{protocolData}};
      \code{c(as.list(assayData(data)), fData(data), sData(data))}.
    }
  }
}

\author{Patrick Aboyoun}

\seealso{
  \code{\link{readNanoStringRccSet}},
  \code{\link{writeNanoStringRccSet}},
  \code{\linkS4class{ExpressionSet}}
}

\examples{
# Create NanoStringRccSet from data files
datadir <- system.file("extdata", "3D_Bio_Example_Data",
                       package = "NanoStringBase")
rccs <- dir(datadir, pattern = "SKMEL.*\\\\.RCC$", full.names = TRUE)
rlf <- file.path(datadir, "3D_SolidTumor_Sig.rlf")
pheno <- file.path(datadir, "3D_SolidTumor_PhenoData.csv")
solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = rlf, phenoDataFile = pheno)


# Create a deep copy of a NanoStringRccSet object
deepCopy <- NanoStringRccSet(solidTumor)
all.equal(solidTumor, deepCopy)
identical(solidTumor, deepCopy)


# Accessing sample data and column names
head(sData(solidTumor))
svarLabels(solidTumor)


# Marginal summarizing of NanoStringRccSet assayData matrices
head(summary(solidTumor, 1)) # Marginal summaries along features
head(summary(solidTumor, 2)) # Marginal summaries along samples


# Subsetting NanoStringRccSet objects
# Extract the positive controls for wildtype BRAF
dim(solidTumor)
dim(subset(solidTumor, BarcodeClass == "Positive", BRAFGenotype == "wt/wt"))

# Extract by barcode class
with(solidTumor, table(BarcodeClass))
with(endogenousSet(solidTumor), table(BarcodeClass))
with(housekeepingSet(solidTumor), table(BarcodeClass))
with(negativeControlSet(solidTumor), table(BarcodeClass))
with(positiveControlSet(solidTumor), table(BarcodeClass))
with(controlSet(solidTumor), table(BarcodeClass))
with(nonControlSet(solidTumor), table(BarcodeClass))


# Looping over NanoStringRccSet assayData matrices
log1pCoefVar <- function(x){
  x <- log1p(x)
  sd(x) / mean(x)
}

# Log1p Coefficient of Variation along Features
head(esApply(solidTumor, 1, log1pCoefVar))
head(endogenousApply(solidTumor, 1, log1pCoefVar))
head(housekeepingApply(solidTumor, 1, log1pCoefVar))
head(negativeControlApply(solidTumor, 1, log1pCoefVar))
head(positiveControlApply(solidTumor, 1, log1pCoefVar))

# Log1p Coefficient of Variation along Samples
head(esApply(solidTumor, 2, log1pCoefVar))
head(endogenousApply(solidTumor, 2, log1pCoefVar))
head(housekeepingApply(solidTumor, 2, log1pCoefVar))
head(negativeControlApply(solidTumor, 2, log1pCoefVar))
head(positiveControlApply(solidTumor, 2, log1pCoefVar))


# Transforming NanoSetRccSet assayData matrices
# Subtract max count from each sample
# Create log1p transformation of adjusted counts
thresh <- negativeControlApply(solidTumor, 2, max)
solidTumor2 <-
  transform(solidTumor,
            negCtrlZeroed = sweep(exprs, 2, thresh),
            log1p_negCtrlZeroed = log1p(pmax(negCtrlZeroed, 0)))
assayDataElementNames(solidTumor2)


# Evaluating expression using NanoStringRccSet data
meanLog1pExprs <-
  with(solidTumor,
       {
         means <- split(apply(exprs, 1, function(x) mean(log1p(x))),
                        BarcodeClass)
         means <- means[order(sapply(means, median))]
         boxplot(means, horizontal = TRUE)
         means
       })
}

\keyword{classes}
\keyword{methods}
