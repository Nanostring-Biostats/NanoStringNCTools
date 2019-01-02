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

% show
\alias{show,NanoStringRccSet-method}

% accessing
\alias{sData}
\alias{sData,NanoStringRccSet-method}
\alias{svarLabels}
\alias{svarLabels,NanoStringRccSet-method}
\alias{signatureWeights}
\alias{signatureWeights,NanoStringRccSet-method}
\alias{signatureWeights<-}
\alias{signatureWeights<-,NanoStringRccSet,ANY-method}
\alias{signatureWeights<-,NanoStringRccSet,NULL-method}
\alias{signatureWeights<-,NanoStringRccSet,NumericList-method}
\alias{design,NanoStringRccSet-method}
\alias{design<-,NanoStringRccSet,ANY-method}
\alias{design<-,NanoStringRccSet,NULL-method}
\alias{design<-,NanoStringRccSet,formula-method}

% summarizing
\alias{summary,NanoStringRccSet-method}

% subsetting
\alias{[,NanoStringRccSet-method}
\alias{subset,NanoStringRccSet-method}
\alias{endogenousSubset}
\alias{endogenousSubset,NanoStringRccSet-method}
\alias{housekeepingSubset}
\alias{housekeepingSubset,NanoStringRccSet-method}
\alias{negativeControlSubset}
\alias{negativeControlSubset,NanoStringRccSet-method}
\alias{positiveControlSubset}
\alias{positiveControlSubset,NanoStringRccSet-method}
\alias{controlSubset}
\alias{controlSubset,NanoStringRccSet-method}
\alias{nonControlSubset}
\alias{nonControlSubset,NanoStringRccSet-method}
\alias{signatureSubset}
\alias{signatureSubset,NanoStringRccSet-method}

% looping
\alias{esApply,NanoStringRccSet-method}
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
                 signatureWeights = NumericList(),
                 design = NULL,
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
  \item{signatureWeights}{An optional named \code{\link[IRanges]{NumericList}}
    containing named gene signature weights.}
  \item{design}{An optional one-sided formula representing the experimental
    design based on columns from \code{\link[Biobase]{phenoData}}}
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
    \item{}{
      \code{signatureWeights(object)}: extracts the list of gene signature
      weights.
    }
    \item{}{
      \code{signatureWeights(object) <- value}: replaces the list of gene
      signature weights.
    }
    \item{}{
      \code{design(object)}: extracts the one-sided formula representing
      the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
    \item{}{
      \code{design(object) <- value}: replaces the one-sided formula
      representing the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
  }
}

\section{Summarizing}{
  \describe{
    \item{}{
      \code{summary(object, MARGIN = 2L, GROUP = NULL, elt = "exprs")}:
      Marginal summaries (Geometric Mean with thresholding at 0.5, Robust
      Regression Size Factor, Median Polish-based Size Factor (see below),
      Mean, Standard Deviation, Mean of Log with thresholding at 0.5,
      Standard Deviation of Log with thresholding at 0.5, Minimum, First
      Quartile, Median, Third Quartile, and Maximum) of the \code{elt}
      \code{\link[Biobase]{assayData}} matrix along either the feature
      (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension.

      The Median Polish-based Size Factor is calculated using three steps:
      \enumerate{
        \item Log transform expressions with thresholding at 0.5
        \item Calculate median polish row / column effects
        \item Exponentiate these effects
      }
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
      \code{endogenousSubset(x, subset, select)}: Extracts the endogenous
      barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{housekeepingSubset(x, subset, select)}: Extracts the housekeeping
      barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{negativeControlSubset(x, subset, select)}: Extracts the negative
      control barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{positiveControlSubset(x, subset, select)}: Extracts the positive
      control barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{controlSubset(x, subset, select)}: Extracts the feature subset
      representing the controls of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{nonControlSubset(x, subset, select)}: Extracts the feature subset
      representing the non-controls of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{}{
      \code{signatureSubset(x, subset, select)}: Extracts the feature subset
      representing the genes in the \code{signatureWeights} of \code{x} with
      optional additional subsetting using \code{subset} and \code{select}.
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
      upon \code{name = value} pairs in \code{\ldots}. The transformed matrices
      created in \code{\ldots} are rounded to the nearest integer and the calls
      are appended to the preprocessing list in \code{\link{experimentData}},
      which can be extracted using the \code{\link{preproc}} method.
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


# Set experimental design
design(solidTumor) <- ~ BRAFGenotype + Treatment
design(solidTumor)


# Marginal summarizing of NanoStringRccSet assayData matrices
head(summary(solidTumor, 1)) # Marginal summaries along features
head(summary(solidTumor, 2)) # Marginal summaries along samples


# Subsetting NanoStringRccSet objects
# Extract the positive controls for wildtype BRAF
dim(solidTumor)
dim(subset(solidTumor, BarcodeClass == "Positive", BRAFGenotype == "wt/wt"))

# Extract by barcode class
with(solidTumor, table(BarcodeClass))
with(endogenousSubset(solidTumor), table(BarcodeClass))
with(housekeepingSubset(solidTumor), table(BarcodeClass))
with(negativeControlSubset(solidTumor), table(BarcodeClass))
with(positiveControlSubset(solidTumor), table(BarcodeClass))
with(controlSubset(solidTumor), table(BarcodeClass))
with(nonControlSubset(solidTumor), table(BarcodeClass))


# Looping over NanoStringRccSet assayData matrices
log1pCoefVar <- function(x){
  x <- log1p(x)
  sd(x) / mean(x)
}

# Log1p Coefficient of Variation along Features
head(esApply(solidTumor, 1, log1pCoefVar))

# Log1p Coefficient of Variation along Samples
head(esApply(solidTumor, 2, log1pCoefVar))


# Transforming NanoSetRccSet assayData matrices
# Subtract max count from each sample
# Create log1p transformation of adjusted counts
thresh <- esApply(negativeControlSubset(solidTumor), 2, max)
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
