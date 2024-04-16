\name{NanoStringRccSet-class}
\docType{class}

\alias{class:NanoStringRccSet}
\alias{NanoStringRccSet-class}

% constructors
\alias{NanoStringRccSet}
\alias{NanoStringRccSet,environment-method}
\alias{NanoStringRccSet,ExpressionSet-method}
\alias{NanoStringRccSet,matrix-method}
\alias{NanoStringRccSet,missing-method}
\alias{NanoStringRccSet,NanoStringRccSet-method}

% coercion
\alias{coerce,ExpressionSet,NanoStringRccSet-method}

% show
\alias{show,NanoStringRccSet-method}

% accessing
\alias{sData}
\alias{sData,NanoStringRccSet-method}
\alias{svarLabels}
\alias{svarLabels,NanoStringRccSet-method}
\alias{dimLabels,NanoStringRccSet-method}
\alias{dimLabels<-,NanoStringRccSet,character-method}
\alias{signatures}
\alias{signatures,NanoStringRccSet-method}
\alias{signatures<-}
\alias{signatures<-,NanoStringRccSet,SignatureSet-method}
\alias{signatureScores}
\alias{signatureScores,NanoStringRccSet-method}
\alias{signatureFuncs}
\alias{signatureFuncs,NanoStringRccSet-method}
\alias{signatureGroups}
\alias{signatureGroups,NanoStringRccSet-method}
\alias{setSignatureFuncs<-}
\alias{setSignatureFuncs<-,NanoStringRccSet,character-method}
\alias{setSignatureGroups<-}
\alias{setSignatureGroups<-,NanoStringRccSet,character-method}
\alias{setSignatureGroups<-,NanoStringRccSet,factor-method}
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
\alias{assayDataApply}
\alias{assayDataApply,NanoStringRccSet-method}
\alias{signatureScoresApply}
\alias{signatureScoresApply,NanoStringRccSet-method}
\alias{esBy}
\alias{esBy,NanoStringRccSet-method}

% transforming
\alias{munge}
\alias{munge,NanoStringRccSet-method}
\alias{transform,NanoStringRccSet-method}

% evaluating
\alias{with,NanoStringRccSet-method}

% quality control
\alias{setQCFlags,NanoStringRccSet-method}

% normalization
\alias{normalize,NanoStringRccSet-method}

% plotting
\alias{ggplot.NanoStringRccSet}
\alias{update_geom_params}

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
                 dimLabels = c("GeneName", "SampleID"),
                 signatures = SignatureSet(),
                 design = NULL,
                 \ldots)
}

\arguments{
  \item{assayData}{A \code{matrix} or \code{environment} containing the RCCs.}
  \item{phenoData}{An \code{\link{AnnotatedDataFrame}} containing the phenotypic
    data.}
  \item{featureData}{An \code{\link{AnnotatedDataFrame}} containing columns
    \code{"CodeClass"}, \code{"GeneName"}, \code{"Accession"},
    \code{"IsControl"}, and \code{"ControlConc"}.}
  \item{experimentData}{An optional \code{\link{MIAME}} instance with meta-data
    about the experiment.}
  \item{annotation}{A \code{character} string for the \code{"GeneRLF"}.}
  \item{protocolData}{An \code{\link{AnnotatedDataFrame}} containing columns
    \code{"FileVersion"}, \code{"SoftwareVersion"}, \code{"SystemType"},
    \code{"SampleID"}, \code{"SampleOwner"}, \code{"SampleComments"},
    \code{"SampleDate"}, \code{"SystemAPF"}, \code{"AssayType"},
    \code{"LaneID"}, \code{"FovCount"}, \code{"FovCounted"}, \code{"ScannerID"},
    \code{"StagePosition"}, \code{"BindingDensity"}, \code{"CartridgeID"}, and
    \code{"CartridgeBarcode"}.}
  \item{dimLabels}{A \code{character} vector of length 2 that provides the
    column names to use as labels for the features and samples respectively
    in the \code{autoplot} method.}
  \item{signatures}{An optional \code{\link{SignatureSet}} object containing
    signature definitions.}
  \item{design}{An optional one-sided formula representing the experimental
    design based on columns from \code{\link[Biobase]{phenoData}}}
  \item{\ldots}{Additional arguments for \code{\link{ExpressionSet}}.}
}

\section{Accessing}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} accessor
  methods, \code{NanoStringRccSet} objects have the following:
  \describe{
    \item{sData(object)}{ extracts the \code{data.frame} containing the
      sample data, \code{cbind(pData(object), pData(protocolData(object)))}.
    }
    \item{svarLabels(object)}{ extracts the sample data column names,
      \code{c(varLabels(object), varLabels(protocolData(object)))}.
    }
    \item{dimLabels(object)}{ extracts the column names to use as labels
      for the features and samples in the \code{autoplot} method.
    }
    \item{dimLabels(object) <- value}{ replaces the \code{dimLabels} of
      the \code{object}.
    }
    \item{signatures(object)}{ extracts the \code{\link{SignatureSet}}
      of the \code{object}.
    }
    \item{signatures(object) <- value}{ replaces the
      \code{\link{SignatureSet}} of the \code{object}.
    }
    \item{signatureScores(object, elt = "exprs")}{ extracts the matrix
      of computed signature scores.
    }
    \item{design(object)}{ extracts the one-sided formula representing
      the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
    \item{design(object) <- value}{ replaces the one-sided formula
      representing the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
    \item{setSignatureFuncs(object)}{ returns the signature functions.
    }
    \item{setSignatureFuncs(object) <- value}{ replaces the signature functions.
    }
    \item{setSignatureGroups(object) <- value}{ returns the signature groups.
    }
    \item{setSignatureGroups(object) <- value}{ replaces the signature groups.
    }
  }
}

\section{Summarizing}{
  \describe{
    \item{summary(object, MARGIN = 2L, GROUP = NULL, log2scale = TRUE, elt = "exprs", signatureScores = FALSE)}{
      When \code{signatureScores = FALSE}, the marginal summaries of the
      \code{elt} \code{\link[Biobase]{assayData}} matrix along either the
      feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension.

      When \code{signatureScores = TRUE}, the marginal summaries of the
      \code{elt} \code{signatureScores} matrix along either the
      signature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension.

      When \code{log2scale = FALSE}, the summary statistics are Mean, Standard
      Deviation, Skewness, Excess Kurtosis, Minimum, First Quartile, Median,
      Third Quartile, and Maximum.

      When \code{log2scale = TRUE}, the summary statistics are Geometric Mean
      with thresholding at 0.5, Size Factor
      (\code{2^(`MeanLog2` - mean(`MeanLog2`))}), Mean of Log2 with
      thresholding at 0.5, Standard Deviation of Log2 with thresholding at 0.5,
      Minimum, First Quartile, Median, Third Quartile, and Maximum.
    }
  }
}

\section{Subsetting}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} subsetting
  methods, \code{NanoStringRccSet} objects have the following:
  \describe{
    \item{subset(x, subset, select, \ldots)}{ Subset the feature and sample
      dimensions using the \code{subset} and \code{select} arguments
      respectively. The \code{subset} argument will be evaluated with
      respect to the \code{\link[Biobase]{featureData}}, while the
      \code{select} argument will be evaluated with respect to the
      \code{\link[Biobase]{phenoData}} and \code{\link[Biobase]{protocolData}}.
    }
    \item{endogenousSubset(x, subset, select)}{ Extracts the endogenous
      barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{housekeepingSubset(x, subset, select)}{ Extracts the housekeeping
      barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{negativeControlSubset(x, subset, select)}{ Extracts the negative
      control barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{positiveControlSubset(x, subset, select)}{Extracts the positive
      control barcode class feature subset of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{controlSubset(x, subset, select)}{Extracts the feature subset
      representing the controls of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{nonControlSubset(x, subset, select)}{Extracts the feature subset
      representing the non-controls of \code{x} with optional additional
      subsetting using \code{subset} and \code{select}.
    }
    \item{signatureSubset(x, subset, select)}{Extracts the feature subset
      representing the genes in the \code{signatures} of \code{x} with
      optional additional subsetting using \code{subset} and \code{select}.
    }
  }
}

\section{Looping}{
  \describe{
    \item{assayDataApply(X, MARGIN, FUN, \ldots, elt = "exprs")}{Loop over
      the feature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2}) dimension
      of \code{assayDataElement(X, elt)}.
    }
    \item{signatureScoresApply(X, MARGIN, FUN, \ldots, elt = "exprs")}{Loop
      over the signature (\code{MARGIN = 1}) or sample (\code{MARGIN = 2})
      dimension of \code{signatureScores(X, elt)}.
    }
    \item{esBy(X, GROUP, FUN, \ldots, simplify = TRUE)}{Split
      \code{X} by \code{GROUP} column within \code{featureData},
      \code{phenoData}, or \code{protocolData} and apply \code{FUN}
      to each partition.
    }
  }
}

\section{Transforming}{
  \describe{
    \item{munge(data, mapping = update(design(data), exprs ~ .), extradata = NULL, elt = "exprs", \ldots)}{munge argument \code{data} into a data.frame object for modeling and
      visualization using the \code{mapping} argument. Supplemental data can be
      specified using the \code{extradata} argument.
    }
    \item{transform(`_data`, \ldots)}{Similar to the
      \code{\link[base]{transform}} generic in the \pkg{base} package, creates
      or modifies one or more \code{\link[Biobase]{assayData}} matrices based
      upon \code{name = value} pairs in \code{\ldots}. The expressions in
      \code{\ldots} are appended to the preprocessing list in
      \code{\link{experimentData}}, which can be extracted using the 
      \code{\link{preproc}} method.
    }
  }
}

\section{Evaluating}{
  \describe{
    \item{with(data, expr, \ldots)}{Evaluate expression \code{expr} with
      respect to \code{\link[Biobase]{assayData}},
      \code{\link[Biobase]{featureData}}, \code{\link[Biobase]{phenoData}},
      and \code{\link[Biobase]{protocolData}};
      \code{c(as.list(assayData(data)), fData(data), sData(data))}.
    }
  }
}

\section{Normalizing}{
  \describe{
    \item{normalize(object, type, fromElt = "exprs", toElt = "exprs_norm", \ldots)}{
    }
  }
}

\section{Plotting}{
  \describe{
    \item{ggplot(data, mapping = aes(), \ldots, extradata = NULL, tooltip_digits = 4L, environment = parent.frame())}{
      the \code{NanoStringRccSet} method for \code{ggplot}.
    }
    \item{autoplot(object, type, log2scale = TRUE, elt = "exprs", index = 1L, geomParams = list(), tooltipDigits = 4L, heatmapGroup = NULL, \ldots)}{
    }
  }
}

\value{
  An S4 class containing NanoString Expression Level Assays
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
                       package = "NanoStringNCTools")
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
munge(solidTumor)


# Marginal summarizing of NanoStringRccSet assayData matrices
head(summary(solidTumor, 1)) # Marginal summaries along features
head(summary(solidTumor, 2)) # Marginal summaries along samples


# Subsetting NanoStringRccSet objects
# Extract the positive controls for wildtype BRAF
dim(solidTumor)
dim(subset(solidTumor, CodeClass == "Positive", BRAFGenotype == "wt/wt"))

# Extract by barcode class
with(solidTumor, table(CodeClass))
with(endogenousSubset(solidTumor), table(CodeClass))
with(housekeepingSubset(solidTumor), table(CodeClass))
with(negativeControlSubset(solidTumor), table(CodeClass))
with(positiveControlSubset(solidTumor), table(CodeClass))
with(controlSubset(solidTumor), table(CodeClass))
with(nonControlSubset(solidTumor), table(CodeClass))


# Looping over NanoStringRccSet assayData matrices
log1pCoefVar <- function(x){
  x <- log1p(x)
  sd(x) / mean(x)
}

# Log1p Coefficient of Variation along Features
head(assayDataApply(solidTumor, 1, log1pCoefVar))

# Log1p Coefficient of Variation along Samples
head(assayDataApply(solidTumor, 2, log1pCoefVar))


# Transforming NanoSetRccSet assayData matrices
# Subtract max count from each sample
# Create log1p transformation of adjusted counts
thresh <- assayDataApply(negativeControlSubset(solidTumor), 2, max)
solidTumor2 <-
  transform(solidTumor,
            negCtrlZeroed = sweep(exprs, 2, thresh),
            log1p_negCtrlZeroed = log1p(pmax(negCtrlZeroed, 0)))
assayDataElementNames(solidTumor2)


# Evaluating expression using NanoStringRccSet data
meanLog1pExprs <-
  with(solidTumor,
       {
         means <- split(apply(exprs, 1, function(x) mean(log1p(x))), CodeClass)
         means <- means[order(sapply(means, median))]
         boxplot(means, horizontal = TRUE)
         means
       })
}

\keyword{classes}
\keyword{methods}
