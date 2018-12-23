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

% accessors
\alias{sData}
\alias{sData,NanoStringRccSet-method}
\alias{svarLabels}
\alias{svarLabels,NanoStringRccSet-method}

% summary
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

% transforming
\alias{sweep,NanoStringRccSet-method}
\alias{transform,NanoStringRccSet-method}

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
    \code{"BarcodeClass"}, \code{"GeneName"}, and \code{"Accession"}.}
  \item{experimentData}{An optional \code{\link{MIAME}} instance with meta-data
    about the experiment.}
  \item{annotation}{A character string for the \code{"GeneRLF"}.}
  \item{protocolData}{An \code{\link{AnnotatedDataFrame}} containing columns
    \code{"FileVersion"}, \code{"SoftwareVersion"}, \code{"SystemType"},
    \code{"SampleID"}, \code{"Owner"}, \code{"Comments"}, \code{"Date"},
    \code{"SystemAPF"}, \code{"AssayType"},
    \code{"LaneID"}, \code{"FovCount"}, \code{"FovCounted"}, \code{"ScannerID"},
    \code{"StagePosition"}, \code{"BindingDensity"}, \code{"CartridgeID"}, and
    \code{"CartridgeBarcode"}.}
  \item{\ldots}{Additional arguments for \code{\link{ExpressionSet}}.}
}

\section{Accessors}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} accessor
  methods, NanoStringRccSet objects have the following:
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

\section{Summary}{
  \describe{
    \item{}{
      \code{summary(object, MARGIN, elt = "exprs")}: Summaries the features
      (\code{MARGIN = 1}) or samples (\code{MARGIN = 2}) of the \code{elt}
      assayData matrix.
    }
  }
}

\section{Subsetting}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} subsetting
  methods, NanoStringRccSet objects have the following:
  \describe{
    \item{}{
      \code{subset(x, subset, select, \ldots)}: Subset the features and samples
      using the \code{subset} and \code{select} arguments respectively.
    }
    \item{}{
      \code{endogenousSet(object)}: Extracts the endogenous barcode class
      subset of \code{object}.
    }
    \item{}{
      \code{housekeepingSet(object)}: Extracts the housekeeping barcode class
      subset of \code{object}.
    }
    \item{}{
      \code{negativeControlSet(object)}: Extracts the negative control barcode
      class subset of \code{object}.
    }
    \item{}{
      \code{positiveControlSet(object)}: Extracts the positive control barcode
      class subset of \code{object}.
    }
    \item{}{
      \code{controlSet(object)}: If the \code{"IsControl"} indicator column is
      present in the features, extracts the subset representing the controls of
      \code{object}.
    }
    \item{}{
      \code{nonControlSet(object)}: If the \code{"IsControl"} indicator column
      is present in the features, extracts the subset representing the
      non-controls of \code{object}.
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
      features (\code{MARGIN = 1}) or samples (\code{MARGIN = 2}) for all of
      \code{X}.
    }
    \item{}{
      \code{endogenousApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the features (\code{MARGIN = 1}) or samples (\code{MARGIN = 2})
      for the endogenous barcode class subset of \code{X}.
    }
    \item{}{
      \code{housekeepingApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the features (\code{MARGIN = 1}) or samples (\code{MARGIN = 2})
      for the housekeeping barcode class subset of \code{X}.
    }
    \item{}{
      \code{negativeControlApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the features (\code{MARGIN = 1}) or samples (\code{MARGIN = 2})
      for the negative control barcode class subset of \code{X}.
    }
    \item{}{
      \code{positiveControlApply(X, MARGIN, FUN, \ldots, elt = "exprs")}: Loop
      over the features (\code{MARGIN = 1}) or samples (\code{MARGIN = 2})
      for the positive control barcode class subset of \code{X}.
    }
  }
}

\section{Transforming}{
  Each of the following transformation methods adds or modifies
  \code{assayData} matrices:
  \describe{
    \item{}{
      \code{sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, \ldots,
         fromElt = "exprs", toElt)}: Similar to the \code{\link[base]{sweep}}
         function in the \pkg{base} package, creates a \code{toElt} assayData
         matrix by sweeping out a summary statistic from the \code{fromElt}
         assayData matrix.
    }
    \item{}{
      \code{transform(`_data`, \ldots)}: Similar to the
      \code{\link[base]{transform}} generic in the \pkg{base} package, creates
      or modifies one or more assayData matrices based upon \code{name = value}
      pairs in \code{\ldots}.
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
}
\keyword{classes}
\keyword{methods}
