\name{SignatureSet-class}
\docType{class}

\alias{class:SignatureSet}
\alias{SignatureSet-class}

% constructor
\alias{SignatureSet}

% show
\alias{show,SignatureSet-method}

% utilities
\alias{length,SignatureSet-method}
\alias{lengths,SignatureSet-method}
\alias{names,SignatureSet-method}
\alias{weights,SignatureSet-method}
\alias{weights<-}
\alias{weights<-,SignatureSet,ANY-method}
\alias{weights<-,SignatureSet,NULL-method}
\alias{weights<-,SignatureSet,NumericList-method}
\alias{weights<-,SignatureSet,list-method}
\alias{weights<-,SignatureSet,CompressedNumericList-method}
\alias{getSigFuncs,SignatureSet-method}
\alias{groups,SignatureSet-method}
\alias{groups<-}
\alias{groups<-,SignatureSet,ANY-method}
\alias{groups<-,SignatureSet,NULL-method}
\alias{groups<-,SignatureSet,factor-method}
\alias{version,SignatureSet-method}
\alias{version<-}
\alias{version<-,SignatureSet,ANY-method}
\alias{version<-,SignatureSet,NULL-method}
\alias{version<-,SignatureSet,character-method}

\title{Class to Contain Signature Definitions}

\description{
  The \code{SignatureSet} class defines gene-based signatures.
}

\usage{
  SignatureSet(weights = NumericList(), groups = "factor", func = "character", version = "character", \ldots)
}

\arguments{
  \item{weights}{A named \code{NumericList} defining signatures based on linear
    combinations of genes.}
  \item{groups}{A factor vector indicating groups in the \code{SignatureSet}}
  \item{func}{Character indicating function to use}
  \item{version}{Character indicating version to use}
  \item{\ldots}{Additional arguments for future use.}
}

\value{
  A \code{SignatureSet} object
}

\section{Utilities}{
    \describe{
      \item{}{
        \code{length(x)}: returns the number of signatures in \code{x}.
      }
      \item{}{
        \code{lengths(x, use.names = TRUE)}: returns a named integer vector
        containing the number of genes in each of the signatures in \code{x}.
      }
      \item{}{
        \code{names(x)}: returns a character vector containing the signature
        names in \code{x}.
      }
      \item{}{
        \code{weights(object)}: returns a named \code{\link{NumericList}} that
        defines the linear combination based signatures.
      }
      \item{}{
        \code{weights(object) <- value}: replaces the \code{\link{NumericList}}
        that defines the linear combination based signatures.
      }
    }
}

\author{Patrick Aboyoun}

\seealso{
  \code{\link{NanoStringRccSet}}
}

\examples{
}

\keyword{classes}
\keyword{methods}
