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
\alias{getSigFuncs}
\alias{getSigFuncs,SignatureSet-method}
\alias{groups}
\alias{groups,SignatureSet-method}
\alias{groups<-}
\alias{groups<-,SignatureSet,ANY-method}
\alias{groups<-,SignatureSet,NULL-method}
\alias{groups<-,SignatureSet,factor-method}
\alias{version}
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
  SignatureSet(weights = NumericList(), groups = factor(), func = character(), 
               version = character(), \ldots)
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
      \item{length(x)}{returns the number of signatures in \code{x}.
      }
      \item{lengths(x, use.names = TRUE)}{returns a named integer vector
        containing the number of genes in each of the signatures in \code{x}.
      }
      \item{names(x)}{ returns a character vector containing the signature
        names in \code{x}.
      }
      \item{weights(object)}{returns a named \code{\link{NumericList}} that
        defines the linear combination based signatures.
      }
      \item{weights(object) <- value}{replaces the \code{\link{NumericList}}
        that defines the linear combination based signatures.
      }
      \item{getSigFuncs(object)}{returns the signature functions of an object.
      }
      \item{groups(object)}{returns a factor vector representing the signature groups.
      }
      \item{groups(object) <- value}{replaces the factor vector representing the 
        signature groups.
      }
      \item{version(object)}{
        \code{}: returns the signature version.
      }
    }
}

\author{Patrick Aboyoun}

\seealso{
  \code{\link{NanoStringRccSet}}
}

\examples{
SignatureSet(weights=list(x = c(a = 1),
                          y = c(b = 1/3, d = 2/3),
                          z = c(a = 2, c = 4)),
             groups=factor("x", "y", "z"),
             func = c(x="default", y="default", z="default"))
}

\keyword{classes}
\keyword{methods}
