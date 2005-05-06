\name{summary.SBML}
\alias{summary.SBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get summary data frames from an SBML model}
\description{
  This function extracts species and reaction information from a model of class \code{SBML} and tabularizes them as two R data frames. 
}
\usage{
summary.SBML(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A model object of class \code{SBML} from which information is to be extracted.}
  \item{...}{ For compatibility with \code{summary} of the base package.}
 }
%\details{no details
     
}
\value{ A list containing the following named dataframe elements
 \item{species}{ Species information (i.e. names, ICs, BCs, and compartments) as a data frame.}
  \item{reactions}{ Reaction information tabularized as a dataframe, including string laws and initial fluxes.}
}

\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}

\seealso{ \code{\link{==.SBML}}  }
\examples{
library(SBMLR)  
curto=readSBMLR(file.path(.path.package("SBMLR"), "models/Curto.r"))  
summary(curto)

}
\keyword{ math }% at least one, from doc/KEYWORDS
