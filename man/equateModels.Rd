\name{==.SBML}
\alias{==.SBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Check the equality of the summary data frames of two SBML models}
\description{
  This function tests the equivalence of two models with respect to their summary data frames.  
}
\usage{
model1==model2
"=="(model1,model2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model1}{The first of the two model objects of class SBML which are to be compared.}
  \item{model2}{The second model object.}
 }
% \details{no details}

\value{ A list containing the following named dataframe elements as from \code{summary}, but now with boolean entries
 \item{species}{ The equality of species information tabularized  as a data frame.}
  \item{reactions}{ The equality of reaction information tabularized as a dataframe.}
}

\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}

\seealso{ \code{\link{summary.SBML}}  }
\examples{
library(SBMLR)  
curto1=readSBMLR(file.path(.path.package("SBMLR"), "models/curto.r"))  
curto2=readSBML(file.path(.path.package("SBMLR"), "models/curto.xml"))  
curto1==curto2

}
\keyword{ math }% at least one, from doc/KEYWORDS
