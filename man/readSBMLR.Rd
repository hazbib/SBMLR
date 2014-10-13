\name{readSBMLR}
\alias{readSBMLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert an SBMLR file into an R model object of class SBMLR}
\description{
  This function converts an SBMLR model file into a corresponding SBMLR model object.
  This is more than a source-ing: the file is simpler than the object since things are generated,
  such as, rate law and rule R expressions and functions, and mathML. 
}
\usage{
readSBMLR(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{  An SBMLR model definition file. }
}
\details{
A limited subset of SBML level 2 models is currently supported, e.g. events and 
function definitions are not covered.
}
\value{A corresponding \code{SBMLR} model object. }
\author{ Tom Radivoyevitch  }
\note{This function replaces the use of \code{source} in older versions of SBMLR. 
It converts rate law and rule strings to R functions and expressions and to MathML. }

\seealso{ \code{\link{readSBML}} }
\examples{
library(SBMLR)
curtoX=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
curtoR=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
curtoX==curtoR
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line
