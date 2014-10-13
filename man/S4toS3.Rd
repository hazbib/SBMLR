\name{S4toS3}
\alias{S4toS3}
\title{ Converts an S4 class SBML object created by rsbml into an S3 object of class SBMLR}
\description{
  This function provides a path from rsbml to SBMLR. The latter, being S3, is less cluttered
  with empty fields/slots than the former. The advantage of the S4 object is that it comes from
  more robust SBML reading: rsbml uses libsbml to parse SBML, SBMLR  uses the R package XML.
}
\usage{
S4toS3(dom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dom}{An S4 DOM object of class SBML produced by rsbml.}
}
\details{
Carried over are compartments, species, global parameters, rules and reactions.
}
\value{A corresponding SBMLR model object, i.e. an S3 list of lists type of object. }
\author{ Tom Radivoyevitch }

\examples{

library(rsbml)
(dom <- rsbml_read(file.path(system.file(package="SBMLR"), "models/sod.xml")))
library(SBMLR)
(mod=S4toS3(dom))
summary(mod)
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line
