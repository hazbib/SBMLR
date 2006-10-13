\name{makeLaw}
\alias{makeLaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate an R function for the reaction rate law}
\description{
  This function creates an R function for a rate law given three arguments, the inputs variable, the parameters, and 
  the rate law as an R expression. 
}
\usage{
makeLaw(r,p,e)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ A vector of the reactant and modulator names.}
  \item{p}{ A named numeric vector of the function's local parameters.}
 \item{e}{ An R expression (i.e. nested calls and tokens) of the reaction rate law. }
 }

\details{
}

\value{ An R function that returns the value of \code{e} given \code{r} and \code{p}, e.g. a rate law.}

\note{This function is also used for rules with p=NULL. It is used by 
\code{readSBML} and \code{readSBMLR}. }

\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}

\examples{
library(SBMLR)  
curto=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
curto$reactions[["ada"]]
e=curto$reactions[["ada"]]$exprLaw;e
r=curto$reactions[["ada"]]$reactants;r
m=curto$reactions[["ada"]]$modifiers;m
r=c(r,m)
p=curto$reactions[["ada"]]$parameters;p
makeLaw(r,p,e)

## compare to 

curto$reactions[["ada"]]$law

## indeed, above is how $law functions are now made, and in part why readSBMLR has replaced source.

}
\keyword{ math }% at least one, from doc/KEYWORDS
