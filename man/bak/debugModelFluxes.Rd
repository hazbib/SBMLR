\name{debugModelFluxes}
\alias{debugModelFluxes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ debugModelFluxes}
\description{
  This cranks out the flux vector for the current state vector, i.e. helps debug the reaction rate laws}
\usage{
debugModelFluxes(model, V0ref, nrxns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ \code{model} The SBML-R model}
  \item{V0ref}{ \code{V0ref} the reference fluxes that the model should have}
  \item{nrxns}{ \code{nrxns} is the number of reactions in the model}
}
\details{
  no extra details
}
\value{
  flux vector output}
\references{ See http://www.sbml.org}
\author{ Tomas Radivoyevitch}
\note{ none }


\seealso{ See Also as \code{\link{write.SBML}}}
\examples{


}
\keyword{ math }% at least one, from doc/KEYWORDS
