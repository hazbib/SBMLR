\name{fderivM}
\alias{fderivM}
\title{ fderivM}
\description{
  This function generates the state derivate vector which is passed to lsoda of odesolve. M is passed globally 
  as a modulator of the Vmax values, implementing a step function over time perturbation. 
}
\usage{
fderivM(times, X, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{ \code{times} is the current time}
  \item{X}{ \code{X} is the current state}
  \item{p}{ \code{p} is a parameter vector}
}
\details{
  Since the parameter vector cannot be a list, much is passed globally
}
\value{
  output is a list which goes into lsoda}
\references{ See http://www.sbml.org}
\author{ Tomas Radivoyevitch}
\note{ none }

\seealso{ See Also as lsoda }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(times, X, p)  # state derivative function sent to ODEsolve
{v=rep(0,nrxns)
xp=rep(0,nStates)
St=S0
X[X<0]=0
St[BC==F]=X
if (times<0) m=M[,"control"] else m=M[,patient]

nrules=length(model$rules) 
if (nrules>0) 
    for (j in 1:nrules)
        St[model$rules[[j]]$output]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
for (j in 1:nrxns)
	if (model$rxns[[j]]$rever==F)
	     v[j]=m[rIDs[j]]*model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
xp=incid%*%v
names(xp)<-names(y0)
names(v)<-rIDs
aux=c(v,St[BC==T])
list(xp,aux)}    # ******************  END fderiv function definition
}
\keyword{ math }% at least one, from doc/KEYWORDS
