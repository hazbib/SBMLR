\name{fderiv}
\alias{fderiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fderiv(times, X, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{ ~~Describe \code{times} here~~ }
  \item{X}{ ~~Describe \code{X} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
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
nrules=length(model$rules) 
if (nrules>0) 
    for (j in 1:nrules)
        St[model$rules[[j]]$output]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
for (j in 1:nrxns)
	if (model$rxns[[j]]$rever==F)
	     v[j]=model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
xp=incid%*%v
names(xp)<-names(y0)
names(v)<-rIDs
aux=c(v,St[BC==T])
list(xp,aux)}    # ******************  END fderiv function definition
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
