\name{getModelInfo}
\alias{getModelInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get info from an SBML model}
\description{
  This function extracts basic information from a model of class SBML. 
}
\usage{
getModelInfo(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ A model object of class SBML from which information is to be extracted.}
 }
\details{
     
}
\value{ A list containing the following named elements
 \item{BC}{ A logical vector indicating which species are not state variables, 
  i.e. which species are boundary conditions or auxillary variables.}
  \item{y0}{ The initial state (boundary conditions excluded!). }
  \item{nStates}{ The length of the state vector, i.e. the number of system states.}
 \item{S0}{ The full set of species initial values.}
  \item{nReactions}{ The number of reactions.}
  \item{nSpecies}{ The number of species, including 
  states, boundary conditions and possibly auxillary variables such as the total concentration of 
  dihydofolate reductase in the morrison.r  model. }
  \item{incid}{
  The incidence/stoichiometry matrix. This usually contains ones and minus ones corresponding to fluxes either synthesizing or degrading (respectively) 
  a state variable
  chemical species. This matrix multiplied by the flux vector on its right yields 
  the corresponding concentration state variable time derivatives. }
  
  
}
\note{The list output can be attached to immediately define these model variables of interest. }

\references{
For the example given below: Curto R, Voit EO, Sorribas A, Cascante M: Validation and 
steady-state analysis of a power-law model of purine metabolism in man. 
Biochem J. 1997, 324 ( Pt 3):761-775.}


\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}

\seealso{ \code{\link{summary.SBML}}  }
\examples{
library(SBMLR)  
curto=readSBMLR(file.path(.path.package("SBMLR"), "models/curto.r"))  
mi=getModelInfo(curto);mi
attach(mi)
y0
nReactions

# Metabolic Control Analysis of Curto's purine metabolism model
library(SBMLR)  
curto=readSBMLR(file.path(.path.package("SBMLR"), "models/curto.r"))  

mi=getModelInfo(curto)
attach(mi)
N=incid
S0bk=S0
p0=S0[BC==TRUE]

#qr(N)
# full rank => Nr=N and L=I
J=NULL;
for (j in 1:nReactions) 
J[j]=curto$reactions[[j]]$law(S0[c(curto$reactions[[j]]$reactants,curto$reactions[[j]]$modifiers)],curto$reactions[[j]]$parameters)
names(J)<-rIDs
J
DJ=diag(J)

epsS=matrix(rep(0,nReactions*nStates),nrow=nReactions)
rownames(epsS)<-rIDs
colnames(epsS)<-names(y0)
epsS

for (k in 1:nStates) 
#for (k in 1:1) 
{S0=S0bk
y1=y0
y1[k]=1.01*y0[k]
S0=c(y1,p0)
for (j in 1:nReactions) 
epsS[j,k]=(curto$reactions[[j]]$law(S0[c(curto$reactions[[j]]$reactants,curto$reactions[[j]]$modifiers)],curto$reactions[[j]]$parameters)-J[j])/(.01*J[j])
}
epsS
DS=diag(y0)
dSdVm=-solve(N\%*\%DJ\%*\%epsS)\%*\%N\%*\%DJ
colnames(dSdVm)<-rIDs
dSdVm #  compare this output to Figure 2 of Curto et al 1997
apply(dSdVm,1,sum)
dJdVm=epsS\%*\%dSdVm+diag(rep(1,length(rIDs)));dJdVm # this compares well with Figure 3 of Curto et al 1997
apply(dJdVm,1,sum)





}
\keyword{ math }% at least one, from doc/KEYWORDS
