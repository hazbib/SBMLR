\name{getIncidenceMatrix}
\alias{getIncidenceMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Incidence Matrix from SBMLR Model}
\description{
  This function produces the incidence/stoichiometric matrix of an SBMLR model. 
}
\usage{
getIncidenceMatrix(model, BC, y0, nStates, nrxns, nspcs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ \code{model} is the SBMLR model for which the incidence matrix is being computed.}
  \item{BC}{ \code{BC} is a logical vector indicating which species are not state variables, 
  i.e. which species are boundary conditions or auxillary variables.}
  \item{y0}{ \code{y0} is the initial state. This is passed just to get names for the rows of the output matrix.}
  \item{nStates}{ \code{nStates} is the length of the state vector, i.e. the number of rows in the output.
  matrix [i.e. the length of y0]}
  \item{nrxns}{ \code{nrxns}  is the number of reactions and thus the number of output matrix columns}
  \item{nspcs}{ \code{nspcs} is the size of the species vector. The species vector includes 
  states, boundary conditions and possibly auxillary variables such as the total concentration of 
  dihydofolate reductase in the MorrisonAllegra.r SBMLR model given in the demo directory, see the  
  second rule defined therein.}
}
\details{
  Compared to fderiv, variables are passed to this function explicitly because   there are no
  constraints on function arguments. Although arguments could have been inferred from 
  \code{model} as in the code that precedes the function call in the example given below, users are 
  encouraged to explore these arguments in .GlobalEnv to check their model, and in this case,
  the arguments would 
  have already been computed anyway, so they may as well be passed.   
}
\value{
Matrix of ones and minus ones corresponding to fluxes either synthesizing or degrading (respectively) 
a state variable
chemical species. This matrix multiplied by the flux vector on its right yields 
the corresponding concentration state variable time derivatives. 
}
\references{For the example given below: Curto R, Voit EO, Sorribas A, Cascante M: V
alidation and steady-state analysis of a power-law model of purine metabolism in man. 
Biochem J. 1997, 324 ( Pt 3):761-775.}
\author{ Tomas Radivoyevitch}
\note{For additional details see the manual in the doc directory. }

\seealso{ \code{\link{fderiv}}  }
\examples{
# Metabolic Control Analysis of Curto's purine metabolism model
library(SBMLR)  
setwd(file.path(.path.package("SBMLR"), "demo"))
source("curto.r")  

nrxns=length(model$rxns);nspcs=length(model$species);   # number of reactions and species 
S0=NULL;BC=NULL;rIDs=NULL  # initialize before assignments
for (j in 1:nrxns) rIDs[j]<-model$rxns[[j]]$id
for (i in 1:nspcs){BC[i]=model$species[[i]]$bc; S0[i]=model$species[[i]]$ic}
names(S0)<-names(model$species) 
S0bk=S0
y0=S0[BC==FALSE]
p0=S0[BC==TRUE]
nStates=length(y0)
N=getIncidenceMatrix(model,BC,y0,nStates,nrxns,nspcs)
#qr(N)
# full rank => Nr=N and L=I
J=NULL;
for (j in 1:nrxns) 
J[j]=model$rxns[[j]]$law(S0[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
names(J)<-rIDs
J
DJ=diag(J)

epsS=matrix(rep(0,nrxns*nStates),nrow=nrxns)
rownames(epsS)<-rIDs
colnames(epsS)<-names(y0)
epsS

for (k in 1:nStates) 
#for (k in 1:1) 
{S0=S0bk
y1=y0
y1[k]=1.01*y0[k]
S0=c(y1,p0)
for (j in 1:nrxns) 
epsS[j,k]=(model$rxns[[j]]$law(S0[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)-J[j])/(.01*J[j])
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
