\name{read.SBML}
\alias{read.SBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert SBML model to SBMLR }
\description{
  This function converts an SBML level 2 model in filename.XML into a corresponding SBMLR structure in filenameX.r.
}
\usage{
read.SBML(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{  An SBML level 2 model contained in the file 
  \code{filename}.xml is mapped into an SBMLR model in \code{filename}X.r.}
}
\details{
A limited subset of SBML level 2 models is currently supported, e.g. events and 
function definitions are not covered.
}
\value{
No value. }
\references{
Radivoyevitch, T. A two-way interface between limited Systems Biology Markup Language and R. BMC Bioinformatics 5, 190 (2004).
See http://epbi-radivot.cwru.edu/}
\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}
\note{For additional details see the reference below, also available in the doc directory. }

\seealso{ \code{\link{write.SBML}} }
\examples{
library(SBMLR)
library(odesolve)
setwd(file.path(.path.package("SBMLR"), "models"))
source("Curto.r") # this file defines Curto et al.'s model as an SBML-R model structure 

nrxns=length(model$rxns);nspcs=length(model$species);   # number of reactions and species 
S0=NULL;BC=NULL;rIDs=NULL  # initialize before assignments
for (j in 1:nrxns) rIDs[j]<-model$rxns[[j]]$id
for (i in 1:nspcs){BC[i]=model$species[[i]]$bc; S0[i]=model$species[[i]]$ic}
names(S0)<-names(model$species) 
y0=S0[BC==FALSE]
nStates=length(y0)
my.atol <- rep(1e-4,nStates)
finalT=70
incid=getIncidenceMatrix(model,BC,y0,nStates,nrxns,nspcs)
# NOTE: model,incid, nStates, nrxns, rIDs, S0 and BC are all passed globally to fderiv 
out1=lsoda(y=y0,times=seq(-20,0,1),fderiv,  parms=c(mod=0),  rtol=1e-4, atol= my.atol)
ny0=out1[nrow(out1),2:(nStates+1)]
ny0["PRPP"]=50  # step response to PRPP change from 5 uM to 50 uM 
out2=lsoda(y=ny0,times=seq(0,finalT,1),fderiv,  parms=c(mod=0),  rtol=1e-4, atol= my.atol)
outs=data.frame(rbind(out1,out2));#outs
# the next block plots the dynamic responses in fig. 2 of Curto et al (1998) Math Biosci
attach(outs)
par(mfrow=c(2,1))
plot(time,IMP,type="l")
plot(time,HX,type="l")
par(mfrow=c(1,1))
detach(outs)

write.SBML("Curto") # writes model to SBML (level 2) file Curto.xml
read.SBML("Curto") # converts SBML file Curto.xml into SBMLR model file CurtoX.r 
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line
