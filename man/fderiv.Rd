\name{fderiv}
\alias{fderiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ fderiv}
\description{
  This function generates the state derivate vector which is passed to lsoda of the odesolve package.
}
\usage{
fderiv(times, X, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{ \code{times} is the current time (not used)}
  \item{X}{ \code{X} is the current state}
  \item{p}{ \code{p} is a parameter vector in principle but it is not used, see details.}
}
\details{
  Since the parameter vector cannot be a list, much is passed globally. For example, the incidence
  matrix incid (which must be named "incid") must be defined globally before fderiv is passed to the lsoda call.
}
\value{
  The time derivative of the concentration state variable vector. }
\references{For the folate cycle example given below: Morrison PF, Allegra CJ: Folate cycle kinetics in human 
  breast cancer cells. JBiolChem 1989, 264(18):10552-10566.}
\author{ Tomas Radivoyevitch}
\note{ Rules are implemented through time varying boundary conditions updated at each time point 
as a side effect of fderiv. }

\seealso{ \code{\link{getIncidenceMatrix}} }
\examples{
##---- The following example uses fderiv to generate Morrison's folate system response to 1uM MTX 
library(SBMLR) 
setwd(file.path(.path.package("SBMLR"), "demo"))
library(odesolve)
source("MorrisonAllegra.r")  

nrxns=length(model$rxns);nspcs=length(model$species);   # number of reactions and species 
S0=NULL;BC=NULL;rIDs=NULL  # initialize before assignments
for (j in 1:nrxns) rIDs[j]<-model$rxns[[j]]$id
for (i in 1:nspcs){BC[i]=model$species[[i]]$bc; S0[i]=model$species[[i]]$ic}
names(S0)<-names(model$species) 
names(model$rxns)<-rIDs

sIDs<-names(model$species) 
names(S0)<-sIDs;S0
y0=S0[BC==FALSE]
nStates=length(y0)
my.atol <- rep(1e-4,nStates)


incid=getIncidenceMatrix(model,BC,y0,nStates,nrxns,nspcs)
Inc<-data.frame(incid);Inc # Inc is just for inspection
names(Inc)<-rIDs
Inc

nrules=length(model$rules) # this and much above goes into fderiv implicitly by globals

# *****NOTE!!!!***** incid, nrules, nStates, nrxns, BC, and tvBC all need to be defined globally before the lsoda call

attach(model$parameters)  # this makes Keq globally available 
finalT=30
out1=lsoda(y=y0,times=seq(-20,0,1),fderiv,  parms=c(test1=1),  rtol=1e-4, atol= my.atol)
S0["EMTX"]=1
ny0=out1[nrow(out1),2:(nStates+1)]
out2=lsoda(y=ny0,times=seq(0,finalT,1),fderiv,  parms=c(test1=1),  rtol=1e-4, atol= my.atol)
outs=data.frame(rbind(out1,out2));#outs
detach(model$parameters)  # this makes Keq unavailable 
attach(outs)
par(mfrow=c(3,4))
plot(time,FH2b,type="l",xlab="Hours")
plot(time,FH2f,type="l",xlab="Hours",xlim=c(-1,finalT))
plot(time,DHFRf,type="l",xlab="Hours")
plot(time,DHFRtot,type="l",xlab="Hours")
plot(time,CHOFH4,type="l",xlab="Hours")
plot(time,FH4,type="l",xlab="Hours")
plot(time,CH2FH4,type="l",xlab="Hours")
plot(time,CH3FH4,type="l",xlab="Hours")
plot(time,AICARsyn,type="l",xlab="Hours")
plot(time,MTR,type="l",xlab="Hours")
plot(time,TYMS,type="l",xlab="Hours")
#plot(time,EMTX,type="l",xlab="Hours")
plot(time,DHFReductase,type="l",xlab="Hours")
par(mfrow=c(1,1))
detach(outs)
S0["EMTX"]=0

}
\keyword{ math }% at least one, from doc/KEYWORDS
