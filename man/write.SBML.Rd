\name{write.SBML}
\alias{write.SBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert SBMLR model to SBML}
\description{
  This function converts an SBMLR model definition file into an SBML level 2 file.
}
\usage{
write.SBML(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ This functions takes an SBMLR model in \code{filename}.r and maps it to an SBML 
  model in \code{filename}.xml }
}
\details{
  The XML output file is SBML level 2.  
}
\value{
  No value returned. 
}
\references{
Radivoyevitch, T. A two-way interface between limited Systems Biology Markup Language and R. BMC Bioinformatics 5, 190 (2004).
See http://epbi-radivot.cwru.edu/}
\author{ Tomas Radivoyevitch  (radivot@hal.cwru.edu)}
\note{For additional details see the reference below, also available in the doc directory. }

\section{Warning }{SBML events and 
function definitions are NOT implemented.}

\seealso{ \code{\link{read.SBML}} }
\examples{
library(odesolve)
library(SBMLR)
write.SBML(file.path(.path.package("SBMLR"), "models/MorrisonAllegra")) # writes model to SBML (level 2) file MorrisonAllegra.xml
read.SBML(file.path(.path.package("SBMLR"), "models/MorrisonAllegra")) # converts SBML file MorrisonAllegra.xml into SBMLR model file MorrisonAllegraX.r 
source(file.path(.path.package("SBMLR"), "models/MorrisonAllegraX.r"))  

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
out1=lsoda(y=y0,times=seq(-20,0,1),fderiv,  parms=c(mod=0),  rtol=1e-4, atol= my.atol)
S0["EMTX"]=1
ny0=out1[nrow(out1),2:(nStates+1)]
out2=lsoda(y=ny0,times=seq(0,finalT,1),fderiv,  parms=c(mod=0),  rtol=1e-4, atol= my.atol)
outs=data.frame(rbind(out1,out2));#outs
detach(model$parameters)  # this makes Keq unavailable 
attach(outs)
par(mfrow=c(3,4))
plot(time,FH2b,type="l",xlab="Hours")
plot(time,FH2f,type="l",xlab="Hours")
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
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line
