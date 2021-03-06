\name{sim}
\alias{sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a model of S3 class SBMLR}
\description{
  This function simulates a model given report times and optional modulators. 
  It uses \code{lsoda} of the deSolve package.
}
\usage{
sim(model, times, modulator=NULL,X0=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ The S3 model object to be simulated. Initial conditions are passed through this object.}
  \item{times}{ The sequence of time points to be sampled and provided as rows of the output matrix.}
  \item{modulator}{ Null if there are no modulators (default), 
  a vector of numbers if there are steady state Vmax modulators,  
  and a list of interpolating functions if there are time course Vmax modulators. }
  \item{X0}{ Override model initial conditions in simulations, particularly piece-wise perturbation simulations.} 
  \item{...}{ To pass extra args such as event data frames to deSolve.}
}
\details{
  This is a wrapper for ode. }
\value{ The data frame output that comes out of \code{ode}.}
\references{For the folate cycle example given below: Morrison PF, Allegra CJ: Folate cycle kinetics in human 
  breast cancer cells. JBiolChem 1989, 264(18):10552-10566.}
\author{ Tom Radivoyevitch}
\note{ Rules are implemented through time varying boundary conditions updated at each time point 
as a side effect within the (now internal) function \code{fderiv}. }

\examples{
##---- The following perturbs PRPP from 5 to 50 uM in Curto et al.'s model.
library(SBMLR) 
curto=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=sim(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))

# which should be the same plots as 
curto=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=sim(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))


##---- The following generates Morrison's folate system response to 1uM MTX 
morr=readSBMLR(file.path(system.file(package="SBMLR"), "models/morrison.r"))  
out1=sim(morr,seq(-20,0,1))
morr$species$EMTX$ic=1
out2=sim(morr,0:30)
outs=data.frame(rbind(out1,out2))
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
morr$species$EMTX$ic=0

## Note: This does not work, since EMTX is not a state variable, it is a bc
##(dEMTX1 <- data.frame(var = "EMTX", time = 0, value = 1,method = "add"))
##(out=simulate(morr,times=seq(-20,30,1),events = list(data = dEMTX1) ) )


}
\keyword{ math }% at least one, from doc/KEYWORDS
