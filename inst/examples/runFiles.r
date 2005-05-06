setwd("C:/cwru/active/XMLtest")
library(XML)
library(odesolve)
dir.create(file.path(tempdir(),"SBMLR"))
source(file.path(.path.package("SBMLR"), "models/Curto.r"))
fileName <- system.file("models", "curto.r", package="SBMLR")
write.SBML(file.path(.path.package("SBMLR"), "models/Curto")) 
fileName

source("sbmlHandler.r")
source("mathml2R.r")
source("R2mathml.r")
source("makeLaw.r")
source("SBMLR2model.r")
source("sbml2model.r")
source("model2SBMLR.r")
source("summary.SBML.r");
source("==.SBML.r")
source("modelInfo.r")
source("getIncidenceMatrix.r")
source("simulate.r")
source("SBMLRf2model.r")


source("SBMLR2model.r");morr1<-SBMLR2model("morrison")
source("model2SBMLR.r");model2SBMLR(morr1,"morr1")
source("model2sbml.noDOM.r");model2sbml.noDOM(morr1,"morrNoDom")


#morr0<-SBMLR2model("morrisonAllegra")

source("sbmlHandler.r");source("sbml2model.r");morrND=sbml2model("morrNoDom")
morrND$globalParameters


#morr1==morrND


model2SBMLR(morrND,"morrND")

source("SBMLR2model.r");morrNDM<-SBMLR2model("morrND")

morr=morrND
source("simulate.r");out1=simulate(morr,seq(-20,0,1))
morr$species$EMTX$ic=1
out2=simulate(morr,0:30)
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

morrND$globalParameters


source("model2SBMLR.r");model2SBMLR(morrOld,"morrOld")




source("sbml2SBMLR.r")



help.search(keyword="character") 

source("model2sbml.r");model2sbml(morr,"morr")







#source("modelInfo.r");mi=modelInfo(morr)
#attach(mi)
#rIDs
#y0
#detach(mi)
#





