library(SBMLR)
setwd(file.path(.path.package("SBMLR"), "BMCcancer04")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed

source(file.path(.path.package("SBMLR"), "models/MorrisonAllegra.r"))  


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
par(mfrow=c(3,3))
#plot(time,MTX1+MTX2+MTX3+MTX4+MTX5+MTX1b+MTX2b+MTX3b+MTX4b+MTX5b,type="l",xlab="Hours",ylab="MTX")
plot(time,EMTX,type="l",xlab="Hours",ylab="EMTX")
plot(time,FH2b+FH2f,type="l",xlab="Hours",ylab="DHF")
#plot(time,FH2f,type="l",xlab="Hours")
plot(time,CH2FH4,type="l",xlab="Hours",ylab="CH2THF")
#plot(time,DHFRf,type="l",xlab="Hours")
#plot(time,DHFRtot,type="l",xlab="Hours")
plot(time,CHOFH4,type="l",xlab="Hours",ylab="CHOTHF")
plot(time,FH4,type="l",xlab="Hours",ylab="THF")
plot(time,CH3FH4,type="l",xlab="Hours",ylab="CH3THF")
#plot(time,AICARsyn,type="l",xlab="Hours")
plot(time,MTHFR,type="l",xlab="Hours",ylab="MTHFR")
plot(time,GARFT,type="l",xlab="Hours",ylab="GARFT")
plot(time,TYMS,type="l",xlab="Hours",ylab="TYMS")
#plot(time,DHFReductase,type="l",xlab="Hours")
par(mfrow=c(1,1))
detach(outs)
S0["EMTX"]=0

dev.copy(pdf,file="fig2runMorr.pdf", width = 7, height = 7)
dev.off()


