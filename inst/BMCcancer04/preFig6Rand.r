# This is Morrison's folate model driven by random data
library(Biobase)
library(odesolve)
library(annotate)
library(hgu95av2)
library(SBMLR)  
setwd(file.path(.path.package("SBMLR"), "BMCcancer04")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed
source(file.path(.path.package("SBMLR"), "models/MorrisonAllegra.r"))  
morrsym=c('MTHFD1','GART','ATIC','TYMS','DHFR')
morrsym=c('SHMT1','MTHFR','MTR','MTHFD1','GART','ATIC','TYMS','DHFR')
key=c(GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")
key=c(MTHFR="MTHFR",MTR="MTR",SHMT="SHMT1",SHMTr="SHMT1",GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")

npats=1000
aa=matrix(rnorm(npats*length(morrsym),mean=1,sd=.3),ncol=npats)
aa=cbind(aa,control=rep(1,length(morrsym)) )
rownames(aa)=morrsym
colnames(aa)<-c(paste("r",1:npats,sep=""),"control")


nrxns=length(model$rxns);nspcs=length(model$species);   # number of reactions and species 
S0=NULL;BC=NULL;rIDs=NULL  # initialize before assignments
for (j in 1:nrxns) rIDs[j]<-model$rxns[[j]]$id
rIDs

M=matrix(rep(1,dim(aa)[2]*length(rIDs)),nrow=length(rIDs))
rownames(M)<-rIDs
colnames(M)<-colnames(aa)
M[names(key),]
tmp=as.matrix(aa[key,])
rownames(tmp)<-names(key)
M[names(key),]=tmp
M

for (i in 1:nspcs){BC[i]=model$species[[i]]$bc; S0[i]=model$species[[i]]$ic}
names(S0)<-names(model$species) 
y0=S0[BC==F]
nStates=length(y0)
my.atol <- rep(1e-4,nStates)
patient=npats+1   # plot out dynamics for the "control" patient
incid=getIncidenceMatrix(model,BC,y0,nStates,nrxns,nspcs)
nrules=length(model$rules) # this and much above goes into fderiv implicitly by globals
# NOTE: model,incid, nStates, nrxns, rIDs, S0 and BC are all passed globally to fderiv 


attach(model$parameters)  # this makes Keq globally available 
finalT=30
out1=lsoda(y=y0,times=seq(-20,0,1),fderiv,  parms=c(mod=1),  rtol=1e-4, atol= my.atol)
S0["EMTX"]=1
ny0=out1[nrow(out1),2:(nStates+1)]
out2=lsoda(y=ny0,times=seq(0,finalT,1),fderiv,  parms=c(mod=1),  rtol=1e-4, atol= my.atol)
outs=data.frame(rbind(out1,out2));#outs
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

nFluxes=length(rIDs)
# now make the big flux matrix. This takes time to run!!!!!
flux=matrix(rep(0,nFluxes*(npats+1)),ncol=nFluxes,nrow=(npats+1))
conc=matrix(rep(0,nStates*(npats+1)),ncol=nStates,nrow=(npats+1))
rownames(flux)<-c(paste("r",1:npats,sep=""),"control")
colnames(flux)<-rIDs
rownames(conc)<-c(paste("r",1:npats,sep=""),"control")
colnames(conc)<-names(y0)
flux
conc

for (patient in 1:(npats+1))
{
print(patient)
out1=lsoda(y=y0,times=seq(-20,0,1),fderiv,  parms=c(mod=1),  rtol=1e-4, atol= my.atol)
ny0=out1[nrow(out1),2:(nStates+1)]
out2=lsoda(y=ny0,times=seq(0,finalT,1),fderiv,  parms=c(mod=1),  rtol=1e-4, atol= my.atol)
outs=data.frame(rbind(out1,out2));#outs
conc[patient,]=as.numeric(outs[dim(outs)[1],2:(nStates+1)])
flux[patient,]=as.numeric(outs[dim(outs)[1],(nStates+2):(nStates+nFluxes+1)])

par(mfrow=c(3,1))
plot(FH2f~time,data=outs)
title(main=paste("patient ",patient))
plot(FH4~time,data=outs)
plot(CH2FH4~time,data=outs)
par(mfrow=c(1,1))
}

flux=data.frame(flux)
conc=data.frame(conc)
save(flux,conc,file="FmorrRand.Rdata")# save flux array since it takes much time to recompute
#  END big computation loop

detach(model$parameters)  # this makes Keq unavailable 
#  Now do plotting and stats for the predicted fluxes
# load("FmorrRossBT.Rdata") # uncomment this if you saved the flux array 4 lines up 
attach(flux)
attach(conc)
plot(TYMS,MTHFD/2,pch=1,xlim=range(TYMS),ylim=range(MTHFD/2),xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)")
