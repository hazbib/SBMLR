# MCA perturbation of Morrison's folate model 
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
qr(N)
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
dSdVm=-solve(N%*%DJ%*%epsS)%*%N%*%DJ
colnames(dSdVm)<-rIDs
dSdVm #  compare this output to Figure 2 of Curto et al 1997
apply(dSdVm,1,sum)
dJdVm=epsS%*%dSdVm+diag(rep(1,length(rIDs)));dJdVm # this compares well with Figure 3 of Curto et al 1997
apply(dJdVm,1,sum)

