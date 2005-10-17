"simulate" <-
function(model,times, modulator=NULL,X0=NULL, ...)  # this is a wrapper for lsoda
{



 mi=summary(model)
# print(mi);next
 attach(mi)
 my.atol <- rep(1e-4,nStates)
# print(nRules)
 attach(model$globalParameters)  # in Morrison this makes Keq globally available 
 mod=0
 if (class(modulator)=="numeric") mod=1
 if (class(modulator)=="list") mod=2

#print(mod)
#print(modulator)

fderiv<-function(t, X, p)  # state derivative function sent to ODEsolve
{v=rep(0,nReactions)
xp=rep(0,nStates)
#print(p)
St=S0
X[X<0]=0
St[BC==FALSE]=X
#print(model$rule[[1]]$inputs)
#print(St[model$rule[[1]]$inputs])
if (nRules>0) 
 for (j in 1:nRules)
    St[model$rules[[j]]$idOutput]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
if (p["mod"]==1) {if (t<0) m=rep(1,length(modulator)) else m=modulator
names(m)<-rIDs}

for (j in 1:nReactions)
  if (model$reactions[[j]]$reversible==FALSE)
{
#print(t)
if (p["mod"]==0)   v[j]=model$reactions[[j]]$law(St[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers)],model$reactions[[j]]$parameters)
if (p["mod"]==1)   v[j]=m[rIDs[j]]*model$reactions[[j]]$law(St[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers)],model$reactions[[j]]$parameters)
if (p["mod"]==2)   v[j]=modulator[[rIDs[j]]](t)*model$reactions[[j]]$law(St[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers)],model$reactions[[j]]$parameters)
}
xp=incid%*%v
names(xp)<-names(y0)
names(v)<-rIDs
aux=c(v,St[BC==TRUE])
list(xp,aux)}    # ******************  END fderiv function definition



if (!is.null(X0)) out=lsoda(y=X0,times=times,fderiv,  parms=c(mod=mod),  rtol=1e-4, atol= my.atol) else
                  out=lsoda(y=y0,times=times,fderiv,  parms=c(mod=mod),  rtol=1e-4, atol= my.atol) 

 detach(mi)
 detach(model$globalParameters)
 out
} 
