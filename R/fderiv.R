"fderiv" <-
function(t, X, p)  # state derivative function sent to ODEsolve
{v=rep(0,nrxns)
xp=rep(0,nStates)
#print(p)
St=S0
X[X<0]=0
St[BC==FALSE]=X
nrules=length(model$rules) 
if (nrules>0) 
 for (j in 1:nrules)
    St[model$rules[[j]]$output]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
if (p["mod"]==1) if (t<0) m=M[,"control"] else m=M[,patient]
for (j in 1:nrxns)
  if (model$rxns[[j]]$rever==FALSE)
{
if (p["mod"]==0)   v[j]=model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
if (p["mod"]==1)   v[j]=m[rIDs[j]]*model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
if (p["mod"]==2)   v[j]=Mt[[rIDs[j]]](t)*model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
}
xp=incid%*%v
names(xp)<-names(y0)
names(v)<-rIDs
aux=c(v,St[BC==TRUE])
list(xp,aux)}    # ******************  END fderiv function definition
