"fderivM" <-
function(times, X, p)  # state derivative function sent to ODEsolve
{v=rep(0,nrxns)
xp=rep(0,nStates)
St=S0
X[X<0]=0
St[BC==FALSE]=X
if (times<0) m=M[,"control"] else m=M[,patient]

nrules=length(model$rules) 
if (nrules>0) 
    for (j in 1:nrules)
        St[model$rules[[j]]$output]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
for (j in 1:nrxns)
	if (model$rxns[[j]]$rever==FALSE)
	     v[j]=m[rIDs[j]]*model$rxns[[j]]$law(St[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
xp=incid%*%v
names(xp)<-names(y0)
names(v)<-rIDs
aux=c(v,St[BC==T])
list(xp,aux)}    # ******************  END fderiv function definition

