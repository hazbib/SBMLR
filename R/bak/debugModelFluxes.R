"debugModelFluxes" <-
function(model,V0ref,nrxns)
{# this function compares computed and expected/reference fluxes 
V0s=NULL;
for (j in 1:nrxns) 
if (model$rxns[[j]]$rever==F)
V0s[j]=model$rxns[[j]]$law(S0[c(model$rxns[[j]]$reacts,model$rxns[[j]]$mods)],model$rxns[[j]]$params)
n=min(length(V0ref),length(V0s))
dataf<-data.frame(V0ref=V0ref[1:n], computed=V0s[1:n] )
rownames(dataf)<-rIDs
print(dataf)
} # end debug function

