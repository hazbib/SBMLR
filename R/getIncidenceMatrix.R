"getIncidenceMatrix" <-
function(model,BC,y0,nStates,nrxns,nspcs)
{
incid=matrix(rep(0,nStates*nrxns),nrow=nStates)
indx=(1:nspcs)[BC==FALSE]
for (i in 1:nStates)
  for (j in 1:nrxns)
	{if ( is.element(model$species[[indx[i]]]$id, model$rxns[[j]]$prods)) 
	     incid[i,j] = summary(factor(model$rxns[[j]]$prods))[[model$species[[indx[i]]]$id]]
	if ( is.element(model$species[[indx[i]]]$id, model$rxns[[j]]$reacts)) 
	     incid[i,j] = incid[i,j]-summary(factor(model$rxns[[j]]$reacts))[[model$species[[indx[i]]]$id]]  }     
rownames(incid)<-names(y0)      
incid}  # end getIncidenceMatrix function definition

