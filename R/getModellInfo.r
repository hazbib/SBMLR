"getModelInfo"<-function(model)
{
sIDs=names(model$species)
rIDs=names(model$reactions)
ruleIDs=names(model$rules)
nReactions=length(model$reactions);nSpecies=length(model$species);nRules=length(model$rules) 

# Species
S0=NULL;BC=NULL # initialize 
for (i in 1:nSpecies){
  BC[i]=model$species[[i]]$bc; 
  S0[i]=model$species[[i]]$ic
  }
names(S0)<-sIDs 
names(BC)<-sIDs 
y0=S0[BC==FALSE]
nStates=length(y0)

# Reactions
rLaws=NULL;V0=NULL # initialize
for (j in 1:nReactions) {
   rLaws[j]<-model$reactions[[j]]$strLaw
   V0[j]=model$reactions[[j]]$law(S0[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers)], model$reactions[[j]]$parameters)
   }
names(rLaws)<-rIDs
names(V0)<-rIDs

# Incidence Matrix
incid=matrix(rep(0,nStates*nReactions),nrow=nStates)
indx=(1:nSpecies)[BC==FALSE]
for (i in 1:nStates)
  for (j in 1:nReactions)
	{if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$products)) 
	     incid[i,j] = summary(factor(model$reactions[[j]]$products))[[model$species[[indx[i]]]$id]]
	if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$reactants)) 
	     incid[i,j] = incid[i,j]-summary(factor(model$reactions[[j]]$reactants))[[model$species[[indx[i]]]$id]]  }     
rownames(incid)<-names(y0)

# return a list of model information
list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,nStates=nStates,y0=y0,nReactions=nReactions,rIDs=rIDs,rLaws=rLaws, V0=V0, incid=incid,nRules=nRules,ruleIDs=ruleIDs)
}



