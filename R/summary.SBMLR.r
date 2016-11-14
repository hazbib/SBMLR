"summary.SBMLR"<-function(object, ...)
{
	model=object
	sIDs=names(model$species)
	rIDs=names(model$reactions)
	ruleIDs=names(model$rules)
	nReactions=length(model$reactions);nSpecies=length(model$species);nRules=length(model$rules) 
	
#	##### by Vishak Venkateswaran (VV)
#	pIDs=names(model$ParametersList)
#	nParams = length(model$ParametersList)
#	#Parameters
#	P0 = NULL; VP = NULL
#	if(nParams > 0)
#		for (i in 1:nParams){
#			P0[i]=model$ParametersList[[i]]$id; 
#			VP[i]=model$ParametersList[[i]]$value
#		}
#	names(VP)<-pIDs
#	################  
# Species
	S0=NULL;
	BC=NULL # initialize 
	for (i in 1:nSpecies){
		if ( "bc" %in% model$species) {BC[i]=model$species[[i]]$bc}; 
#	  BC[i]=model$species[[i]]$bc;
		S0[i]=model$species[[i]]$ic
	}
	names(S0)<-sIDs 
#	names(BC)<-sIDs 
	  if ( "bc" %in% model$species) {names(BC)<-sIDs}
#	y0=S0[BC==FALSE]
	y0=S0
	nStates=length(y0)

# these 4 lines replace the block commented after it  
  rLaws=sapply(model$reactions,function(x) x$strLaw)
  globals=model$globalParameters
	attach(globals)  # e.g. for global coordination of k5 in SOD2012
#   V0=with(model$globalParameters, # not sure why this didn't work in SOD2012
# 			sapply(model$reactions,function(x) x$law(S0[c(x$reactants,x$modifiers)],x$parameters)))
	V0=sapply(model$reactions,function(x) x$law(S0[c(x$reactants,x$modifiers)],x$parameters))
	names(V0)<-rIDs
	detach(globals)  
	
# 	attach(model$globalParameters)  # e.g. for Vmax global coordination in RNR
# 	
# 	# Reactions
# 	rLaws=NULL;V0=NULL # initialize
# 	for (j in 1:nReactions) {
# 	  rLaws[j]<-model$reactions[[j]]$strLaw  	#this gives you null which is wrong
# 	  
# 	  V0[j]=model$reactions[[j]]$law(S0[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers,model$reactions[[j]]$products)], model$reactions[[j]]$parameters)
# 	}
# 	names(rLaws)<-rIDs
# 	names(V0)<-rIDs
# 	detach(model$globalParameters)  
  
  
# Incidence Matrix
#	incid=matrix(rep(0,nStates*nReactions),nrow=nStates)
#	indx=(1:nSpecies)[BC==FALSE]
#	for (i in 1:nStates)
#		for (j in 1:nReactions)
#		{if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$products)) 
#				incid[i,j] = summary(factor(model$reactions[[j]]$products))[[model$species[[indx[i]]]$id]]
#			if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$reactants)) 
#				incid[i,j] = incid[i,j]-summary(factor(model$reactions[[j]]$reactants))[[model$species[[indx[i]]]$id]]  }     
#	rownames(incid)<-names(y0)
	
	incid=matrix(rep(0,nStates*nReactions),nrow=nStates)
	indx=(1:nSpecies)[BC==FALSE]
	for (i in 1:nStates)
	  for (j in 1:nReactions)
	  {if ( is.element(model$species[[i]]$id, model$reactions[[j]]$products)) 
	    incid[i,j] = summary(factor(model$reactions[[j]]$products))[[model$species[[i]]$id]]
	  if ( is.element(model$species[[i]]$id, model$reactions[[j]]$reactants)) 
	    incid[i,j] = incid[i,j]-summary(factor(model$reactions[[j]]$reactants))[[model$species[[i]]$id]]  }     
	rownames(incid)<-names(y0)
	
# return a list of model information
	options(stringsAsFactors=FALSE)
	DFs=data.frame(index=1:nSpecies,initialConcentrations=S0,boundaryConditions=BC);row.names(DFs)<-sIDs
	DFr=data.frame(index=1:nReactions,Laws=rLaws,initialFluxes=V0);   row.names(DFr)<-rIDs
	list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,nStates=nStates,
			y0=y0,nReactions=nReactions,rIDs=rIDs,rLaws=rLaws, V0=V0,globalVec=unlist(globals), # P0, VP,
			incid=incid,nRules=nRules,ruleIDs=ruleIDs,species=DFs,reactions=DFr) 
}



