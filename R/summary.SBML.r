"summary.SBML"<-function(object, ...)
{

mi=getModelInfo(object)
print(mi)
attach(mi)
DFs=data.frame(initialConcentrations=S0,boundaryConditions=BC,speciesIDs=sIDs); row.names(DFs)<-1:nSpecies
DFr=data.frame(Laws=rLaws,reactionIDs=rIDs,initialFluxes=V0);   row.names(DFr)<-1:nReactions
detach(mi)
list(species=DFs,reactions=DFr)
}



