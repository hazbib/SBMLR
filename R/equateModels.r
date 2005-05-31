"==.SBML"<-function(model1,model2,...)
{
mi=getModelInfo(model1)
attach(mi)
DFs1=data.frame(initialConcentrations=S0,boundaryConditions=BC,speciesIDs=I(sIDs))
DFr1=data.frame(Laws=I(rLaws),reactionIDs=I(rIDs),initialFluxes=V0)
#print(DFr)
detach(mi)

mi=getModelInfo(model2)
attach(mi)
DFs2=data.frame(initialConcentrations=S0,boundaryConditions=BC,speciesIDs=I(sIDs)) 
DFr2=data.frame(Laws=I(rLaws),reactionIDs=I(rIDs),initialFluxes=V0);

species=DFs1==DFs2
reactions=DFr1==DFr2
rownames(species)<-sIDs
rownames(reactions)<-rIDs
detach(mi)

return(list(species=species, reactions=reactions))
}



