#"==.SBML"<-function(model1,model2)
"equateModels"<-function(model1,model2)
{
mi=summary(model1)
DFs1=data.frame(initialConcentrations=mi$S0,boundaryConditions=mi$BC,speciesIDs=I(mi$sIDs))
DFr1=data.frame(Laws=I(mi$rLaws),reactionIDs=I(mi$rIDs),initialFluxes=mi$V0)
#print(DFr)

mi=summary(model2)
DFs2=data.frame(initialConcentrations=mi$S0,boundaryConditions=mi$BC,speciesIDs=I(mi$sIDs))
DFr2=data.frame(Laws=I(mi$rLaws),reactionIDs=I(mi$rIDs),initialFluxes=mi$V0)

species=DFs1==DFs2
reactions=DFr1==DFr2
rownames(species)<-mi$sIDs
rownames(reactions)<-mi$rIDs

return(list(species=species, reactions=reactions))
}



