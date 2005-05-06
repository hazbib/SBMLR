"model2sbml" <-
function(model,filename)
{  # takes model into an SBML file in filename.xml 

sbml=model[["sbml"]]
id=model[["id"]]
notes=model[["notes"]]
htmlNotes=model[["htmlNotes"]]
comps=model[["comps"]]
species=model[["species"]]
parameters=model[["parameters"]]
rules=model[["rules"]]
rxns=model[["rxns"]]
units=model[["units"]]

nnotes=length(notes)
ncomps=length(comps)
nrxns=length(rxns)
nspcs=length(species)
nparameters=length(parameters)
nrules=length(rules)
nunits=length(units)

 con <- xmlOutputDOM()
 con$addTag("sbml",attrs=c(xmlns="http://www.sbml.org/sbml/level2",level="2",version="1"),close=F)
 con$addTag("model",attrs=c(id=id),close=F)
 con$addNode(htmlNotes)

if(ncomps>0){
 con$addTag("listOfCompartments",close=F)
for (i in 1:ncomps) {con$addNode(xmlNode("compartment",attrs=c(id=comps[[i]][["id"]],size=comps[[i]][["size"]])))} 
con$closeTag()
}

if(nspcs>0){
 con$addTag("listOfSpecies",close=F)
for (i in 1:nspcs) {
names(species[[i]])<-c("id","initialConcentration","compartment","boundaryCondition") 
con$addTag("species",attrs=species[[i]])}
con$closeTag()
}
if(nparameters>0){
 con$addTag("listOfParameters",close=F)
for (i in 1:nparameters)  con$addTag("parameter",
                          attrs=c(id=names(parameters[i]), value=parameters[[i]]))
con$closeTag()
}


if(nrules>0){
 con$addTag("listOfRules",close=F)
for (i in 1:nrules)  
{con$addTag("assignmentRule",attrs=c(variable=rules[[i]][["id"]]) ,close=F)
math=xmlNode("math",attrs=c(xmlns="http://www.w3.org/1998/Math/MathML"))
math2=rules[[i]][["mathmlLaw"]]
math3=append.xmlNode(math,math2)
#print(math3)
#print(class(math3))
con$addNode(math3)  # add and close mathml
con$closeTag()  # close rule
} # for loop on reactions
con$closeTag()  # close list of rules (if there are any!)
}  # end if nrules> 0 


if(nrxns>0){
 con$addTag("listOfReactions",close=F)
for (i in 1:nrxns)  
{con$addTag("reaction",attrs=c(id=rxns[[i]][["id"]], reversible=rxns[[i]][["rever"]]) ,close=F)
math=xmlNode("math",attrs=c(xmlns="http://www.w3.org/1998/Math/MathML"))
math2=rxns[[i]][["mathmlLaw"]]
math3=append.xmlNode(math,math2)
#print(math3)
#print(class(math3))

params=rxns[[i]][["params"]]
reacts=rxns[[i]][["reacts"]]
mods=rxns[[i]][["mods"]]
prods=rxns[[i]][["prods"]]

nr=length(reacts)
if (nr>0){ con$addTag("listOfReactants",close=F)
for (k in 1:nr) con$addTag("speciesReference",attrs=c(species=reacts[k]) )
con$closeTag()
}

np=length(prods)
if (np>0){ con$addTag("listOfProducts",close=F)
for (k in 1:np) con$addTag("speciesReference",attrs=c(species=prods[k]) )
con$closeTag()
}

nm=length(mods)
if (nm>0){
con$addTag("listOfModifiers",close=F)
for (k in 1:nm) con$addTag("modifierSpeciesReference",attrs=c(species=mods[k]) )
con$closeTag()
}
con$addTag("kineticLaw",close=F)
con$addNode(math3)  # add and close mathml
npa=length(params)
if (npa>0){
con$addTag("listOfParameters",close=F)
for (k in 1:npa) con$addTag("parameter",attrs=c(id=names(params[k]), value=params[[k]]) )
con$closeTag()  # close params
}
con$closeTag()  # close kinetic law
con$closeTag()  # close reaction
} # for loop on reactions
con$closeTag()  # close list of reactions (if there are any!)
}  # end if rxns> 0 

con$closeTag() # close model
con$closeTag() # close sbml



 dom=con$value()[[1]];dom
 dom$id=dom$id[[1]]
 saveXML(dom,file=paste(filename,".xml",sep=""),prefix ='<?xml version="1.0" encoding="UTF-8"?>\n')
 }  # end write.sbml function defition

