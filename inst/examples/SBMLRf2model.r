"SBMLRf2model"<-function(filename)
{# SBMLR file input, SBML model object output
source(paste(filename,".R",sep=""),local=TRUE)

nrxns=length(model$rxns);nspcs=length(model$species);nrules=length(model$rules) 
sIDs=NULL;rIDs=NULL;ruleIDs=NULL  # initialize before assignments
for (j in 1:nrxns) rIDs[j]<-model$rxns[[j]]$id
for (i in 1:nspcs) sIDs[i]<-model$species[[i]]$id
for (i in 1:nrules) ruleIDs[i]<-model$rules[[i]]$id
names(model$species)<-sIDs
names(model$rxns)<-rIDs
names(model$rules)<-ruleIDs
notes=model$notes
nnotes=length(notes)
 con <- xmlOutputDOM()
 con$addTag("notes", close=F)
 con$addTag("body",attrs=c(xmlns="http://www.w3.org/1999/xhtml"), close=F  )
for (i in 1:nnotes)  
{
con$addTag("p",close=F  )
con$addNode(xmlTextNode(notes[i]))
con$closeTag()
}
con$closeTag();con$closeTag()
dom=con$value()[[1]]
model$htmlNotes=dom

#print(model)

for (i in 1:nrules)
{
bod=body(model$rules[[i]]$law)
nbod=length(bod)
model$rules[[i]]$exprLaw=bod[[nbod]]
model$rules[[i]]$strLaw=gsub(" ","",toString(bod[nbod]))
e=model$rules[[i]]$exprLaw
model$rules[[i]]$mathmlLaw=R2MathML(e)$value()[[1]]
}



for (i in 1:nrxns)
{
bod=body(model$rxns[[i]]$law)
nbod=length(bod)
model$rxns[[i]]$exprLaw=bod[[nbod]]
model$rxns[[i]]$strLaw=gsub(" ","",toString(bod[nbod]))
#print("here")
e=model$rxns[[i]]$exprLaw
#print(r)
#print("here")
model$rxns[[i]]$mathmlLaw=R2MathML(e)$value()[[1]]
}

class(model)<-"SBML"
return(model)}



