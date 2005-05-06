"saveSBMLR" <-
function(model,filename)
{  # takes model and maps it to a SBMLR file

#               fid <- file(paste(filename,".r",sep=""), "w")  # open the output file connection

fid <- file(filename,"w")  # open the output file connection
cat("#",filename,"\n", file=fid, sep="")  
cat("model=list( \n", file=fid, sep="\n")

# read in the R file 

#sbml=model[["sbml"]]
id=model[["id"]]
notes=model[["notes"]]
#htmlNotes=model[["htmlNotes"]]
compartments=model[["compartments"]]
species=model[["species"]]
globalParameters=model[["globalParameters"]]
rules=model[["rules"]]
reactions=model[["reactions"]]
units=model[["units"]]

nNotes=length(notes)
nCompartments=length(compartments)
nReactions=length(reactions)
nSpecies=length(species)
nGlobalParameters=length(globalParameters)
nRules=length(rules)
nUnits=length(units)

#cat("sbml=c(level=", sbml[2],", version=",sbml[3],"),\n\n", file=fid,sep="" )
cat("id=\"",id,"\",\n\n", file=fid, sep=""  )
#cat("htmlNotes=c(\"",toString(htmlNotes),"\"),\n", file=fid, sep="\n")

if(nNotes>0){
cat("notes=c(", file=fid, sep="\n")
for (i in 1:nNotes) 
cat(sprintf("\"%s\"%s", notes[[i]],ifelse(i==nNotes,"\n),\n",",")), file=fid, sep="\n")
}



if(nCompartments>0){
cat("compartments=list(", file=fid, sep="\n")
for (i in 1:nCompartments) 
cat(sprintf("list(id=\"%s\", size = %g)%s", compartments[[i]][["id"]],compartments[[i]][["size"]],ifelse(i==nCompartments,"\n),\n",",")), 
file=fid, sep="\n")
}

if(nSpecies>0){
cat("species=list(", file=fid, sep="\n")
for (i in 1:nSpecies) 
cat(sprintf(paste("list(id=%-20s,ic=%8g,    compartment=%-10s,bc=%-7s)%s",sep=""),dQuote(species[[i]][["id"]]),species[[i]][["ic"]],
dQuote(species[[i]][["compartment"]]),ifelse(species[[i]][["bc"]],"TRUE","FALSE"),ifelse(i==nSpecies,"\n),\n",",") ), file=fid, sep="\n")
}


if(nGlobalParameters>0){
cat("globalParameters=list(", file=fid, sep="")
for (i in 1:nGlobalParameters) 
{
#print(globalParameters[i])
#cat(names(globalParameters[i]),globalParameters[i],  sep="\n")
#cat(sprintf("%s=%g", names(globalParameters[i]),globalParameters[[i]] ),  sep="\n")
cat(sprintf("%s=%g%s", names(globalParameters[i]),globalParameters[[i]],ifelse(i==nGlobalParameters,"),\n",",") ), file=fid, sep="\n")
}
}

if(nRules>0){
cat("\n\n#BEGIN the Rules Section\nRules=list(", file=fid, sep="\n")
for (i in 1:nRules)
{
idOutput=rules[[i]][["idOutput"]]
inputs=rules[[i]][["inputs"]]
strLaw=rules[[i]][["strLaw"]]
cat(sprintf("list(idOutput=%s,",dQuote(idOutput)),file=fid, sep="\n")
#print("here")
ni=length(inputs)
if (ni>0){cat("inputs=c(", file=fid, sep="")
for (k in 1:ni) cat(sprintf("\"%s\"%s",inputs[k],ifelse(k==ni,"),\n",",")  ), file=fid, sep="")}
cat(sprintf("strLaw = %s",dQuote(strLaw)),file=fid, sep="\n")
cat(sprintf("%s",ifelse(i==nRules,")),\n# END of the Rules Section","),\n")  ), file=fid, sep="\n")

}  # end for loop on rules
} # end if there are rules
#print("here")

if(nReactions>0){
cat("\n\n#BEGIN the Reaction Section \n reactions=list(", file=fid, sep="\n")
for (i in 1:nReactions)
{
cat(sprintf("list( id=\"%s\", reversible=%s,",reactions[[i]][["id"]],ifelse(reactions[[i]][["rever"]],"TRUE","FALSE")),file=fid, sep="\n")
law=reactions[[i]][["law"]]
parameters=reactions[[i]][["parameters"]]
reactants=reactions[[i]][["reactants"]]
modifiers=reactions[[i]][["modifiers"]]
products=reactions[[i]][["products"]]
strLaw=reactions[[i]][["strLaw"]]

nr=length(reactants)
if (nr>0){cat("reactants=c(", file=fid, sep="")
for (k in 1:nr) cat(sprintf("\"%s\"%s",reactants[k], ifelse(k==nr,"),\n",",")  ), file=fid, sep="")}

nm=length(modifiers)
if (nm>0){cat("modifiers=c(", file=fid, sep="")
for (k in 1:nm) cat(sprintf("\"%s\"%s",modifiers[k], ifelse(k==nm,"),\n",",")  ), file=fid, sep="")}

np=length(products)
if (np>0){cat("products=c(", file=fid, sep="")
for (k in 1:np) cat(sprintf("\"%s\"%s",products[k], ifelse(k==np,"),\n",",")  ), file=fid, sep="")}

npa=length(parameters)
if (npa>0){cat("parameters=c(", file=fid, sep="")
for (k in 1:npa) cat(sprintf("%s = %g%s",names(parameters)[k],as.numeric(parameters[k]), ifelse(k==npa,"),\n",",")  ), file=fid, sep="")}

# always have a law
cat(sprintf("strLaw=\"%s\"\n",strLaw), file=fid, sep="")


#if (funcLaw) {
#cat("law=", file=fid, sep=" ")
#tt=deparse(law,width=500)
#write(tt,fid)
#} 

cat(sprintf("%s",ifelse(i==nReactions,")) #END of the Reaction Section\n\n","),\n")  ), file=fid, sep="\n")
}  # end for loop on reactions
} # end if there are reactions

# State Units in Notes from here on out!
#
#if(nUnits>0){
#cat("units=c(", file=fid, sep="")
#for (i in 1:nUnits) 
#cat(sprintf("\"%s\"%s", units[i],ifelse(i==nUnits,")\n",",") ), file=fid, sep="")
#}
#  

cat(")  # END of the Model!", file=fid, sep="\n")

close(fid)
} # end read.SBML function definition




#if(nRules>0){
#cat("rules=list(", file=fid, sep="\n")
#for (i in 1:nRules) 
#{
#output=rules[[i]][["output"]][[1]]
#fcn=rules[[i]][["law"]][[2]]
#last=fcn[[length(fcn)]]
#recurs(last)
#}
#}
#
#for (i in 1:n )
#{
#math=rules[[i]][["math"]][[1]]
##print(math)
#leaves<-getRuleLeaves(math)
##print(leaves)
#inputs=setdiff(leaves,parameters)
##print(inputs)
#
##print(paste(getRuleLeaves(math),collapse=","))
#cat(sprintf("list(output=\"%s\", inputs=c(", xmlAttrs(rules[[i]])["variable"]),file=fid, sep="")
#ni=length(inputs)
#for (j in 1:ni) 
#  cat(sprintf("\"%s\"%s", inputs[j],ifelse(j==ni,"),\n",",")), file=fid, sep="")
## always have a law
#cat("law=function(r){", file=fid, sep="\n")
#if (ni>0) for (k in 1:ni) cat(sprintf("%s=r[\"%s\"];",inputs[k],inputs[k]), file=fid, sep="")
#cat(" ", file=fid, sep="\n")
#cat(paste(recurs(math),collapse=""), file=fid, sep="")
#cat(sprintf("}  \n %s",ifelse(i==n,")\n), \n","),")  ), file=fid, sep="\n")
#
#} # end for loop through the rules
#} # end if rule
#}
#