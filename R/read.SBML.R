"read.SBML" <-
function(filename)
{  # takes SBML in filename.xml and maps it to a model definition filenameX.r
if(!require(XML)) print("Error in Read.SBML(): First Install the XML package http://www.omegahat.org/RSXML/XML_0.95-6.zip")

ML2R<- function(type)   # map MathML operator symbols into R symbols
       switch(type,
             "times" = "*",
             "divide" = "/",
             "plus" = "+",
             "minus" = "-",
             "power" = "^",
             "exp" = "exp",
              "ln" = "log",
              "not found") # end definition of ML2R

recurs<-function(math) 
{ n=length(math)
 S=c("(")
 op=ML2R(xmlName(math[[1]]))
 for (j in 2:n )
 if (xmlName(math[[j]])=="ci") S=c(S,as.character(xmlValue(math[[j]])), ifelse(j==n,"",op)) else 
 S=c(S,Recall(math[[j]]), ifelse(j==n,"",op))
S=c(S,")")
S
} # end recursive function definition


getRuleLeaves<-function(math) 
{ n=length(math)
 S=c(NULL)
 op=ML2R(xmlName(math[[1]]))
 for (j in 2:n )
 if (xmlName(math[[j]])=="ci") S=c(S,as.character(xmlValue(math[[j]]))) else 
 S=c(S,Recall(math[[j]])  ) 
S
} # end recursive function definition


# start main part of read.SBML functino definition
fid <- file(paste(filename,"X.r",sep=""), "w")  # open the output file connection
cat("#",filename,"X.r\n", file=fid, sep="")  # X indicates R via XML
cat("model=list( ", file=fid, sep="\n")

# now read in a SBML file 
doc <- xmlTreeParse(paste(filename,".xml",sep=""),getDTD=FALSE)

cat("notes=c(", file=fid, sep="\n")
notes=doc$children$sbml[["model"]][["notes"]][["body"]]
n=length(xmlChildren(notes))
for (i in 1:n)
cat(sprintf("\"%s\"%s", xmlValue(xmlChildren(notes)[[i]]),ifelse(i==n,"\n),\n",",")), file=fid, sep="\n")

cat("comps=list(", file=fid, sep="\n")
comps=xmlChildren(doc$children$sbml[["model"]][["listOfCompartments"]])#;comps
n=length(comps)
for (i in 1:n )
cat(sprintf("list(id=\"%s\", vol = %g)%s", xmlAttrs(comps[[i]])["id"],as.numeric(xmlAttrs(comps[[i]])["size"]),ifelse(i==n,"\n),\n",",")), 
file=fid, sep="\n")

cat("species=list(", file=fid, sep="\n")
species=xmlChildren(doc$children$sbml[["model"]][["listOfSpecies"]])#;species
n=length(species)
for (i in 1:n)
cat(sprintf("%s =list( id=\"%s\", ic=%g,  comp=\"%s\", bc=%s)%s",xmlAttrs(species[[i]])[["id"]],xmlAttrs(species[[i]])[["id"]],
as.numeric(xmlAttrs(species[[i]])[["initialConcentration"]]),xmlAttrs(species[[i]])[["compartment"]],
ifelse(xmlAttrs(species[[i]])[["boundaryCondition"]]=="true","TRUE","FALSE"),ifelse(i==n,"\n),\n",",")),file=fid, sep="\n")

params=NULL
if(!is.null(doc$children$sbml[["model"]][["listOfParameters"]]))
{
parameters=xmlChildren(doc$children$sbml[["model"]][["listOfParameters"]])
n=length(parameters)
if (n>0){
cat("parameters=list(", file=fid, sep="")
for (i in 1:n )
{cat(sprintf("%s=%g%s", xmlAttrs(parameters[[i]])["id"],as.numeric(xmlAttrs(parameters[[i]])["value"]),ifelse(i==n,"),\n",",") ), file=fid, sep="\n")
params=c(params,xmlAttrs(parameters[[i]])["id"])}
}
}

if(!is.null(doc$children$sbml[["model"]][["listOfRules"]]))  # TR 5/27/04. catch missing rules before xmlChildren call
{
rules=xmlChildren(doc$children$sbml[["model"]][["listOfRules"]])
n=length(rules)
if (n>0){
cat("rules=list(", file=fid, sep="\n")
for (i in 1:n )
{
math=rules[[i]][["math"]][[1]]
#print(math)
leaves<-getRuleLeaves(math)
#print(leaves)
inputs=setdiff(leaves,params)
#print(inputs)

#print(paste(getRuleLeaves(math),collapse=","))
cat(sprintf("list(output=\"%s\", inputs=c(", xmlAttrs(rules[[i]])["variable"]),file=fid, sep="")
ni=length(inputs)
for (j in 1:ni) 
  cat(sprintf("\"%s\"%s", inputs[j],ifelse(j==ni,"),\n",",")), file=fid, sep="")
# always have a law
cat("law=function(r){", file=fid, sep="\n")
if (ni>0) for (k in 1:ni) cat(sprintf("%s=r[\"%s\"];",inputs[k],inputs[k]), file=fid, sep="")
cat(" ", file=fid, sep="\n")
cat(paste(recurs(math),collapse=""), file=fid, sep="")
cat(sprintf("}  \n %s",ifelse(i==n,")\n), \n","),")  ), file=fid, sep="\n")

} # end for loop through the rules
} # end if rule
}


cat("rxns=list(", file=fid, sep="\n")
rxns=xmlChildren(doc$children$sbml[["model"]][["listOfReactions"]])
n=length(rxns)
for (i in 1:n)
{
cat(sprintf("list( id=\"%s\", rever=%s,",xmlAttrs(rxns[[i]])[["id"]],ifelse(xmlAttrs(rxns[[i]])[["reversible"]]=="true","TRUE","FALSE")),file=fid, sep="\n")
law=rxns[[i]]["kineticLaw"][[1]]
params=law[["listOfParameters"]]
math=law[["math"]][[1]]
reacts=rxns[[i]]["listOfReactants"][[1]]
mods=rxns[[i]]["listOfModifiers"][[1]]
prods=rxns[[i]]["listOfProducts"][[1]]

nr=length(reacts)
if (nr>0){cat("reacts=c(", file=fid, sep="")
for (k in 1:nr) cat(sprintf("\"%s\"%s",xmlAttrs(reacts[[k]]), ifelse(k==nr,"),\n",",")  ), file=fid, sep="")}

nm=length(mods)
if (nm>0){cat("mods=c(", file=fid, sep="")
for (k in 1:nm) cat(sprintf("\"%s\"%s",xmlAttrs(mods[[k]]), ifelse(k==nm,"),\n",",")  ), file=fid, sep="")}

np=length(prods)
if (np>0){cat("prods=c(", file=fid, sep="")
for (k in 1:np) cat(sprintf("\"%s\"%s",xmlAttrs(prods[[k]]), ifelse(k==np,"),\n",",")  ), file=fid, sep="")}

npa=length(params)
if (npa>0){cat("params=c(", file=fid, sep="")
for (k in 1:npa) cat(sprintf("%s = %g%s",xmlAttrs(params[[k]])[["id"]],as.numeric(xmlAttrs(params[[k]])[["value"]]), ifelse(k==npa,"),\n",",")  ), file=fid, sep="")}
# always have a law
cat("law=function(r,p){", file=fid, sep="\n")
if (npa>0) for (k in 1:npa) cat(sprintf("%s = p[\"%s\"];",xmlAttrs(params[[k]])[["id"]],xmlAttrs(params[[k]])[["id"]]  ), file=fid, sep="")
cat(" ", file=fid, sep="\n")
if (nr>0) for (k in 1:nr) cat(sprintf("%s = r[\"%s\"];",xmlAttrs(reacts[[k]]),xmlAttrs(reacts[[k]])), file=fid, sep="")
if (nm>0) for (k in 1:nm) cat(sprintf("%s = r[\"%s\"];",xmlAttrs(mods[[k]]),xmlAttrs(mods[[k]])), file=fid, sep="")
#if (np>0) for (k in 1:np) cat(sprintf("%s = r[\"%s\"];",xmlAttrs(prods[[k]]),xmlAttrs(prods[[k]])), file=fid, sep="")
# the line above may be needed for systems with reversible reactions since product concentrations would then enter the rate laws.
cat(" ", file=fid, sep="\n")
cat(paste(recurs(math),collapse=""), file=fid, sep="")
cat(sprintf("}  \n %s",ifelse(i==n,"))) \n","),")  ), file=fid, sep="\n")
}  # end for loop on reactions
close(fid)
} # end read.SBML function definition

