"write.SBML" <-
function(filename)
{  # takes R model in filename.r and maps it to an SBML model in filename.xml 

r2ml<- function(type)   # maps R operator symbols into MathML 
       switch(type,
             "*" = "<times/>",
             "/" = "<divide/>",
             "+" = "<plus/>",
             "-" = "<minus/>",
             "^" = "<power/>",
             "exp" = "<exp/>",
              "log" = "<ln/>",
              "not found")    # end of r2ml sub-function definition. 

fid <- file(paste(filename,".xml",sep=""), "w")  # open the output file connection
           
recurs<-function(last)  
{ if(last[[1]]=="(") {last=last[[2]]} # remove parentheses  
  cat("    <apply>", file=fid, sep="\n")
  cat(sprintf("      %s",r2ml(as.character(last[[1]])) ), file=fid, sep="\n")
  for (j in 2:length(last))
   if ((class(last[[j]])=="name")|(class(last[[j]])=="numeric"))  
    cat(sprintf("     <ci>%s</ci>",as.character(last[[j]]) ) , file=fid, sep="\n")  else Recall(last[[j]]) 
  cat("    </apply>", file=fid, sep="\n")
} # end of recurs sub-function definition
    
    
# now start the main part of write.SBML 
# read in the R file 
e=parse(file=paste(filename,".r",sep=""),n=-1)
notes=e[[1]][[3]][["notes"]]
comps=e[[1]][[3]][["comps"]]
species=e[[1]][[3]][["species"]]
parameters=e[[1]][[3]][["parameters"]]
rules=e[[1]][[3]][["rules"]]

rxns=e[[1]][[3]][["rxns"]]
units=e[[1]][[3]][["units"]]
nnotes=length(notes)
ncomps=length(comps)
nrxns=length(rxns)
nspcs=length(species)
nparameters=length(parameters)
nrules=length(rules)

cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file=fid, sep="\n")
cat("<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">", file=fid, sep="\n")
#print(sprintf("<model id=\"%s\">",filename))
cat(sprintf("<model id=\"%s\">",filename), file=fid, sep="\n")
if(nnotes>1){
cat("<notes>", file=fid, sep="\n")
cat(" <body xmlns=\"http://www.w3.org/1999/xhtml\">", file=fid, sep="\n")
for (i in 2:nnotes) cat(sprintf("   <p> %s  </p>",notes[[i]]), file=fid, sep="\n")
cat(" </body>", file=fid, sep="\n")
cat("</notes>", file=fid, sep="\n")
}
if(ncomps>1){
cat("<listOfCompartments>", file=fid, sep="\n")
for (i in 2:ncomps) cat(sprintf("   <compartment id=\"%s\"  size=\"%g\" />",comps[[i]][[2]],comps[[i]][[3]]), file=fid, sep="\n")
cat("</listOfCompartments>", file=fid, sep="\n")
}
if(nspcs>1){
cat("<listOfSpecies>", file=fid, sep="\n")
for (i in 2:nspcs) cat(sprintf("   <species id=\"%s\"  initialConcentration=\"%g\"  compartment=\"%s\" boundaryCondition=\"%s\"/>",
     species[[i]][[2]],species[[i]][[3]],species[[i]][[4]],ifelse(species[[i]][[5]],"true","false")), file=fid, sep="\n")
cat("</listOfSpecies>", file=fid, sep="\n")
}

if(nparameters>1){
cat("<listOfParameters>", file=fid, sep="\n")
for (i in 2:nparameters) cat(sprintf("   <parameter id=\"%s\"  value=\"%g\" />",names(parameters[i]),parameters[[i]]), file=fid, sep="\n")
cat("</listOfParameters>", file=fid, sep="\n")
}
if(nrules>1){
cat("<listOfRules>", file=fid, sep="\n")
for (i in 2:nrules) 
{
output=rules[[i]][["output"]][[2]]
cat(sprintf("  <assignmentRule variable=\"%s\">",output), file=fid, sep="\n")
#inputs=rules[[i]][["inputs"]]
#ni=length(inputs)
#cat("     <listOfInputs>", file=fid, sep="\n")
#     for (j in 2:ni) 
#	cat(sprintf("      <Input species = \"%s\" />", inputs[[j]]), file=fid, sep="\n")
#cat("     </listOfInputs>", file=fid, sep="\n")
cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
fcn=rules[[i]][["law"]][[3]]
last=fcn[[length(fcn)]]
recurs(last)
cat("    </math>", file=fid, sep="\n")     
cat("  </assignmentRule>", file=fid, sep="\n")     
}
cat("</listOfRules>", file=fid, sep="\n")
}

cat("<listOfReactions>", file=fid, sep="\n")
for (i in 2:nrxns) 
{cat(sprintf("  <reaction id=\"%s\"  reversible=\"%s\">",
rxns[[i]][["id"]],ifelse(rxns[[i]][["rever"]],"true","false")), file=fid, sep="\n")
reacts=rxns[[i]][["reacts"]]
if (!is.null(reacts[[2]])) {
cat("    <listOfReactants>", file=fid, sep="\n")
     for (j in 2:length(reacts)) 
	cat(sprintf("      <speciesReference species = \"%s\" />", reacts[[j]]), file=fid, sep="\n")
cat("    </listOfReactants>", file=fid, sep="\n")}

# Just switched the order of these two blocks to fix the errors
prods=rxns[[i]][["prods"]]
if (!is.null(prods[[2]])) {
cat("    <listOfProducts>", file=fid, sep="\n")
     for (j in 2:length(prods)) 
	cat(sprintf("      <speciesReference species = \"%s\" />", prods[[j]]), file=fid, sep="\n")
cat("    </listOfProducts>", file=fid, sep="\n")     }

mods=rxns[[i]][["mods"]]
if (!is.null(mods[[2]])) {
cat("    <listOfModifiers>", file=fid, sep="\n")
     for (j in 2:length(mods)) 
	cat(sprintf("      <modifierSpeciesReference species = \"%s\" />", mods[[j]]), file=fid, sep="\n")
cat("    </listOfModifiers>", file=fid, sep="\n")}


cat("  <kineticLaw>", file=fid, sep="\n")
fcn=rxns[[i]][["law"]][[3]]
last=fcn[[length(fcn)]]
cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
recurs(last)
cat("    </math>", file=fid, sep="\n")     
cat("    <listOfParameters>", file=fid, sep="\n")
params=rxns[[i]][["params"]]
for (j in 2:length(params)) 
if (length(params[[j]])==1) 
cat(sprintf("      <parameter id = \"%s\" value = \"%g\"/>", names(params)[j],params[[j]]), file=fid, sep="\n") else
cat(sprintf("      <parameter id = \"%s\" value = \"%g\"/>", names(params)[j],-params[[j]][[2]]), file=fid, sep="\n")
cat("    </listOfParameters>", file=fid, sep="\n")
cat("    </kineticLaw>", file=fid, sep="\n")
cat("  </reaction>", file=fid, sep="\n")
}
cat("</listOfReactions>", file=fid, sep="\n")
cat("</model>", file=fid, sep="\n")
cat("</sbml>", file=fid, sep="\n")
close(fid)
}  # end write.sbml function defition

