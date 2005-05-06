"readSBML" <-
function(filename)
{  # takes SBML in filename.xml and maps it to a SBML class model 
# using both Sax and DOM (for mathml) based parsing.

sbmlHandler <- function () 
{
 sbml<-"x"
 modelid<-"x"
 lnotes<-NULL
 compartments <- list()
 reactLaws <- list()
 species <- list()
 rules<- list()
 reactions<- list()
 reactants=NULL
 products=NULL
 modifiers=NULL
 currRxnID=NULL
 parameters=NULL   # local to rate law
 parameterIDs=NULL   # local to rate law
 globalParameters=list()   
 globalParameterIDs=NULL   

 notes=FALSE; reactant=FALSE; product=FALSE
 law=FALSE; parameter=FALSE; math=FALSE

 startElement <- function(name, atts, ...) {
#   cat("Start: Name =",name," ",paste(names(atts),atts,sep=" = "),"\n")
   if(name=="sbml")  {sbml<<-atts}
   if(name=="model")  {modelid<<-atts[[1]]}
   if(name=="compartment")  {compartments[[atts[1]]]<<-atts}
   if(name=="species")  {species[[atts[1]]]<<-atts }
   if(name=="assignmentRule")  {rules[[atts[1]]]$idOutput<<-atts[[1]] }
   if(name=="reaction")  {reactions[[atts[1]]]$id<<-atts[[1]]
                  reactions[[atts[1]]]$reversible<<-as.logical(atts[[2]])
                  currRxnID<<-atts[1]}
   if(name=="listOfReactants")  {reactant<<-TRUE}
   if(name=="listOfProducts")  {product<<-TRUE}
   if(name=="kineticLaw")  {law<<-TRUE}
   if(name=="math")  {math<<-TRUE}
   if((name=="speciesReference")&reactant){
          reactants<<-c(reactants,species=atts[[1]])}
   if((name=="speciesReference")&product){
          products<<-c(products,species=atts[[1]])}
   if(name=="modifierSpeciesReference"){
          modifiers<<-c(modifiers,species=atts[[1]])}
   if((name=="parameter")&law){
          parameterIDs<<-c(parameterIDs,atts[[1]])
          parameters<<-c(parameters,atts[[2]])}
   if((name=="parameter")&(!law)){
          globalParameterIDs<<-c(globalParameterIDs,atts[[1]])
          globalParameters<<-c(globalParameters,as.numeric(atts[[2]]))}
  }  
 
endElement <- function(name) {
   if(name=="listOfReactants")  {reactant<<-FALSE  }
   if(name=="listOfProducts")  {product<<-FALSE   }
   if(name=="kineticLaw")  {law<<-FALSE}
   if(name=="math")  {math<<-FALSE}
   if((name=="listOfParameters")&(!law)) {names(globalParameters)<<-globalParameterIDs } 
   if(name=="reaction")  {
    names(reactants)<<-NULL
    names(modifiers)<<-NULL
    names(products)<<-NULL
    reactions[[currRxnID]]$reactants<<-reactants
    reactions[[currRxnID]]$modifiers<<-modifiers
    reactions[[currRxnID]]$products<<-products
    parameters<<-as.numeric(parameters)
    names(parameters)<<-parameterIDs
    reactions[[currRxnID]]$parameters<<-parameters
    reactants<<-NULL;products<<-NULL
    modifiers<<-NULL;parameters<<-NULL;parameterIDs<<-NULL      }
 }
 
 text <- function(x, ...) {
  if (!math) lnotes<<-c(lnotes,x)
#  cat("Txt:", x,"\n")
  }

getModel<- function() 
{ 
fixComps=function(x) {lst=list(x[[1]],as.numeric(x[[2]]) );  names(lst)<-names(x); lst }
fixSpecies=function(x) {lst=list(x[[1]],as.numeric(x[[2]]),x[[3]],as.logical(x[[4]])); names(lst)<-c("id","ic","comp","bc"); lst }
#
compartments=sapply(compartments,fixComps, simplify = FALSE)
#species=t(sapply(species,fixSpecies, simplify = TRUE)[2:4,]) # this changes the species model structure for better looks in R dumps
species=sapply(species,fixSpecies, simplify = FALSE)     # this keeps the better looks in the SBMLR model definition file

list(sbml=sbml,id=modelid[[1]], notes=lnotes,compartments=compartments,species=species,globalParameters=globalParameters, rules=rules,reactions=reactions)}

    list(startElement = startElement, endElement = endElement, 
        text = text,   # , dom = function() {con}
         getModel = getModel     
        )
}

#  END handler definition

# NOTE: though handlers are neat, one must question if the added baggage is worth it, i.e. compare to read.SBML of older versions 

# *********************************************************************************
# The next block of three functions converts mathML XMLnode objects into R expression objects
# This approach is better than the old read.SBML approach in the parenthesis overkill is avoided!
mathml2R <-function(node)  {UseMethod("mathml2R", node)}

mathml2R.XMLDocument <-function(doc) {return(mathml2R(doc$doc$children))}

mathml2R.default<-function(children) 
{  expr <- expression()  # this gets used when a "list" of children nodes are sent in
  for(i in children) {    expr <- c(expr, mathml2R(i))  }
 return(expr)
}

mathml2R.XMLNode <-function(node){
 nm <- xmlName(node) 
# cat("XMLNode: node name is ",nm," and the node class is",class(node),"\n")
 if(nm=="power"||nm == "divide"||nm =="times"||nm=="plus"||nm=="minus") {
    op <- switch(nm, power="^", divide="/",times="*",plus="+",minus="-")
     val <- as.name(op)
     } else if((nm == "ci")|(nm == "cn")) {
       if(nm == "ci") val <- as.name(node$children[[1]]$value)
      if(nm == "cn") val <- as.numeric(node$children[[1]]$value)
       }  else if(nm == "apply") {
       val <- mathml2R(node$children)
       mode(val) <- "call"
      } else  {cat("error: nm =",nm," not in set!\n")}
 return(as.expression(val))
}
# ********** END the mathML2R block of method based on node type codes  *************************




#  These next to are used by rules and were taken straight from read.SBML
# The idea is that SBML doesn't provide a list of atoms/leavs with rules, so we have to create them
# to place them in their model slots, and to use them to create the R function definition for the rule
# using makeLaw with a null for parameters, since they are passed global for rules.
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


getRuleLeaves<-function(math) 
{ n=length(math)
 S=c(NULL)
 op=ML2R(xmlName(math[[1]]))
 for (j in 2:n )
 if ((xmlName(math[[j]])=="ci")|(xmlName(math[[j]])=="cn"))  S=c(S,as.character(xmlValue(math[[j]]))) else 
 S=c(S,Recall(math[[j]])  ) 
S
} 


if(!require(XML)) print("Error in Read.SBML(): First Install the XML package http://www.omegahat.org/RSXML")

edoc <- xmlEventParse(filename,handlers=sbmlHandler(),ignoreBlanks = TRUE)
model=edoc$getModel()
#print(model);next
doc <- xmlTreeParse(filename,ignoreBlanks = TRUE)
model$htmlNotes=doc$doc$children$sbml[["model"]][["notes"]]
rules=doc$doc$children$sbml[["model"]][["listOfRules"]]
reactions=doc$doc$children$sbml[["model"]][["listOfReactions"]]

#print(model$rules)
globalParameters=names(model$globalParameters)
#print(globalParameters)
#print("before rules")
#print(rules[[1]][[1]][[1]])

nSpecies=length(model$species)
if (nSpecies>0){
sIDs=NULL;
for (i in 1:nSpecies) sIDs[i]<-model$species[[i]]$id
# This is for indexing by names/IDs of reactions, species and rules
names(model$species)<-sIDs
}


nRules=length(rules)
#print(nRules)
#nRules=0
if (nRules>0){
ruleIDs=NULL
for (i in 1:nRules)
{  # assume they are assignment rules
mathml<-rules[[i]][["math"]][[1]]
model$rules[[i]]$mathmlLaw=mathml
e<-mathml2R(mathml)
model$rules[[i]]$exprLaw<-e[[1]]
model$rules[[i]]$strLaw<-gsub(" ","",toString(e[1]))
leaves<-getRuleLeaves(mathml)
model$rules[[i]]$inputs=setdiff(leaves,globalParameters)
r=model$rules[[i]]$inputs
#model$rules[[i]]$idOutput=xmlAttrs(rules[[i]])[["variable"]][[1]] # moved to handler to preset the rules list length!
model$rules[[i]]$law=makeLaw(r,NULL,model$rules[[i]]$exprLaw)
ruleIDs[i]<-model$rules[[i]]$idOutput
}
# This is for indexing by names/IDs of rules
names(model$rules)<-ruleIDs
} 

#next
nReactions=length(reactions)
if (nReactions>0){
rIDs=NULL;  
for (i in 1:nReactions)
{
model$reactions[[i]]$mathmlLaw=reactions[[i]][["kineticLaw"]][["math"]][[1]]
e=mathml2R(reactions[[i]][["kineticLaw"]][["math"]][[1]])
model$reactions[[i]]$exprLaw=e[[1]]
#print(e[[1]])
#print(toString(e[1]))
model$reactions[[i]]$strLaw=gsub(" ","",toString(e[1]))
#paste(as.character(model$reactions[[i]]$expr)[c(2,1,3)],collapse=""))
r=model$reactions[[i]]$reactants
p=names(model$reactions[[i]]$parameters)
m=model$reactions[[i]]$modifiers
r=c(r,m)
e=model$reactions[[i]]$exprLaw
model$reactions[[i]]$law=makeLaw(r,p,e)
rIDs[i]<-model$reactions[[i]]$id
}
# This is for indexing by names/IDs of reactions
names(model$reactions)<-rIDs
}





class(model)<-"SBML"
model}

