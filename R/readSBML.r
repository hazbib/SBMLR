"readSBML" <-
  function(filename)
  {  # takes SBML in filename.xml and maps it to a SBML class model 
    # using both Sax and DOM (for mathml) based parsing.
    
    sbmlHandler <- function ()   
    { # first block here sets up the parent environment used by all handler functions
      sbml<-"x"     # "x" is just a starting string value
      modelid<-"x"  				#storing model id
      lnotes<-NULL
      compartments <- list()
      reactLaws <- list()
      species <- list()
      rules<- list()
      reactions<- list()
      globalParameters=list()   
      reactants=NULL
      products=NULL
      modifiers=NULL
      currRxnID=NULL
      parameters=NULL   # local to rate law
      parameterIDs=NULL   # local to rate law
      globalParameterIDs=NULL   
      
      notes=FALSE; reactant=FALSE; product=FALSE
      law=FALSE; parameter=FALSE; math=FALSE
      
      # 		# VV's additions
      modelname<-"x"  				#storing the model name if separately given
      paramException <- 1
      globalParamException <- 1
      ParametersList<- list()  			#list of parameter objects	
      
      
      .startElement <- function(name, atts, ...) {
        #   cat("Start: Name =",name," ",paste(names(atts),atts,sep=" = "),"\n")
        if(name=="sbml")
        { 
          sbml<<-atts
        }
        
        #			if(name=="annotation")  print("skipping annotation") 
        
        # TR version
        #			if(name=="model")  {modelid<<-atts[["id"]]} # VV replaces this with ...
        
        # VV version
        # if(name=="model")  
        # {   
        #   numitems <- length(atts)
        #   if(numitems < 1)			#if model does not contain a name/id, we give it an arbitrary one.
        #   {  
        #     modelid[[1]]<<-"BioModel"	  
        #   }
        #   else if(numitems == 1)				# if only one attribute supplied
        #   {
        #     if(is.character(atts[[1]]))	#if the attribute is a string, it must be model name. Copy to also id.
        #     {
        #       modelname<<-atts[[1]]   # store model name
        #       modelid<<-atts[[1]]     # store as model id (copy Essentially)
        #     }
        #   }   
        #   else if(numitems > 1)			#both Id and name of model are supplied, read both
        #   {
        #     modelname<<-atts[["name"]]   # first element is model name
        #     modelid <<-atts[["id"]]    # second element has to be model id
        #   }
        # }
        
        
        # VP version        
        if(name == "model") {
          modelname <<- "no name"
          modelid <<- "no id"
          # the expected attributes are: name and id
          if("name" %in% names(atts))
            modelname <<- atts[["name"]]
          if("id" %in% names(atts))
            modelid <<- atts[["id"]]
        }

        
        
        #       if(name=="compartment")  {compartments[[atts[1]]]<<-atts}
        if(name=="compartment")  
        {
          values <<- names(atts)
          #			  
          if( "id" %in% values) compartments[[atts["id"]]]<<-atts
          
          #				if( "id" %in% names(atts)) compartments[[atts["id"]]]<<-atts
        }
        
        #       if(name=="species")  {species[[atts[1]]]<<-atts }
        if(name=="species") 
        {
          values <<- names(atts)
          if( "id" %in% values) species[[atts["id"]]]<<-atts 
          #				if( "id" %in% names(atts)) species[[atts["id"]]]<<-atts 
        }
        
        #       if(name=="assignmentRule")  rules[[atts[1]]]$idOutput<<-atts[[1]] 
        if(name=="assignmentRule")  
        {  
          rules[[atts["variable"]]]$idOutput<<-atts[["variable"]] 
        }
        #       if(name=="reaction")  {reactions[[atts[1]]]$id<<-atts[[1]]
        #         reactions[[atts[1]]]$reversible<<-as.logical(atts[[2]])
        #         currRxnID<<-atts[1]}
        if(name=="reaction")  
        {
          lstnames <- names(atts)
          numitems <- length(lstnames)
          nameslist <- list()
          id <- "x"
          reverse <- FALSE
          name <- "x"
          count <- 1
          while( count <= numitems )
          {
            switch(lstnames[[count]],
                   "id" = { id = atts[[count]]; nameslist[[length(nameslist)+1]] <- "id"}, 
                   "reversible" = { reverse = as.logical(atts[[count]]) ;nameslist[[length(nameslist)+1]] <- "reversible" },
                   "name" = { name = as.character(atts[[count]]); nameslist[[length(nameslist)+1]] <- "name"}
            )
            count <- count + 1
          } 
          
          reactions[[atts["id"]]]$id<<-id
          reactions[[atts["id"]]]$reversible<<-reverse  
          #           if(reverse) reactions[[atts["id"]]]$reversible<<-reverse #carry in R only  if true
          currRxnID<<-atts["id"]
          #reactions[[atts["id"]]]$id<<-atts[["id"]]
          #reactions[[atts[1]]]$reversible<<-as.logical(atts[[2]])
          #currRxnID<<-atts[1]
        }
        
        #		  print ("error2")
        if(name=="listOfReactants")  reactant<<-TRUE
        if(name=="listOfProducts")  product<<-TRUE
        if(name=="kineticLaw")  law<<-TRUE
        if(name=="math")  math<<-TRUE
        if((name=="speciesReference")&reactant)
          reactants<<-c(reactants,species=atts[["species"]])
        #         reactants<<-c(reactants,species=atts[[1]])
        if((name=="speciesReference")&product)
          products<<-c(products,species=atts[["species"]])
        #         products<<-c(products,species=atts[[1]])
        if(name=="modifierSpeciesReference")
          modifiers<<-c(modifiers,species=atts[["species"]])
        
        #			      if((name=="parameter")&law){
        #			        parameterIDs<<-c(parameterIDs,atts[["id"]])
        #			        parameters<<-c(parameters,atts[["value"]])}
        if((name=="parameter")&law)  		#parameter encountered within a kinetic law definition
        {
          values <- names(atts)
          if( "id" %in% values) parameterIDs<<-c(parameterIDs,atts[["id"]])
          else {
            cat('Parameter parsed without id. Setting default id', '\n')
            parameterIDs<<-c(parameterIDs, paste("default", paramException))
            paramException <- paramException + 1
          }
          if( "value" %in% values) parameters<<-c(parameters,atts[["value"]])
          else {
            cat('Warning..Parsing parameter without value. Setting it as 0.' ,'\n')
            parameters<<-c(parameters, as.numeric(0))
          }
        }
        
        #			      if((name=="parameter")&(!law)){
        #			        globalParameterIDs<<-c(globalParameterIDs,atts[["id"]])
        #			        globalParameters<<-c(globalParameters,as.numeric(atts[["value"]]))}
        
        if((name=="parameter")&!law)  		#parameter encountered outside a kinetic law definition - So in globalparamslist
        {
          #cat("within parameters:", atts[["id"]], atts[["value"]], "\n")
          values <- names(atts)
          if( "id" %in% values) {
            globalParameterIDs<<-c(globalParameterIDs,atts[["id"]])
            ParametersList[[atts["id"]]] <<- atts		#our new list of Parameter Objects
          }  
          else {
            cat('Global Parameter parsed without id. Setting default id', '\n')
            tempParamId <- paste("Globaldefault", globalParamException)
            globalParameterIDs<<-c(globalParameterIDs, tempParamId)
            ParametersList[[tempParamId]] <<- atts
            globalParamException <- globalParamException + 1
          }
          if( "value" %in% values) globalParameters<<-c(globalParameters,as.numeric(atts[["value"]]))
          else {
            cat('Warning..Parsing Global parameter without value. Setting it to 0.', '\n')
            globalParameters<<-c(globalParameters, as.numeric(0))
          }
          # 			} # end if param in law
        }
      } # end .startElement()  
      
      
      .endElement <- function(name) {
        if(name=="listOfReactants")  reactant<<-FALSE  
        if(name=="listOfProducts")  product<<-FALSE   
        if(name=="kineticLaw")  law<<-FALSE
        if(name=="math")  math<<-FALSE
        if((name=="listOfParameters")&(!law)) names(globalParameters)<<-globalParameterIDs  
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
      
      .text <- function(x, ...) {
        if (!math) lnotes<<-c(lnotes,x)
        #  cat("Txt:", x,"\n")
      }
      
      getModel<- function()       
      { 
        #       fixComps=function(x) {
        #         lst=list(x[[1]],as.numeric(x[[2]]) ); 
        #         names(lst)<-names(x); 
        #         lst 
        #       }
        #  VV replaces fixComps with the following:
        #			fixComps=function(x) 
        #			{
        #				lstnames <- names(x)
        #				count <- 1
        #				numit <- length(lstnames)
        #				id <- "x"
        #				size <- 0
        #				name <- "x"
        #				nameslist <- list()
        #				while( count <= numit )
        #				{
        #					switch(lstnames[[count]],
        #							"id" = { id = x[[count]]; nameslist[[length(nameslist)+1]] <- "id"}, 
        #							"size" = { size = as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "size" },
        #							"name" = { name = as.character(x[[count]]); nameslist[[length(nameslist)+1]] <- "name"}
        #					)
        #					count = count + 1
        #				}
        #				#cat("Compartment id:", id,"\n", "Compartment size:" , size, "\n","Compartment name:", name, "\n")
        #				if(numit == 2)						# only 2 attributes present. We need to find them.
        #				{ if(id == "x")		#id not set but name and size are.
        #						id <- "default"		
        #					else if(name == "x")    #name not set, we copy the id.
        #						name <- id
        #					else if(size== "0")	#size not set
        #						size <- 1 	#arbitrary setting as 1
        #					lst = list(id,size,name)	
        #					names(lst)<-c("id","size","name")
        #					lst
        #				} else if(numit == 3)					# 3 attributes/items present.
        #				{ lst = list(id,size,name)
        #					names(lst)<-c("id","size","name")
        #					lst 
        #				}
        #			}
        
        fixComps <- function(x){
          out <- as.list(x) 
          out$size <- as.numeric(out$size)  # convert size to numeric 
        }
        
#               fixSpecies=function(x) {
#                 lst=list(x[[1]],as.numeric(x[[2]]),x[[3]],as.logical(x[[4]])); 
#                 names(lst)<-c("id","ic","compartment","bc"); 
#                 lst 
#               }
        #  VV replaces fixSpecies with the following
#        fixSpecies=function(x) 
#        {
#          #cat (names(x), "\n")
#          #cat(toString(x) , "\n")
#          numitems <- length(x)
#          lstnames <- names(x)
#          count <-1
#          id <- "x"			#species Id
#          ic <- 0				#species initial concentration
#          compart <- "def"		#species compartment
#          bc <- FALSE			#species boundary condition
#          name <- "def"
#          nameslist <- list()
#          while( count <= numitems)
#          {
#            switch(lstnames[[count]],
#                   "id" = { id <- x[[count]]; nameslist[[length(nameslist)+1]] <- "id"},
#                   "name" = { name <- x[[count]]; nameslist[[length(nameslist)+1]] <- "name"},
#                   "initialConcentration" = { ic <- as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "ic" },
#                   "compartment" = { compart <- as.character(x[[count]]); nameslist[[length(nameslist)+1]] <- "compartment"},
#                   "boundaryCondition" = { bc <- as.logical(x[[count]]); nameslist[[length(nameslist)+1]] <- "bc"}
#            )
#            count = count + 1
#          }
#          #lst = list(id,ic,compart,bc, name)
#          lst = list(id,as.numeric(ic), compart, as.logical(bc))
#          names(lst) <- c("id","ic","compartment","bc")
#          #names(lst)<-c("id","ic","compartment","bc", "name"); 
#          lst 
#        }
        
        # VP version 
        
        fixSpecies = function(x){
           out <- as.list(x)
#          
          # sorting out amount and concentration
          if("initialAmount" %in% names(out))
            out$initialAmount <- as.numeric(out$initialAmount)
         if("initialConcentration" %in% names(out)){
            out$initialConcentration <- as.numeric(out$initialConcentration)      #!!!!! Error: Object 'initialconcentration' not found 
          }else{
          out$initialConcentration <- out$initialAmount/compartments[[out$compartment]]
            
          }
          out$ic <- out$initialConcentration
          
          # boundary condition
          # by default boundary condition should be false 
          if(!("boundaryCondition" %in% names(out)))
            out$boundaryCondition <- "false"
          
          out$bc <- c("true"=TRUE, "false"=FALSE)[out$boundaryCondition]
          
          return(out)
          
        }
        
        
        # and VV adds in fixParams
        fixParams=function(x) 
        {
          numitems <- length(x)
          lstnames <- names(x)
          count <-1
          id <- "x"			#Parameter Id
          value <- 0			#Parameter value
          name <- "def"
          constant <- FALSE
          nameslist <- list()
          while( count <= numitems)
          {
            switch(lstnames[[count]],
                   "id" = { id <- x[[count]]; nameslist[[length(nameslist)+1]] <- "id"},
                   "name" = { name <- x[[count]]; nameslist[[length(nameslist)+1]] <- "name"},
                   "value" = { value <- as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "value" },
                   "constant" = { constant <- as.logical(x[[count]]) ; nameslist[[length(nameslist)+1]] <- "constant"}
            )
            count = count + 1
          }
          
          lst = list(id,as.numeric(value))
          names(lst) <- c("id","value")
          lst 
        }		
        
#        compartments <- sapply(compartments,fixComps, simplify = FALSE)
        #species=t(sapply(species,fixSpecies, simplify = TRUE)[2:4,]) # this changes the species model structure for better looks in R dumps
#        species=sapply(species,fixSpecies, simplify = FALSE)     # this keeps the better looks in the SBMLR model definition file
        # 			ParametersList = sapply(ParametersList, fixParams, simplify = FALSE)  #VV: building params list 
        compartments <- lapply(compartments, fixComps)
        species <- lapply(species, fixSpecies)
        
        
        
        # 			list(sbml=sbml,id=modelid[[1]], notes=lnotes,compartments=compartments, # VV, not clear how ParametersList differs from globalParameters
        # 					species=species,globalParameters=globalParameters, ParametersList=ParametersList, rules=rules,reactions=reactions)
        list(sbml=sbml,id=modelid[[1]], notes=lnotes,compartments=compartments, # TR may revert to this??
             species=species,globalParameters=globalParameters, rules=rules,reactions=reactions) # returns values accrued in parent env
      } # end of getModel
      
      list(.startElement = .startElement, .endElement = .endElement, 
           .text = .text,   # , dom = function() {con}
           getModel = getModel     
      ) # function returns a list of functions, each with a common parent environment = stuff before function definitions
    }
    
    #  END handler definition
    
    # NOTE: though handlers are neat, one must question if the added baggage is worth it, i.e. compare to read.SBML of older versions 
    
    # *********************************************************************************
    # The next block of three functions converts mathML XMLnode objects into R expression objects
    # This approach is better than the old read.SBML approach in that the parenthesis overkill is avoided!
    mathml2R <-function(node)  {UseMethod("mathml2R", node)}
    
    mathml2R.XMLDocument <-function(doc) {return(mathml2R(doc$doc$children))}
    
    #     mathml2R.default<-function(children) 
    #     {  expr <- expression()  # this gets used when a "list" of children nodes are sent in
    #        for(i in children)  expr=c(expr, mathml2R(i)) 
    #        return(expr)
    #     }
    
    mathml2R.default<-function(children) 
    {  expr <- expression()  # this gets used when a "list" of children nodes are sent in
    n=length(children)
    #    cat("into default length n is ",n,"\n")
    #    for(i in children)  expr=c(expr, mathml2R(i)) 
    for(i in 1:n)  expr=c(expr, mathml2R(children[[i]])) 
    if (n>3) {#print("n>3")  # this fixes libsbml problem that times is not binary
      if (expr[[1]]=="*") expr[[1]]=as.name("prod") # in R, prod takes arb # of args
      if (expr[[1]]=="+") expr[[1]]=as.name("sum")  # similary for sum
    }
    #    print(children)
    #    print(expr)
    #    print("leaving default")
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
    
    
    # The next two functions are used by rules and were taken straight from read.SBML
    # The idea is that SBML doesn't provide a list of atoms/leaves with rules, so we have to create them
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
    model=edoc$getModel() # SAX approach using the handler. Output of getModel() in edoc list is what we want.
    doc <- xmlTreeParse(filename,ignoreBlanks = TRUE)  # use DOM just for rules and reactions
    model$htmlNotes=doc$doc$children$sbml[["model"]][["notes"]] 
    rules=doc$doc$children$sbml[["model"]][["listOfRules"]]
    reactions=doc$doc$children$sbml[["model"]][["listOfReactions"]]
    
    globalParameters=names(model$globalParameters)
    
    
    nRules=length(rules)
    #print(nRules)
    #nRules=0
    if (nRules>0){
      #    ruleIDs=NULL
      for (i in 1:nRules)    #  for( i in 1:(nRules-1))   # VV stops 1 shy of end????
      {  # assume they are assignment rules
        mathml<-rules[[i]][["math"]][[1]]
        model$rules[[i]]$mathmlLaw=mathml
        e<-mathml2R(mathml)
        model$rules[[i]]$exprLaw<-e[[1]]
        model$rules[[i]]$strLaw<-gsub(" ","",toString(e[1]))
        leaves<-getRuleLeaves(mathml)
        r<-model$rules[[i]]$inputs<-setdiff(leaves,globalParameters) # must deduce inputs by substracting global params
        #      r=model$rules[[i]]$inputs
        #model$rules[[i]]$idOutput=xmlAttrs(rules[[i]])[["variable"]][[1]] # moved to handler to preset the rules list length!
        model$rules[[i]]$law=makeLaw(r,NULL,model$rules[[i]]$exprLaw)
        #      ruleIDs[i]<-model$rules[[i]]$idOutput
      }
      #    names(model$rules)<-sapply(model$rules,function(x) x$idOutput)
    } 
    
    nReactions=length(reactions)
    if (nReactions>0){
      #    rIDs=NULL;  
      for (i in 1:nReactions)
      {
        model$reactions[[i]]$mathmlLaw=reactions[[i]][["kineticLaw"]][["math"]][[1]]
        e=mathml2R(reactions[[i]][["kineticLaw"]][["math"]][[1]])
        model$reactions[[i]]$exprLaw=e[[1]]
        #print(e[[1]])
        #print(toString(e[1]))
        model$reactions[[i]]$strLaw=gsub(" ","",toString(e[1]))
        #paste(as.character(model$reactions[[i]]$expr)[c(2,1,3)],collapse=""))
        # r=model$reactions[[i]]$reactants
        #  VV wants to add products in here, perhaps for reversible reactions
        # VP - yes, for reversible reaction we need values of products too
        r <- c(model$reactions[[i]]$reactants, model$reactions[[i]]$products)	#build using both reactants and products objects. TODO - add compartments
        p=names(model$reactions[[i]]$parameters)
        m=model$reactions[[i]]$modifiers
        e=model$reactions[[i]]$exprLaw
        model$reactions[[i]]$law=makeLaw(c(r,m),p,e, compartments=model$compartments)
#        model$reactions[[i]]$law = makeLaw(c(r, m), p, e)
        #      rIDs[i]<-model$reactions[[i]]$id
      }
      # This is for indexing by names/IDs of reactions
      #    names(model$reactions)<-rIDs
      #    names(model$reactions)<-sapply(model$reactions,function(x) x$id)
    }
    
    #---DEBUG CODE
    cat("Number of species: " , length(model$species), "\n")
    cat("Number of rules: ", nRules, "\n")
    cat("Number of Global Parameters: " , length(globalParameters), "\n")
    cat("Number of reactions: " , nReactions, "\n")
    cat("Parsing Successful !" , "\n")
    #-----------
    
    class(model)<-"SBMLR"
    model
  }

# the following is called by both readSBML and readSBMLR so it outside where both can reach it.
# Note that keeing it here instead of in a separate file => no need to document it

#"makeLaw"<-function(r,p,e, compartments = NULL){
#    attach(compartments)
#  # takes reactant list r, parameter list p and rate law R expression e 
#  # and makes a reaction rate law function out of them.
#  lawTempl=function(r,p=NULL){ }
#  i=2
#  for (j in seq(along=p)){
#    #		if(!is.null(p))
#    #		for (j in 1:length(p)){
#    body(lawTempl)[[i]]<-call("=",as.name(p[j]),call("[",as.name("p"),p[j]))
#    i=i+1}
#  #   for (j in 1:length(r)){ 
#  for (j in seq(along=r)){
#    body(lawTempl)[[i]]<-call("=",as.name(r[j]),call("[",as.name("r"),r[j]))
#    i=i+1}
#  body(lawTempl)[[i]]<-e
#  lawTempl
#}

function (r, p, e, compartments = NULL) 
{
  attach(compartments)
  .. <- c(r, p)
  lawTempl = function(r, p = NULL) {
  }
  i = 2
  body(lawTempl)[[i]] <- call("<-", as.name(".."), call("c", as.name("p"), as.name("r")))
  for (j in seq_along(..)) {
    body(lawTempl)[[i]] <- call("=", as.name(..[j]), call("[", as.name(".."), ..[j]))
    i <- i + 1
  }
  body(lawTempl)[[i]] <- e
  return(lawTempl)
}



#
#"makeLaw"<-function(inputs,pars,lawCall)
#	function(inputs,lawCall,pars=NULL)
#	with(as.list(c(inputs,pars)),lawCall)
#



