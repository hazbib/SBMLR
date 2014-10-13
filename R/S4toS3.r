"S4toS3" <-
  function(dom)
  { 
	 # takes S4 SBML object produced by rsbml and converts to S3 SBMLR object
	require(rsbml)
	compartments=compartments(model(dom))
	for (i in seq(along=compartments))
		compartments[[i]]=list(
				id=id(compartments[[i]]), 
				size=size(compartments[[i]])
	) 
#dput(compartments$cell )

	species=species(model(dom))
	for (i in seq(along=species))
		species[[i]]=list(
				id=id(species[[i]]), 
				ic=initialConcentration(species[[i]]), 
				compartment=compartment(species[[i]]),
				bc=boundaryCondition(species[[i]])
		) 
  
	parameters=parameters(model(dom))
  globals=lapply(parameters,value)
  names(globals)<-sapply(parameters,id)

	rules=rules(model(dom))
	for (i in seq(along=rules)) {
	  e=math(rules[[i]])
	  lawS=as.character(e)
    idOutput=variable(rules[[i]])
#     mathml=SBMLR:::R2MathML(e)
	  mathml=SBMLR:::R2MathML(e[[1]])$value()[[1]]
	  leaves=SBMLR:::getRuleLeaves(mathml)
    inputs<-setdiff(leaves,names(globals))
	  rules[[i]]=list(
	    idOutput=idOutput, 
      inputs=inputs,
	    strLaw=lawS,
	    exprLaw=e[[1]],
	    law=SBMLR:::makeLaw(inputs,NULL,e[[1]]),
      mathmlLaw=mathml
	  )
	  kill=which(sapply(rules[[i]],length)==0)
	  if (length(kill)>0) rules[[i]]=rules[[i]][-kill]
	}
	
  
  
	reactions=reactions(model(dom))
	for (i in seq(along=reactions)) {
		pnms=sapply(parameters(kineticLaw(reactions[[i]])),id)
		pars=sapply(parameters(kineticLaw(reactions[[i]])),value)
		names(pars)=pnms
		reacts=names(reactants(reactions[[i]]))
		mods=names(modifiers(reactions[[i]]))
		prods=names(products(reactions[[i]]))
		e=math(kineticLaw(reactions[[i]]))
		reactions[[i]]=list(
				id=id(reactions[[i]]), 
				reversible=reversible(reactions[[i]]),
				reactants=reacts, 
				modifiers=mods,
				products=prods,
				parameters=pars,
				strLaw=as.character(e),
		    exprLaw=e[[1]],
		   law=SBMLR:::makeLaw(c(reacts,mods),names(pars),e[[1]]),
				mathmlLaw=SBMLR:::R2MathML(e[[1]])$value()[[1]]
		)
		kill=which(sapply(reactions[[i]],length)==0)
		if (length(kill)>0) reactions[[i]]=reactions[[i]][-kill]
	}
	
    mod=list(
			 id=id(model(dom)),
			 notes= strsplit(notes(model(dom)),split="\n")[[1]],
			 compartments=compartments,
			 species=species,
			 reactions=reactions,
			 globalParameters=globals,
			 rules=rules
	 )
	 
	 class(mod)<-"SBMLR"
    mod
  }

