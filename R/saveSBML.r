"saveSBML" <-
  function(model,filename)
  {  # takes model object of class SBML ,and maps it to an SBML model in filename.xml 
    fid <- file(filename, "w")  # open the output file connection
    sbml=model[["sbml"]]
    id=model[["id"]]
    notes=model[["notes"]]
    htmlNotes=model[["htmlNotes"]]
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
    
    cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file=fid, sep="\n")
    cat("<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">", file=fid, sep="\n")
    #print(sprintf("<model id=\"%s\">",filename))
    cat(sprintf("<model id=\"%s\">",id), file=fid, sep="\n")
    if(nNotes>0){
      cat("<notes>", file=fid, sep="\n")
      cat(" <body xmlns=\"http://www.w3.org/1999/xhtml\">", file=fid, sep="\n")
      for (i in 1:nNotes) cat(sprintf("   <p> %s  </p>",notes[[i]]), file=fid, sep="\n")
      cat(" </body>", file=fid, sep="\n")
      cat("</notes>", file=fid, sep="\n")
    }
    if(nCompartments>0){
      cat("<listOfCompartments>", file=fid, sep="\n")
      for (i in 1:nCompartments) 
        cat(sprintf("   <compartment id=\"%s\"  size=\"%g\" name=\"%s\"/>",  # VV 
                    compartments[[i]][["id"]],compartments[[i]][["size"]],compartments[[i]][["name"]]), file=fid, sep="\n")
      #       cat(sprintf("   <compartment id=\"%s\"  size=\"%g\" />",compartments[[i]][[1]],compartments[[i]][[2]]), file=fid, sep="\n")
      cat("</listOfCompartments>", file=fid, sep="\n")
    }
    if(nSpecies>0){
      cat("<listOfSpecies>", file=fid, sep="\n")
      for (i in 1:nSpecies) 
        cat(sprintf("   <species id=\"%s\"  initialConcentration=\"%g\"  compartment=\"%s\" boundaryCondition=\"%s\"/>",
                    species[[i]][["id"]],species[[i]][["ic"]],species[[i]][["compartment"]],
                    ifelse(species[[i]][["bc"]],"true","false")), file=fid, sep="\n") #VV
      #               species[[i]][[1]],species[[i]][[2]],species[[i]][[3]],ifelse(species[[i]][[4]],"true","false")), file=fid, sep="\n")
      cat("</listOfSpecies>", file=fid, sep="\n")
    }
    if(nGlobalParameters>0){
      cat("<listOfParameters>", file=fid, sep="\n")
      for (i in 1:nGlobalParameters) 
        cat(sprintf("   <parameter id=\"%s\"  value=\"%g\" />",names(globalParameters[i]),globalParameters[[i]]), file=fid, sep="\n")
      cat("</listOfParameters>", file=fid, sep="\n")
    } 
    if(nRules>0){
      cat("<listOfRules>", file=fid, sep="\n")
      for (i in 1:nRules) 
      {
        idOutput=rules[[i]][["idOutput"]][[1]]
        cat(sprintf("  <assignmentRule variable=\"%s\">",idOutput), file=fid, sep="\n")
        cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
        ml=saveXML(rules[[i]]$mathmlLaw,prefix=NULL,file=fid)  
        cat("    </math>", file=fid, sep="\n")           
        cat("  </assignmentRule>", file=fid, sep="\n")   
# VV replaces the last 3 lines with this, but I don't get where newtry comes from????  
# In case of dirty xml file with messed up rules, we shouldn't write
#         if(!(is.null(newtry$rules[[i]]$mathmlLaw)))
#         {							# the messed up parts.
#           ml=saveXML(rules[[i]]$mathmlLaw,prefix=NULL,file=fid)  
#           cat("    </math>", file=fid, sep="\n")     
#           cat("  </assignmentRule>", file=fid, sep="\n")     
#         }
#         else
#         {
#           cat("    </math>", file=fid, sep="\n")     
#           cat("  </assignmentRule>", file=fid, sep="\n")     
#           next;						#move onto next rule
#         }
        
      }
      cat("</listOfRules>", file=fid, sep="\n")
    }
    
    cat("<listOfReactions>", file=fid, sep="\n")
    for (i in 1:nReactions) 
    {cat(sprintf("  <reaction id=\"%s\"  reversible=\"%s\">",
                 reactions[[i]][["id"]],ifelse(reactions[[i]][["reversible"]],"true","false")), file=fid, sep="\n")
     reactants=reactions[[i]][["reactants"]]
     if (!is.null(reactants[[1]])) {
       cat("    <listOfReactants>", file=fid, sep="\n")
       for (j in 1:length(reactants)) 
         cat(sprintf("      <speciesReference species = \"%s\" />", reactants[[j]]), file=fid, sep="\n")
       cat("    </listOfReactants>", file=fid, sep="\n")}
     
     # Just switched the order of these two blocks to fix the errors
     products=reactions[[i]][["products"]]
     if (!is.null(products[[1]])) {
       cat("    <listOfProducts>", file=fid, sep="\n")
       for (j in 1:length(products)) 
         cat(sprintf("      <speciesReference species = \"%s\" />", products[[j]]), file=fid, sep="\n")
       cat("    </listOfProducts>", file=fid, sep="\n")     }
     
     modifiers=reactions[[i]][["modifiers"]]
     if (!is.null(modifiers[[1]])) {
       cat("    <listOfModifiers>", file=fid, sep="\n")
       for (j in 1:length(modifiers)) 
         cat(sprintf("      <modifierSpeciesReference species = \"%s\" />", modifiers[[j]]), file=fid, sep="\n")
       cat("    </listOfModifiers>", file=fid, sep="\n")}
     
     cat("  <kineticLaw>", file=fid, sep="\n")
     cat("    <math xmlns=\"http://www.w3.org/1998/Math/MathML\">", file=fid, sep="\n")
     ml=saveXML(reactions[[i]]$mathmlLaw,prefix=NULL,file=fid) # annoying warnings were coming from here without file=fid
     cat("    </math>", file=fid, sep="\n")     
     #     cat("    <listOfParameters>", file=fid, sep="\n")
     #     parameters=reactions[[i]][["parameters"]]
     #     for (j in 1:length(parameters)) 
     #       if (length(parameters[[j]])==1) 
     #         cat(sprintf("      <parameter id = \"%s\" value = \"%g\"/>", names(parameters)[j],parameters[[j]]), file=fid, sep="\n") else
     #         cat(sprintf("      <parameter id = \"%s\" value = \"%g\"/>", names(parameters)[j],-parameters[[j]][[1]]), file=fid, sep="\n")
     #     cat("    </listOfParameters>", file=fid, sep="\n") # I forgot what the else is for, i.e. when is it not 1? is it >1 or 0???
     # VV replaces block above with this block. 
     parameters=reactions[[i]][["parameters"]]
     nlocalParameters = length(parameters)
     if(nlocalParameters > 0)			#if local parameters exist, we write it. Else we avoid it.
     { cat("    <listOfParameters>", file=fid, sep="\n")
       for (j in 1:nlocalParameters)		# Write each and every local parameter.
         cat(sprintf("      <parameter id = \"%s\" value = \"%g\"/>", names(parameters)[j],parameters[[j]]), file=fid, sep="\n")
       cat("    </listOfParameters>", file=fid, sep="\n")   }
     
     cat("    </kineticLaw>", file=fid, sep="\n")
     cat("  </reaction>", file=fid, sep="\n")
    }
    cat("</listOfReactions>", file=fid, sep="\n")
    cat("</model>", file=fid, sep="\n")
    cat("</sbml>", file=fid, sep="\n")
    close(fid)
  }  # end write.sbml function defition

