#sod.R
model=list( 
  
  id="sod",
  
  notes=c(
    "This is the SOD model of Kowald et al, JTB 238 (2006) 828-840.", 
    "Below, L is lipid radical, SO is superoxide, OH=hydroxyl radical,",
    "LOO is a radical, and LOOH is stable oxidized lipid.",
    "SOD is the total SOD, SODII=metal oxidized, SODII=reduced.",
    "HSO is protonated SO, a fast reaction; HSO neutral so it crosses membranes.",
    "Catalase. T and SODT are externally manipulated boundary conditions (bc).",
    "HSO and SODI are rule/algebraically determined. They are cloaked here as bc.",
    "Units are all in molar (M) and in seconds."
  ),
  
  
  compartments=list(
    list(id="cell", size = 1)
  ),
  
  
  species=list(
    list(id="SO"		,ic=       0,    compartment="cell",bc=  FALSE),
    list(id="SODII"	,ic=   5e-06,    compartment="cell",bc=  FALSE),
    list(id="H2O2"	,ic=       0,    compartment="cell",bc=  FALSE),
    list(id="OH"		,ic=       0,    compartment="cell",bc=  FALSE),
    list(id="L"		  ,ic=       0,    compartment="cell",bc=  FALSE),
    list(id="LOO"   ,ic=       0,    compartment="cell",bc=  FALSE),
    list(id="LOOH"	,ic=       0,    compartment="cell",bc=  FALSE),
#     list(id="HSO"		,ic=  0.0024,	   compartment="cell",bc=	 TRUE), # from BioModel
#     list(id="SODI"	,ic=    0.64,	   compartment="cell",bc=	 TRUE), #
    list(id="HSO"  	,ic=       0,	   compartment="cell",bc=	 TRUE), # from paper
    list(id="SODI"	,ic=    5e-6,	   compartment="cell",bc=	 TRUE), #
    list(id="SOD" 	,ic=   1e-05,    compartment="cell",bc=  TRUE),
    list(id="catalase"		,ic=   1e-05,    compartment="cell",bc=  TRUE)
  ),
  
   globalParameters=list(Keq=100), 
  
  
  #BEGIN the Rules Section
  rules=list(
    list(idOutput="HSO",
         inputs=c("SO"),
          strLaw = "SO/Keq"
#          strLaw = "SO/100"
    ),
    
    list(idOutput="SODI",
         inputs=c("SOD","SODII"),
         strLaw = "SOD-SODII"
    )),
  # END of the Rules Section
  
  
  #BEGIN the Reaction Section 
  reactions=list(
    list( id="etcV1",  
          parameters=c(k1=6.6e-7),
          products=c("SO"),
          strLaw="k1*1" # avoid a bug in my code by giving this an operator
    ),
    
    list( id="sooV2",  
          reactants=c("SO","SODII"),
          parameters=c(k2=1.6e9),
          strLaw="k2*SO*SODII"
    ),
    
    list( id="sorV3",  
          reactants=c("SO"),
          modifiers=c("SODI"),
          products=c("H2O2","SODII"),
          parameters=c(k3=1.6e9),
          strLaw="k3*SO*SODI"
    ),
    
    list( id="loorV4",  
          reactants=c("SO","LOO"),
          products=c("LOOH"),
          parameters=c(k4=1e5),
          strLaw="k4*SO*LOO"
    ),
    
    list( id="habWeiV5",  
          reactants=c("SO","H2O2"),
          products=c("OH"), # SO becomes O2, makes -OH anion also,  
          parameters=c(k5=2e4),   # neither of which are book keeped
          strLaw="k5*SO*H2O2"
    ),
    
    list( id="sodFenV6",  
          reactants=c("H2O2"),
          modifiers=c("SODII"),
          products=c("OH","OH"),
          parameters=c(k6=1),
          strLaw="k6*H2O2*SODII"
    ),
    
    list( id="catGpxV7",  
          reactants=c("H2O2"),
          modifiers=c("catalase"),
          parameters=c(k7=3.4e7),
          strLaw="k7*H2O2*catalase"
    ),
    
    list( id="ohSinkV9",  
          reactants=c("OH"),
          parameters=c(k9=1e6),
          strLaw="k9*OH"
    ),
    
    list( id="hsoLipPeroxV10",  
          products=c("L","H2O2"),
          modifiers=c("HSO"),
          parameters=c(k10=1000),
          strLaw="k10*HSO"
    ),

    
    list( id="hoLipPeroxV11",  
          reactants=c("OH"),
          products=c("L"),
          parameters=c(k11=2.5e+08),
          strLaw="k11*OH"
    ),
    
    list( id="loohDeOxV12",  
          reactants=c("LOOH"),
          parameters=c(k12=0.38),
          strLaw="k12*LOOH"
    ),
    
    list( id="sodoV13a",  
          products=c("SODII"),
          modifiers=c("SODI"),
          parameters=c(k13a=0.0087),
          strLaw="k13a*SODI"
    ),
    
    list( id="sodrV13b",  
          reactants=c("SODII"),
          parameters=c(k13b=0.0087),
          strLaw="k13b*SODII"
    ),

#     list( id="hsoSynFastV16",  # v16 also not shown in Fig. 1 
#           reactants=c("SO"),
#           parameters=c(k10=1000),
#           strLaw="k10*SO"
#     ),
    
  
    
    list( id="lipRadOxV17",  
          reactants=c("L"),
          products=c("LOO"),
          parameters=c(k17=30000),
          strLaw="k17*L"
    ),
    
    list( id="chainRxV18",  
          reactants=c("LOO"),
          products=c("L","LOOH"),
          parameters=c(k18=7),
          strLaw="k18*LOO"
    ),
    
  
    
    list( id="looSelfV19",  
          reactants=c("LOO","LOO"),
          parameters=c(k19=88000),
          strLaw="k19*LOO^2"
    )
    
    
  ) #END of the Reaction Section
  
  
)  # END of the Model!
