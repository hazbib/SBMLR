# The Morrison and Allegra Folate model (MorrisonAllegra.r)
model=list(
notes=c("This is a folate model that includes folate polyglutamation.", "Morrison and Allegra, JBC:264,10552-10566 (1989)",
        "Folate cycle kinetics in breast cancer cells",
        "Note: two flow BCs were converted into two downstream concentration BCs, thus removing the GAR and dUMP state variables.",
        "This dropped the number of ODEs from 21 to 19."),

comps=list(list(id= "cell",vol=1),list(id= "ext",vol=1)),

species=list( 
FH2f= 		list(id="FH2f", 	ic=0.0012,	comp="cell",bc=FALSE),
FH2b= 		list(id="FH2b", 	ic=0.0024,	comp="cell",bc=TRUE), # FH2 free and bound to DHFR are in rapid equilibrium
DHFRf= 		list(id="DHFRf", 	ic=0.64,	comp="cell",bc=FALSE),
DHFRtot= 	list(id="DHFRtot", 	ic=0.64,	comp="cell",bc=TRUE),
FH4= 		list(id="FH4", 		ic=0.46,	comp="cell",bc=FALSE),
CH2FH4=		list(id="CH2FH4", 	ic=0.26,	comp="cell",bc=FALSE),
CH3FH4= 	list(id="CH3FH4", 	ic=1.63,	comp="cell",bc=FALSE),
CHOFH4=		list(id="CHOFH4", 	ic=1.00,	comp="cell",bc=FALSE),
FFH2= 		list(id="FFH2", 	ic=0.000332,	comp="cell",bc=FALSE),
HCHO=		list(id="HCHO", 	ic=0.0074,	comp="cell",bc=FALSE),
FGAR= 		list(id="FGAR", 	ic=16.49,	comp="cell",bc=FALSE),
AICAR=		list(id="AICAR", 	ic=3.695,	comp="cell",bc=FALSE),

MTX1=		list(id="MTX1",	 	ic=0,		comp="cell",bc=FALSE),
MTX2=		list(id="MTX2",	 	ic=0,		comp="cell",bc=FALSE),
MTX3=		list(id="MTX3",	 	ic=0,		comp="cell",bc=FALSE),
MTX4=		list(id="MTX4",	 	ic=0,		comp="cell",bc=FALSE),
MTX5=		list(id="MTX5",	 	ic=0,		comp="cell",bc=FALSE),

MTX1b=		list(id="MTX1b", 	ic=0,		comp="cell",bc=FALSE),  # MTX bound to DHFR
MTX2b=		list(id="MTX2b", 	ic=0,		comp="cell",bc=FALSE),
MTX3b=		list(id="MTX3b", 	ic=0,		comp="cell",bc=FALSE),
MTX4b=		list(id="MTX4b", 	ic=0,		comp="cell",bc=FALSE),
MTX5b=		list(id="MTX5b", 	ic=0,		comp="cell",bc=FALSE),

EMTX=		list(id="EMTX",	 	ic=0,		comp="ext", bc=TRUE),
dUMP=		list(id="dUMP",		ic=20.76,	comp="cell",bc=TRUE), 
GAR=		list(id="GAR",		ic=689.6,	comp="cell",bc=TRUE), 
serine=		list(id="serine",	ic=123.3,	comp="cell",bc=TRUE), 
formate=	list(id="formate",	ic=500,		comp="cell",bc=TRUE), 
ATP=		list(id="ATP",		ic=2980,	comp="cell",bc=TRUE), 
glutamine=	list(id="glutamine",	ic=7170,	comp="cell",bc=TRUE),
glycine=	list(id="glycine",	ic=1600,	comp="cell",bc=TRUE), 
NADP=		list(id="NADP",		ic=6.73,	comp="cell",bc=TRUE),
NADPH=		list(id="NADPH",	ic=294,		comp="cell",bc=TRUE),
homocysteine=	list(id="homocysteine",	ic=10,		comp="cell",bc=TRUE)
),

parameters=list(Keq=0.32), # attach this list to make it global

rules=list(
list(id="FH2.binding",output=c("FH2b"),inputs=c("FH2f","DHFRf"),law=function(r){FH2f=r["FH2f"]; DHFRf=r["DHFRf"]; FH2f*DHFRf/Keq}),
#,list(id="TotalDHFR",output=c("DHFRtot"),inputs=c("FH2b","DHFRf","MTX1b","MTX2b","MTX3b","MTX4b","MTX5b"),law=function(r) sum(r))
# code above works fine in R alone, but for SBML export we need the rule's function explicitly, i.e.
list(id="TotalDHFR",output=c("DHFRtot"),inputs=c("FH2b","DHFRf","MTX1b","MTX2b","MTX3b","MTX4b","MTX5b"),law=function(r) 
{FH2b=r["FH2b"];DHFRf=r["DHFRf"];MTX1b=r["MTX1b"];MTX2b=r["MTX2b"];MTX3b=r["MTX3b"];MTX4b=r["MTX4b"];MTX5b=r["MTX5b"]
FH2b+DHFRf+MTX1b+MTX2b+MTX3b+MTX4b+MTX5b})
),
            
rxns=list(
list(	id="SHMT",		rever=FALSE,   # v1
	reacts=c("FH4","serine"),
	mods=c(NULL),
	prods =c("CH2FH4"),
	params=c(Vm=18330, Km1=1.7, Km2=210),
	law   = function(r,p) 
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"]
	FH4=r["FH4"];serine=r["serine"]
	Vm*( (serine/Km2)/(1+(serine/Km2)) ) * (FH4/Km1)/(1+(FH4/Km1)) }),

list(	id="SHMTr",		rever=FALSE,  # v1r
	reacts=c("CH2FH4"),
	mods=c("glycine"),
	prods =c("FH4"),
	params=c(Vm=12.2e6, Km1=3200, Km2=1e4),
	law   = function(r,p)
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"]
	CH2FH4=r["CH2FH4"];glycine=r["glycine"]
	Vm*( (glycine/Km2)/(1+(glycine/Km2)) ) * (CH2FH4/Km1)/(1+(CH2FH4/Km1)) }),

list(	id="HCHOtoCH2FH4",		rever=FALSE,  # v2f
	reacts=c("FH4","HCHO"),
	mods=c(NULL),
	prods =c("CH2FH4"),
	params=c(hp=23.2),
	law   = function(r,p)
	{hp=p["hp"]
	FH4=r["FH4"];HCHO=r["HCHO"]
 	hp*FH4*HCHO} ),
	
list(	id="CH2FH4toHCHO",		rever=FALSE,  # v2r
	reacts=c("CH2FH4"),
	mods=c(NULL),
	prods =c("FH4","HCHO"),
	params=c(hl=0.3),
	law   = function(r,p)
	{hl=p["hl"]
	CH2FH4=r["CH2FH4"]
	hl*CH2FH4} ),

list(	id="MTHFR",	rever=FALSE,   # v3
	reacts=c("CH2FH4","NADPH"),
	mods=c("FH2f","MTX1","MTX2","MTX3","MTX4","MTX5"),
	prods =c("CH3FH4"),
	params=c(Vm=224.8,Km1=50,Km2=50,Ki1=.40,Ki21=59,Ki22=21.3,Ki23=7.68,Ki24=2.77,Ki25=1.00), # Vm=k2
	law   = function(r,p) 
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Ki1=p["Ki1"];Ki21=p["Ki21"];Ki22=p["Ki22"];Ki23=p["Ki23"];Ki24=p["Ki24"];Ki25=p["Ki25"];
	CH2FH4=r["CH2FH4"];NADPH=r["NADPH"];FH2f=r["FH2f"];MTX1=r["MTX1"];MTX2=r["MTX2"];MTX3=r["MTX3"];MTX4=r["MTX4"];MTX5=r["MTX5"]
	Vm*CH2FH4*NADPH/(NADPH*CH2FH4+CH2FH4*Km2+(NADPH+Km2)*Km1*(1+(MTX1/Ki21)+(MTX2/Ki22)+(MTX3/Ki23)+(MTX4/Ki24)+(MTX5/Ki25)+FH2f/Ki1))}),

list(	id="MTR",	rever=FALSE,   # v4   # methionine synthase
	reacts=c("CH3FH4","homocysteine"),
	mods=c(NULL),
	prods =c("FH4"),
	params=c(Vm=22600, Km1=125, Km2=2900), 
	law   = function(r,p) 
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"]
	CH3FH4=r["CH3FH4"];homocysteine=r["homocysteine"]
	Vm*( (homocysteine/Km2)/(1+(homocysteine/Km2)) ) * (CH3FH4/Km1)/(1+(CH3FH4/Km1)) }),

list(	id="HCOOHtoCHOFH4",	rever=FALSE,  # v5
#	r5=FTS[4]/((1+FTS[1]/FH4)*(1+FTS[2]/ATP)*(1+FTS[3]/formate));
#	FTS=c(230, 56, 1600, 3600);
	reacts=c("FH4","formate","ATP"),
	mods=c(NULL),
	prods=c("CHOFH4"),  # = out of sustem
	params=c(Vm=3600, Km1=230,Km2=56,Km3=1600),
	law   = function(r,p) 
 	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Km3=p["Km3"]
	ATP=r["ATP"];formate=r["formate"];FH4=r["FH4"]
	Vm/((1+Km1/FH4)*(1+Km2/ATP)*(1+Km3/formate)) }),

list(	id="GARFT",	rever=FALSE,   # v6
	reacts=c("CHOFH4","GAR"),
	mods=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
	prods =c("FGAR","FH4"),
	params=c(Vm=4126,Km1=4.9,Km2=52,Ki1=5,Ki1f=1,Ki21=84,Ki22=60,Ki23=43,Ki24=31,Ki25=22), # Vm=k2
	law   = function(r,p) 
{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Ki1=p["Ki1"];Ki1f=p["Ki1f"];Ki21=p["Ki21"];Ki22=p["Ki22"];Ki23=p["Ki23"];Ki24=p["Ki24"];Ki25=p["Ki25"];
CHOFH4=r["CHOFH4"];GAR=r["GAR"];NADPH=r["NADPH"];FH2f=r["FH2f"];FFH2=r["FFH2"];MTX1=r["MTX1"];MTX2=r["MTX2"];MTX3=r["MTX3"];MTX4=r["MTX4"];MTX5=r["MTX5"]
Vm*CHOFH4*GAR/(GAR*CHOFH4+CHOFH4*Km2+(GAR+Km2)*Km1*(1+(MTX1/Ki21)+(MTX2/Ki22)+(MTX3/Ki23)+(MTX4/Ki24)+(MTX5/Ki25)+FH2f/Ki1+FFH2/Ki1f))}),

list(	id="ATIC7",	rever=FALSE,   # v7
	reacts=c("CHOFH4","AICAR"),
	mods=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
	prods =c("FH4"),  # + "FAICAR"
	params=c(Vm=31675,Km1=5.5,Km2=24,Ki1=2.89,Ki1f=5.3,Ki21=40,Ki22=31.5,Ki23=2.33,Ki24=3.61,Ki25=5.89), 
	law   = function(r,p) 
{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Ki1=p["Ki1"];Ki1f=p["Ki1f"];Ki21=p["Ki21"];Ki22=p["Ki22"];Ki23=p["Ki23"];Ki24=p["Ki24"];Ki25=p["Ki25"];
CHOFH4=r["CHOFH4"];AICAR=r["AICAR"];NADPH=r["NADPH"];FH2f=r["FH2f"];FFH2=r["FFH2"];MTX1=r["MTX1"];MTX2=r["MTX2"];MTX3=r["MTX3"];MTX4=r["MTX4"];MTX5=r["MTX5"]
Vm*CHOFH4*AICAR/(AICAR*CHOFH4+CHOFH4*Km2+(AICAR+Km2)*Km1*(1+(MTX1/Ki21)+(MTX2/Ki22)+(MTX3/Ki23)+(MTX4/Ki24)+(MTX5/Ki25)+FH2f/Ki1+FFH2/Ki1f))}),

list(	id="MTHFD",	rever=FALSE,   # v8
	reacts=c("CH2FH4","NADP"),
	mods=c(NULL),
	prods =c("CHOFH4"),
	params=c(Vm=68500, Km1=3, Km2=21.8), 
	law   = function(r,p) 
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"]
	CH2FH4=r["CH2FH4"];NADP=r["NADP"]
	Vm*((CH2FH4/Km1)/(1+(CH2FH4/Km1))) * ((NADP/Km2)/(1+(NADP/Km2)) )  }),

list(	id="TYMS",	rever=FALSE,   # v9
	reacts=c("CH2FH4","dUMP"),
	mods=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
	prods =c("FH2f"),
	params=c(Vm=58,Km1=2.5,Km2=1.8,Ki1=3,Ki1f=1.6,Ki21=13,Ki22=.08,Ki23=.07,Ki24=.065,Ki25=0.047), 
	law   = function(r,p) 
{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Ki1=p["Ki1"];Ki1f=p["Ki1f"];Ki21=p["Ki21"];Ki22=p["Ki22"];Ki23=p["Ki23"];Ki24=p["Ki24"];Ki25=p["Ki25"];
CH2FH4=r["CH2FH4"];dUMP=r["dUMP"];FH2f=r["FH2f"];FFH2=r["FFH2"];MTX1=r["MTX1"];MTX2=r["MTX2"];MTX3=r["MTX3"];MTX4=r["MTX4"];MTX5=r["MTX5"]
Vm*CH2FH4*dUMP/(dUMP*CH2FH4*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1) + 
Km1*dUMP*( (FFH2/Ki1f)*(MTX1/Ki21) + (1+FFH2/Ki1f)* (1+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1)) +
Km1*Km2*(1+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1))}),

list(	id="DHFReductase",	rever=FALSE,   # v10
	reacts=c("FH2f"),
	mods=c("FH2b"), 
	prods =c("FH4"),
	params=c(kter=2109.4), 
	law   = function(r,p) 
	{kter=p["kter"];FH2b=r["FH2b"]
	kter*FH2b}),

list(	id="FFH2syn",		rever=FALSE,  # v11
	reacts=c("FH2f"),
	mods=c(NULL),
	prods =c("FFH2"),
	params=c(Vm=65),
	law   = function(r,p)
	{Vm=p["Vm"]
	FH2f=r["FH2f"]
	Vm*FH2f} ),

list(	id="ATIC12",	rever=FALSE,   # v12
	reacts=c("FFH2","AICAR"),
	mods=c("FH2f","MTX1","MTX2","MTX3","MTX4","MTX5"),
	prods =c("FH2f"),  #  and FAICAR
	params=c(Vm=9503,Km1=5.3,Km2=24,Ki1=2.89,Ki1f=5.5,Ki21=40,Ki22=31.5,Ki23=2.33,Ki24=3.61,Ki25=5.89), 
	law   = function(r,p) 
{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"];Ki1=p["Ki1"];Ki1f=p["Ki1f"];Ki21=p["Ki21"];Ki22=p["Ki22"];Ki23=p["Ki23"];Ki24=p["Ki24"];Ki25=p["Ki25"];
AICAR=r["AICAR"];FH2f=r["FH2f"];FFH2=r["FFH2"];MTX1=r["MTX1"];MTX2=r["MTX2"];MTX3=r["MTX3"];MTX4=r["MTX4"];MTX5=r["MTX5"]
Vm*FFH2*AICAR/(AICAR*FFH2+FFH2*Km2+(AICAR+Km2)*Km1*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1+FFH2/Ki1f))}),

list(	id="AICARsyn",		rever=FALSE,  # v13  =15th reaction!
	reacts=c("FGAR"),
	mods=c("glutamine"),
	prods =c("AICAR"),
	params=c(Vm=4656,Km1=100,Km2=100),
	law   = function(r,p)
	{Vm=p["Vm"];Km1=p["Km1"];Km2=p["Km2"]
	FGAR=r["FGAR"];glutamine=r["glutamine"]
	Vm * ( (glutamine/Km1)/(1+(glutamine/Km1)) ) * ( (FGAR/Km2)/(1+(FGAR/Km2)) ) }),

#   End of reactions explicitly shown in Morrison's graph

list(	id="FPGS12",		rever=FALSE,  
	reacts=c("MTX1"),
	mods=c(NULL),
	prods =c("MTX2"),
	params=c(Vm=0.129),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX1=r["MTX1"]
	Vm * MTX1 }),
	
list(	id="FPGS23",		rever=FALSE,  
	reacts=c("MTX2"),
	mods=c(NULL),
	prods =c("MTX3"),
	params=c(Vm=0.369),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2=r["MTX2"]
	Vm * MTX2 }),
	
list(	id="FPGS34",		rever=FALSE,  
	reacts=c("MTX3"),
	mods=c(NULL),
	prods =c("MTX4"),
	params=c(Vm=0.118),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3=r["MTX3"]
	Vm * MTX3 }),

list(	id="FPGS45",		rever=FALSE,     # v19
	reacts=c("MTX4"),
	mods=c(NULL),
	prods =c("MTX5"),
	params=c(Vm=0.185),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4=r["MTX4"]
	Vm * MTX4 }),

list(	id="GGH21",		rever=FALSE,   # v20
	reacts=c("MTX2"),
	mods=c(NULL),
	prods =c("MTX1"),
	params=c(Vm=0.195),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2=r["MTX2"]
	Vm * MTX2 }),
	
list(	id="GGH32",		rever=FALSE,  
	reacts=c("MTX3"),
	mods=c(NULL),
	prods =c("MTX2"),
	params=c(Vm=0.025),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3=r["MTX3"]
	Vm * MTX3 }),
	
list(	id="GGH43",		rever=FALSE,  
	reacts=c("MTX4"),
	mods=c(NULL),
	prods =c("MTX3"),
	params=c(Vm=0.031),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4=r["MTX4"]
	Vm * MTX4 }),

list(	id="GGH54",		rever=FALSE,  
	reacts=c("MTX5"),
	mods=c(NULL),
	prods =c("MTX4"),
	params=c(Vm=0.191),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX5=r["MTX5"]
	Vm * MTX5 }),

list(	id="RFC",		rever=FALSE,    
	reacts=c("EMTX"),
	mods=c(NULL),
	prods =c("MTX1"),
	params=c(Vm=82.2,Km=8.2),
	law   = function(r,p)
	{Vm=p["Vm"];Km=p["Km"]
	EMTX=r["EMTX"]
	Vm * EMTX/(Km+EMTX) }),

list(	id="MTX1export",		rever=FALSE,  # v25
	reacts=c("MTX1"),
	mods=c(NULL),
	prods =c(NULL), 
	params=c(Vm=4.65),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX1=r["MTX1"]
	Vm * MTX1 }),

list(	id="MTX2export",		rever=FALSE,  
	reacts=c("MTX2"),  # NOTE: MTX1 efflux is so strong that
	mods=c(NULL),      # GGH21 overwhelms this process  
	prods =c(NULL),  #  and approximates it with
	params=c(Vm=0.0), # zero!!!!!
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2=r["MTX2"]
	Vm * MTX2 }),

list(	id="MTX3export",		rever=FALSE,  
	reacts=c("MTX3"),
	mods=c(NULL),
	prods =c(NULL), 
	params=c(Vm=0.063),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3=r["MTX3"]
	Vm * MTX3 }),

list(	id="MTX4export",		rever=FALSE,  
	reacts=c("MTX4"),
	mods=c(NULL),
	prods =c(NULL),  
	params=c(Vm=0.063),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4=r["MTX4"]
	Vm * MTX4 }),
	
list(	id="MTX5export",		rever=FALSE,  
	reacts=c("MTX5"),
	mods=c(NULL),
	prods =c(NULL), 
	params=c(Vm=0.063),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX5=r["MTX5"]
	Vm * MTX5 }),

# MTX binding to DHFR reactions
list(	id="MTX1on",		rever=FALSE,           #v30  
	reacts=c("MTX1","DHFRf"),
	mods=c(NULL),
	prods =c("MTX1b"),
	params=c(Vm=.231e5),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX1=r["MTX1"];DHFRf=r["DHFRf"]
	Vm * DHFRf* MTX1 }),

list(	id="MTX2on",		rever=FALSE,  
	reacts=c("MTX2","DHFRf"),
	mods=c(NULL),
	prods =c("MTX2b"),
	params=c(Vm=.443e5),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2=r["MTX2"];DHFRf=r["DHFRf"]
	Vm * DHFRf* MTX2 }),

list(	id="MTX3on",		rever=FALSE,  
	reacts=c("MTX3","DHFRf"),
	mods=c(NULL),
	prods =c("MTX3b"),
	params=c(Vm=.851e5),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3=r["MTX3"];DHFRf=r["DHFRf"]
	Vm * DHFRf* MTX3 }),

list(	id="MTX4on",		rever=FALSE,  
	reacts=c("MTX4","DHFRf"),
	mods=c(NULL),
	prods =c("MTX4b"),
	params=c(Vm=1.63e5),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4=r["MTX4"];DHFRf=r["DHFRf"]
	Vm * DHFRf* MTX4 }),

list(	id="MTX5on",		rever=FALSE,  
	reacts=c("MTX5","DHFRf"),
	mods=c(NULL),
	prods =c("MTX5b"),
	params=c(Vm=3.14e5),
	law   = function(r,p) 
	{Vm=p["Vm"]; 
	MTX5=r["MTX5"];DHFRf=r["DHFRf"]
	Vm * DHFRf* MTX5 }),

list(	id="MTX1off",		rever=FALSE,  # v35
	reacts=c("MTX1b"),
	mods=c(NULL),
	prods =c("MTX1","DHFRf"),
	params=c(Vm=.42),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX1b=r["MTX1b"]
	Vm * MTX1b }),
	
list(	id="MTX2off",		rever=FALSE,  
	reacts=c("MTX2b"),
	mods=c(NULL),
	prods =c("MTX2","DHFRf"),
	params=c(Vm=.42),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2b=r["MTX2b"]
	Vm * MTX2b }),

list(	id="MTX3off",		rever=FALSE,  
	reacts=c("MTX3b"),
	mods=c(NULL),
	prods =c("MTX3","DHFRf"),
	params=c(Vm=.42),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3b=r["MTX3b"]
	Vm * MTX3b }),

list(	id="MTX4off",		rever=FALSE,  
	reacts=c("MTX4b"),
	mods=c(NULL),
	prods =c("MTX4","DHFRf"),
	params=c(Vm=.42),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4b=r["MTX4b"]
	Vm * MTX4b }),

list(	id="MTX5off",		rever=FALSE,  
	reacts=c("MTX5b"),
	mods=c(NULL),
	prods =c("MTX5","DHFRf"),
	params=c(Vm=.42),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX5b=r["MTX5b"]
	Vm * MTX5b }),

#   DHFR protein synthesis 
list(	id="DHFRfsyn",	rever=FALSE,   # v40
	reacts=c(NULL),
	mods=c("EMTX"),
	prods =c("DHFRf"),
	params=c(k0=.0192,k1=.04416), # figure these backwards to get SS concentration right with kdeg=.03
	law   = function(r,p)         # and also to get 1uM MTX to ramp it by a factor of 3.3
	{k0=p["k0"];k1=p["k1"]
	DHFRf=r["DHFRf"];
	EMTX=r["EMTX"];
	k0+k1*EMTX}),

#   These are the DHFR protein degradation reactions: independent of binding.
list(	id="DHFRdeg",		rever=FALSE,  
	reacts=c("DHFRf"),# can't pull mass out of BC therefore pull it form free DHFR
	mods=c("FH2b"),
	prods =c(NULL),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	DHFRf=r["DHFRf"];FH2b=r["FH2b"]
	Vm * (DHFRf+FH2b) }),  

list(	id="FH2bdeg",		rever=FALSE,  
	reacts=c(NULL),
	mods=c("FH2b"),
	prods =c("FH2f"),  # create back the FH2 missing from the protein degradation above
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	FH2b=r["FH2b"]
	Vm * FH2b }),  
	
list(	id="MTX1deg",		rever=FALSE,  
	reacts=c("MTX1b"),
	mods=c(NULL),
	prods =c("MTX1"),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX1b=r["MTX1b"]
	Vm * MTX1b }),
	
list(	id="MTX2deg",		rever=FALSE,  
	reacts=c("MTX2b"),
	mods=c(NULL),
	prods =c("MTX2"),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX2b=r["MTX2b"]
	Vm * MTX2b }),
	
list(	id="MTX3deg",		rever=FALSE,  # V45
	reacts=c("MTX3b"),
	mods=c(NULL),
	prods =c("MTX3"),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX3b=r["MTX3b"]
	Vm * MTX3b }),
	
list(	id="MTX4deg",		rever=FALSE,  
	reacts=c("MTX4b"),
	mods=c(NULL),
	prods =c("MTX4"),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX4b=r["MTX4b"]
	Vm * MTX4b }),
	
list(	id="MTX5deg",		rever=FALSE,  #47
	reacts=c("MTX5b"),
	mods=c(NULL),
	prods =c("MTX5"),
	params=c(Vm=.03),
	law   = function(r,p)
	{Vm=p["Vm"]
	MTX5b=r["MTX5b"]
	Vm * MTX5b })	
),  # end rxn list
 
units=c("micromolar","hour") 
) # end model list
 
 