# Curto's purine model with reaction names modified to lower case old notation
model=list(
notes=c("This is a purine metabolism model that is geared toward studies of gout.", 
"The model is fully described in Curto et al., MBSC 151 (1998) pp 1-49",
"The model uses Generalized Mass Action (GMA;i.e. power law) descriptions of reaction rate laws.",
"Such descriptions are local approximations that assume independent substrate binding."),

# vden= 2.39;  this in umole/min/KG which is 2.4*60=144 uM/h (if we let each Kg be a liter of water)
# Morrison and Allegra (JBC, 1989) have vden at 650 uM/h (model) and 415 (exp)

comps=list(list(id= "cell",vol=1)),

species=list( # IC's have already been run out to the system's steady state. 
PRPP	=list(id="PRPP", ic=5.01742447344584,	comp="cell",	 bc=FALSE),	 
IMP	=list(id="IMP",	 ic=98.2633505727176,	comp="cell",	 bc=FALSE),	 
SAMP 	=list(id="SAMP", ic=.198189051721450,	comp="cell",	 bc=FALSE),	 
ATP	=list(id="ATP",	 ic=2475.35163831273,	comp="cell",	 bc=FALSE),	 
SAM	=list(id="SAM",	 ic=3.99186964914143,	comp="cell",	 bc=FALSE),	 
Ade	=list(id="Ade",	 ic=0.984730215691558,	comp="cell",	 bc=FALSE),	 
XMP	=list(id="XMP",	 ic=24.7929534893554,	comp="cell",	 bc=FALSE),	 
GTP	=list(id="GTP",	 ic=410.222754152001,	comp="cell",	 bc=FALSE),	 
dATP 	=list(id="dATP", ic=6.01412837487651,	comp="cell",	 bc=FALSE),	 
dGTP 	=list(id="dGTP", ic=3.02581107154342,	comp="cell",	 bc=FALSE),	 
RNA	=list(id="RNA",	 ic=28680.4684823855,	comp="cell",	 bc=FALSE),	 
DNA	=list(id="DNA",	 ic=5179.3359175299,	comp="cell",	 bc=FALSE),	 
HX 	=list(id="HX",	 ic=9.51785484271038,	comp="cell",	 bc=FALSE),	 
Xa 	=list(id="Xa",	 ic=5.05940999539495,	comp="cell",	 bc=FALSE),	 
Gua	=list(id="Gua",	 ic=5.50638133716038,	comp="cell",	 bc=FALSE),	 
UA 	=list(id="UA",	 ic=100.293296074217,	comp="cell",	 bc=FALSE),	 
R5P	=list(id="R5P",	 ic=18,	 		comp="cell",	 bc=TRUE),	 
Pi 	=list(id="Pi",	 ic=1400,		comp="cell",	 bc=TRUE)	 
),
      
rxns=list(
list(	id="ada",		rever=FALSE,   # v1
	reacts=c("ATP"),
	prods =c("HX"),
	params=c(aada =0.001062, fada4 =0.97),
	law   = function(r,p) 
	{aada=p["aada"];fada4=p["fada4"]
	ATP=r["ATP"]
	aada*ATP^fada4 }),

list(	id="ade",		rever=FALSE,   # v2
	reacts=c("Ade"),
	params=c(aade =0.01, fade6 =0.55),
	law   = function(r,p) 
	{aade=p["aade"];fade6=p["fade6"]
	Ade=r["Ade"]
	aade   * Ade^fade6 }),

list(	id="adna",		rever=FALSE,   # v3
	reacts=c("dATP"),
	mods=c("dGTP"),
	prods =c("DNA"),
	params=c(aadna=3.2789, fdnap9 =0.42, fdnap10 =0.33),
	law   = function(r,p) 
	{aadna=p["aadna"];fdnap9=p["fdnap9"];fdnap10=p["fdnap10"]
	dATP=r["dATP"];dGTP=r["dGTP"]
	aadna  * dATP^fdnap9 * dGTP^fdnap10}),

list(	id="adrnr",		rever=FALSE,   # v4
	reacts=c("ATP"),
	mods=c("dGTP","dATP"),
	prods =c("dATP"),
	params=c(aadrnr =0.0602, fadrnr4 =0.1, fadrnr9= -0.3, fadrnr10=0.87),
	law   = function(r,p) 
	{aadrnr=p["aadrnr"];fadrnr4=p["fadrnr4"];fadrnr9=p["fadrnr9"];fadrnr10=p["fadrnr10"]
	ATP=r["ATP"];dATP=r["dATP"];dGTP=r["dGTP"]
	aadrnr * ATP^fadrnr4 * dATP^fadrnr9 * dGTP^fadrnr10}),

list(	id="ampd",		rever=FALSE,   # v5
	reacts=c("ATP"),
	mods=c("GTP","Pi"),
	prods =c("IMP"),
	params=c(aampd =0.02688, fampd4 =0.8, fampd8 =-0.03,fampd18 =-0.1),
	law   = function(r,p) 
	{aampd=p["aampd"];fampd4=p["fampd4"];fampd8=p["fampd8"];fampd18=p["fampd18"]
	ATP=r["ATP"];GTP=r["GTP"];Pi=r["Pi"]
	aampd*ATP^fampd4*GTP^fampd8*Pi^fampd18}),

list(	id="aprt",		rever=FALSE,   # v6
	reacts=c("PRPP","Ade"),
	mods=c("ATP"),
	prods =c("ATP"),
	params=c(aaprt =233.8, faprt1 =0.5, faprt4 =-0.8,faprt6 =0.75),
	law   = function(r,p) 
	{aaprt=p["aaprt"];faprt1=p["faprt1"];faprt4=p["faprt4"];faprt6=p["faprt6"]
	ATP=r["ATP"];PRPP=r["PRPP"];Ade=r["Ade"]
	aaprt*PRPP^faprt1*ATP^faprt4*Ade^faprt6}),

list(	id="arna",		rever=FALSE,   # v33
	reacts=c("ATP"),
	mods=c("GTP"),
	prods =c("RNA"),
	params=c(aarna =614.5,frnap4 =0.05, frnap8 =0.13),
	law   = function(r,p) 
	{aarna=p["aarna"];frnap4=p["frnap4"];frnap8=p["frnap8"]
	ATP=r["ATP"];GTP=r["GTP"]
	aarna *  ATP^frnap4* GTP^frnap8}),

list(	id="asuc",		rever=FALSE,   # v7
	reacts=c("IMP"),
	mods=c("ATP","GTP","Pi"),
	prods =c("SAMP"),
	params=c(aasuc =3.5932, fasuc2 =0.4, fasuc4 =-.24,fasuc8 =0.2, fasuc18 =-.05),
	law   = function(r,p) 
	{aasuc=p["aasuc"];fasuc2=p["fasuc2"];fasuc4=p["fasuc4"];fasuc8=p["fasuc8"];fasuc18=p["fasuc18"]
	IMP=r["IMP"];ATP=r["ATP"];GTP=r["GTP"];Pi=r["Pi"]
	aasuc*IMP^fasuc2*ATP^fasuc4*GTP^fasuc8*Pi^fasuc18}),

list(	id="asli",		rever=FALSE,   # v8
	reacts=c("SAMP"),
	mods=c("ATP"),
	prods =c("ATP"),
	params=c(aasli =66544, fasli3 =0.99, fasli4 =-.95),
	law   = function(r,p) 
	{aasli=p["aasli"];fasli3=p["fasli3"];fasli4=p["fasli4"]
	SAMP=r["SAMP"];ATP=r["ATP"]
	aasli*SAMP^fasli3*ATP^fasli4}),

list(	id="dada",		rever=FALSE,   # v9
	reacts=c("dATP"),
	prods =c("HX"),
	params=c(adada =0.03333, fdada9 =1),
	law   = function(r,p) 
	{adada=p["adada"];fdada9=p["fdada9"]
	dATP=r["dATP"]
	adada*dATP^fdada9}),

list(	id="den",		rever=FALSE,   # v10
	reacts=c("PRPP"),
	mods=c("dGTP","IMP","ATP","GTP","Pi"),
	prods =c("IMP"),
	params=c(aden =5.2728, fden1 =2, fden2 =-.06, fden4 =-.25, fden8 =-.2, fden18 =-.08),
	law   = function(r,p) 
	{aden=p["aden"];fden1=p["fden1"];fden2=p["fden2"];fden4=p["fden4"];fden8=p["fden8"];fden18=p["fden18"]
	PRPP=r["PRPP"];IMP=r["IMP"];ATP=r["ATP"];GTP=r["GTP"];Pi=r["Pi"]
	aden*PRPP^fden1*IMP^fden2*ATP^fden4*GTP^fden8*Pi^fden18}),

list(	id="dgnuc",		rever=FALSE,   # v11
	reacts=c("dGTP"),
	prods =c("Gua"),
	params=c(adgnuc=0.03333, fdgnuc10=1),
	law   = function(r,p) 
	{adgnuc=p["adgnuc"];fdgnuc10=p["fdgnuc10"]
	dGTP=r["dGTP"]
	adgnuc*dGTP^fdgnuc10}),

list(	id="dnaa",		rever=FALSE,   # v12
	reacts=c("DNA"),
	prods =c("dATP"),
	params=c(adnaa =0.001938, fdnan12=1),
	law   = function(r,p) 
	{adnaa=p["adnaa"];fdnan12=p["fdnan12"]
	DNA=r["DNA"]
	adnaa*DNA^fdnan12}),

list(	id="dnag",		rever=FALSE,   # v13
	reacts=c("DNA"),
	prods =c("dGTP"),
	params=c(adnag =0.001318, fdnan12=1),
	law   = function(r,p) 
	{adnag=p["adnag"];fdnan12=p["fdnan12"]
	DNA=r["DNA"]
	adnag*DNA^fdnan12}),

list(	id="gdna",		rever=FALSE,   # v14
	reacts=c("dGTP"),
	mods=c("dATP"),
	prods =c("DNA"),
	params=c(agdna=2.2296, fdnap9 =0.42, fdnap10 =0.33),
	law   = function(r,p) 
	{agdna=p["agdna"];fdnap9=p["fdnap9"];fdnap10=p["fdnap10"]
	dATP=r["dATP"];dGTP=r["dGTP"]
	agdna  * dATP^fdnap9 * dGTP^fdnap10}),

list(	id="gdrnr",		rever=FALSE,   # v15
	reacts=c("GTP"),
	mods=c("dATP","dGTP"),
	prods =c("dGTP"),
	params=c(agdrnr =0.1199, fgdrnr8 =0.4, fgdrnr9 =-1.2,fgdrnr10=-.39),
	law   = function(r,p) 
	{agdrnr=p["agdrnr"];fgdrnr8=p["fgdrnr8"];fgdrnr9=p["fgdrnr9"];fgdrnr10=p["fgdrnr10"]
	dATP=r["dATP"];GTP=r["GTP"];dGTP=r["dGTP"]
	agdrnr * GTP^fgdrnr8 * dATP^fgdrnr9 * dGTP^fgdrnr10}),

list(	id="gmpr",		rever=FALSE,   # v16
	reacts=c("GTP"),
	mods=c("XMP","ATP","IMP"),
	prods =c("IMP"),
	params=c(agmpr =0.3005, fgmpr2 =-.15, fgmpr4 =-.07,fgmpr7 =-.76,fgmpr8 =0.7),
	law   = function(r,p) 
	{agmpr=p["agmpr"];fgmpr2=p["fgmpr2"];fgmpr4=p["fgmpr4"];fgmpr7=p["fgmpr7"];fgmpr8=p["fgmpr8"]
	ATP=r["ATP"];XMP=r["XMP"];IMP=r["IMP"];GTP=r["GTP"]
	agmpr * IMP^fgmpr2 * ATP^fgmpr4 * XMP^fgmpr7 * GTP^fgmpr8}),

list(	id="gmps",		rever=FALSE,   # v17
	reacts=c("XMP"),
	mods=c("ATP"),
	prods =c("GTP"),
	params=c(agmps=0.3738, fgmps4 =0.12, fgmps7 =0.16),
	law   = function(r,p) 
	{agmps=p["agmps"];fgmps4=p["fgmps4"];fgmps7=p["fgmps7"]
	ATP=r["ATP"];XMP=r["XMP"]
	agmps * ATP^fgmps4 * XMP^fgmps7 }),

list(	id="gnuc",		rever=FALSE,   # v18
	reacts=c("GTP"),
	mods=c("Pi"),
	prods =c("Gua"),
	params=c(agnuc=0.2511, fgnuc8 =0.9, fgnuc18 =-.34),
	law   = function(r,p) 
	{agnuc=p["agnuc"];fgnuc8=p["fgnuc8"];fgnuc18=p["fgnuc18"]
	GTP=r["GTP"];Pi=r["Pi"]
	agnuc * GTP^fgnuc8 * Pi^fgnuc18 }),

list(	id="gprt",		rever=FALSE,   # v19
	reacts=c("Gua","PRPP"),
	mods=c("GTP"),
	prods =c("GTP"),
	params=c(agprt =361.69, fgprt1 =1.2, fgprt8 =-1.2,fgprt15 =0.42),
	law   = function(r,p) 
	{agprt=p["agprt"];fgprt1=p["fgprt1"];fgprt8=p["fgprt8"];fgprt15=p["fgprt15"]
	PRPP=r["PRPP"];GTP=r["GTP"];Gua=r["Gua"]
	agprt * PRPP^fgprt1* GTP^fgprt8 * Gua^fgprt15 }),

list(	id="grna",		rever=FALSE,   # v32
	reacts=c("GTP"),
	mods=c("ATP"),
	prods =c("RNA"),
	params=c(agrna =409.6,frnap4 =0.05, frnap8 =0.13),
	law   = function(r,p) 
	{agrna=p["agrna"];frnap4=p["frnap4"];frnap8=p["frnap8"]
	ATP=r["ATP"];GTP=r["GTP"]
	agrna *  ATP^frnap4* GTP^frnap8}),

list(	id="gua",		rever=FALSE,   # v21
	reacts=c("Gua"),
	prods =c("Xa"),
	params=c(agua =0.4919, fgua15 =0.5),
	law   = function(r,p) 
	{agua=p["agua"];fgua15=p["fgua15"]
	Gua=r["Gua"]
	agua *  Gua^fgua15 }),

list(	id="hprt",		rever=FALSE,   # v20
	reacts=c("HX","PRPP"),
	mods=c("IMP"),
	prods =c("IMP"),
	params=c(ahprt =12.569, fhprt1 =1.1,fhprt2 =-.89, fhprt13 =0.48),
	law   = function(r,p) 
	{ahprt=p["ahprt"];fhprt1=p["fhprt1"];fhprt2=p["fhprt2"];fhprt13=p["fhprt13"]
	PRPP=r["PRPP"];IMP=r["IMP"];HX=r["HX"]
	ahprt * PRPP^fhprt1* IMP^fhprt2 * HX^fhprt13 }),

list(	id="hx",		rever=FALSE,   # v22
	reacts=c("HX"),
	mods=c(NULL),
	prods =c(NULL),
	params=c(ahx =0.003793, fhx13 =1.12),
	law   = function(r,p) 
	{ahx=p["ahx"];fhx13=p["fhx13"]
	HX=r["HX"]
	ahx *  HX^fhx13 }),

list(	id="hxd",		rever=FALSE,   # v23
	reacts=c("HX"),
	mods=c(NULL),
	prods =c("Xa"),
	params=c(ahxd =0.2754, fhxd13 =0.65),
	law   = function(r,p) 
	{ahxd=p["ahxd"];fhxd13=p["fhxd13"]
	HX=r["HX"]
	ahxd * HX^fhxd13 }),

list(	id="impd",		rever=FALSE,   # v24
	reacts=c("IMP"),
	mods=c("GTP","XMP"),
	prods =c("XMP"),
	params=c(aimpd =1.2823, fimpd2 =0.15,fimpd7 =-.09, fimpd8 =-.03),
	law   = function(r,p) 
	{aimpd=p["aimpd"];fimpd2=p["fimpd2"];fimpd7=p["fimpd7"];fimpd8=p["fimpd8"]
	IMP=r["IMP"];XMP=r["XMP"];GTP=r["GTP"]
	aimpd*IMP^fimpd2 * XMP^fimpd7 *  GTP^fimpd8  }),

list(	id="inuc",		rever=FALSE,   # v25
	reacts=c("IMP"),
	mods=c("Pi"),
	prods =c("HX"),
	params=c(ainuc =0.9135, finuc2 =0.8, finuc18 =-.36),
	law   = function(r,p) 
	{ainuc=p["ainuc"];finuc2=p["finuc2"];finuc18=p["finuc18"]
	IMP=r["IMP"];Pi=r["Pi"]
	ainuc *  IMP^finuc2 *  Pi^finuc18  }),

list(	id="mat",		rever=FALSE,   # v26
	reacts=c("ATP"),
	mods=c("SAM"),
	prods =c("SAM"),
	params=c(amat =7.2067,fmat4=0.2, fmat5=-.6), 
	law   = function(r,p) 
	{amat=p["amat"];fmat4=p["fmat4"];fmat5=p["fmat5"]
	SAM=r["SAM"];ATP=r["ATP"]
	amat *  ATP^fmat4 *  SAM^fmat5  }),

list(	id="polyam",		rever=FALSE,   # v27
	reacts=c("SAM"),
	prods =c("Ade"),
	params=c(apolyam=0.29, fpolyam5=0.9), 
	law   = function(r,p) 
	{apolyam=p["apolyam"];fpolyam5=p["fpolyam5"]
	SAM=r["SAM"]
	apolyam *  SAM^fpolyam5}),

list(	id="prpps",		rever=FALSE,   # v28
	reacts=c("R5P"),
	mods=c("ATP","GTP","Pi","PRPP"),
	prods =c("PRPP"),
	params=c(aprpps=0.9, fprpps1 =-.03, fprpps4 =-.45, fprpps8 =-.04, fprpps17=0.65, fprpps18 =0.7),
	law   = function(r,p) 
	{aprpps=p["aprpps"];fprpps1=p["fprpps1"];fprpps4=p["fprpps4"];fprpps8=p["fprpps8"];fprpps17=p["fprpps17"];fprpps18=p["fprpps18"]
	PRPP=r["PRPP"];ATP=r["ATP"];GTP=r["GTP"];R5P=r["R5P"];Pi=r["Pi"]
	aprpps * PRPP^fprpps1 * ATP^fprpps4 * GTP^fprpps8 * R5P^fprpps17 * Pi^fprpps18}),

list(	id="pyr",		rever=FALSE,   # v29
	reacts=c("PRPP"),
	params=c(apyr =1.2951, fpyr1 =1.27),
	law   = function(r,p) 
	{apyr=p["apyr"];fpyr1=p["fpyr1"]
	PRPP=r["PRPP"]
	apyr *  PRPP^fpyr1}),

list(	id="rnaa",		rever=FALSE,   # v30
	reacts=c("RNA"),
	prods =c("ATP"),
	params=c(arnaa =0.06923,frnan11 =1),
	law   = function(r,p) 
	{arnaa=p["arnaa"];frnan11=p["frnan11"]
	RNA=r["RNA"]
	arnaa*RNA^frnan11}),

list(	id="rnag",		rever=FALSE,   # v31
	reacts=c("RNA"),
	prods =c("GTP"),
	params=c(arnag =0.04615,frnan11 =1),
	law   = function(r,p) 
	{arnag=p["arnag"];frnan11=p["frnan11"]
	RNA=r["RNA"]
	arnag*RNA^frnan11}),

list(	id="trans",		rever=FALSE,   # v34
	reacts=c("SAM"),
	prods =c("ATP"),
	params=c(atrans =8.8539, ftrans5 =0.33),
	law   = function(r,p) 
	{atrans=p["atrans"];ftrans5=p["ftrans5"]
	SAM=r["SAM"]
	atrans *  SAM^ftrans5}),

list(	id="ua",		rever=FALSE,   # v35
	reacts=c("UA"),
	params=c(aua =0.00008744,fua16 =2.21),
	law   = function(r,p) 
	{aua=p["aua"];fua16=p["fua16"]
	UA=r["UA"]
	aua * UA^fua16}),

list(	id="x",		rever=FALSE,   # v36
	reacts=c("Xa"),
	params=c(ax =0.0012,fx14 =2.0),
	law   = function(r,p) 
	{ax=p["ax"];fx14=p["fx14"]
	Xa=r["Xa"]
	ax *  Xa^fx14}),

list(	id="xd",		rever=FALSE,   # v37
	reacts=c("Xa"),
	prods =c("UA"),
	params=c(axd =0.949,fxd14 =0.55),
	law   = function(r,p) 
	{axd=p["axd"];fxd14=p["fxd14"]
	Xa=r["Xa"]
	axd *  Xa^fxd14}) ), 
units=c("micromolar","minutes") 
) # end model list

