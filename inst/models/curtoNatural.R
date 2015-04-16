# This script, curtoNatural.R, shows how Curto's model
# might be implemented in R naturally (without knowledge of SBML).
library(deSolve)

# first define the powers f of the rate laws 
fada4 =0.97; fade6 =0.55; fadrnr4 =0.1; fadrnr9 =-0.3
fadrnr10=0.87; fampd4 =0.8; fampd8 =-0.03; fampd18 =-0.1
faprt1 =0.5; faprt4 =-0.8; faprt6 =0.75; fasuc2 =0.4
fasuc4 =-.24; fasuc8 =0.2; fasuc18 =-.05; fasli3 =0.99
fasli4 =-.95; fdada9 =1; fden1 =2; fden2 =-.06
fden4 =-.25; fden8 =-.2; fden18 =-.08; fdgnuc10=1
fdnan12=1; fdnap9 =0.42; fdnap10 =0.33; fgdrnr8 =0.4
fgdrnr9 =-1.2; fgdrnr10=-.39; fgmpr2 =-.15; fgmpr4 =-.07
fgmpr7 =-.76; fgmpr8 =0.7; fgmps4 =0.12; fgmps7 =0.16
fgnuc8 =0.9; fgnuc18 =-.34; fgprt1 =1.2; fgprt8 =-1.2
fgprt15 =0.42; fgua15 =0.5; fhprt1 =1.1; fhprt2 =-.89
fhprt13 =0.48; fhx13 =1.12; fhxd13 =0.65; fimpd2 =0.15
fimpd7 =-.09; fimpd8 =-.03; finuc2 =0.8; finuc18 =-.36
fmat4=0.2; fmat5=-.6; fpolyam5=0.9; fprpps1 =-.03
fprpps4 =-.45; fprpps8 =-.04; fprpps17=0.65; fprpps18 =0.7
fpyr1 =1.27; frnan11 =1; frnap4 =0.05; frnap8 =0.13
ftrans5 =0.33; fua16 =2.21; fx14 =2.0; fxd14 =0.55

# now define the rate law amplitude parameters
aada =0.001062; aade =0.01; aadna=3.2789; aadrnr =0.0602
aampd =0.02688; aaprt =233.8; aarna =614.5; aasuc =3.5932
aasli =66544; adada =0.03333; aden =5.2728; adgnuc=0.03333
adnaa =0.001938; adnag =0.001318; agdna=2.2296; agdrnr =0.1199
agmpr =0.3005; agmps=0.3738; agnuc=0.2511; agprt =361.69
agrna =409.6; agua =0.4919; ahprt =12.569; ahx =0.003793
ahxd =0.2754; aimpd =1.2823; ainuc =0.9135; amat =7.2067
apolyam=0.29; aprpps=0.9; apyr =1.2951; arnaa =0.06923
arnag =0.04615; atrans =8.8539; aua =0.00008744; ax =0.0012;axd =0.949

# now define the two boundary conditions
R5P=18; Pi=1400;

# yes, it would be better to stick all of the above in p and pass directly
# rather than globally as is being done here. 

# Now define the ODE right hand side
fpur <- function(t, X, p)
{ # first define the fluxes
  vada =aada * X[4]^fada4;
  vade =aade * X[6]^fade6;
  vadna =aadna * X[9]^fdnap9 * X[10]^fdnap10;
  vadrnr =aadrnr * X[4]^fadrnr4 * X[9]^fadrnr9 * X[10]^fadrnr10;
  vampd =aampd * X[4]^fampd4 * X[8]^fampd8 * Pi^fampd18;
  vaprt =aaprt * X[1]^faprt1 * X[4]^faprt4 * X[6]^faprt6;
  varna =aarna * X[4]^frnap4* X[8]^frnap8;
  vasuc =aasuc * X[2]^fasuc2 * X[4]^fasuc4 * X[8]^fasuc8 * Pi^fasuc18;
  vasli =aasli * X[3]^fasli3 * X[4]^fasli4;
  vdada =adada * X[9]^fdada9;
  vden =aden * X[1]^fden1 * X[2]^fden2 * X[4]^fden4 * X[8]^fden8 * Pi^fden18;
  vdgnuc =adgnuc * X[10]^fdgnuc10;
  vdnaa =adnaa * X[12]^fdnan12;
  vdnag =adnag * X[12]^fdnan12;
  vgdna =agdna * X[9]^fdnap9 * X[10]^fdnap10;
  vgdrnr =agdrnr * X[8]^fgdrnr8 * X[9]^fgdrnr9 * X[10]^fgdrnr10;
  vgmpr =agmpr * X[2]^fgmpr2 * X[4]^fgmpr4 * X[7]^fgmpr7 * X[8]^fgmpr8;
  vgmps =agmps * X[4]^fgmps4 * X[7]^fgmps7;
  vgnuc =agnuc * X[8]^fgnuc8 * Pi^fgnuc18;
  vgprt =agprt * X[1]^fgprt1 * X[8]^fgprt8* X[15]^fgprt15;
  vgrna =agrna * X[8]^frnap8* X[4]^frnap4;
  vgua =agua * X[15]^fgua15;
  vhprt =ahprt * X[1]^fhprt1 * X[2]^fhprt2 * X[13]^fhprt13;
  vhx =ahx * X[13]^fhx13;
  vhxd =ahxd * X[13]^fhxd13;
  vimpd =aimpd * X[2]^fimpd2 * X[7]^fimpd7 * X[8]^fimpd8;
  vinuc =ainuc * X[2]^finuc2 * Pi^finuc18;
  vmat =amat * X[4]^fmat4 * X[5]^fmat5;
  vpolyam=apolyam* X[5]^fpolyam5;
  vprpps =aprpps * X[1]^fprpps1 * X[4]^fprpps4 * X[8]^fprpps8 * R5P^fprpps17 * Pi^fprpps18;
  vpyr =apyr * X[1]^fpyr1;
  vrnaa =arnaa * X[11]^frnan11;
  vrnag =arnag * X[11]^frnan11;
  vtrans =atrans * X[5]^ftrans5;
  vua =aua * X[16]^fua16;
  vx =ax * X[14]^fx14;
  vxd =axd * X[14]^fxd14;
  # now define dX/dt (Xi prime) as fluxes into a node minus fluxes out
  X1p = vprpps-vgprt-vhprt-vaprt-vden-vpyr;
  X2p = vden+vgmpr+vhprt+vampd-vimpd-vasuc-vinuc;
  X3p = vasuc-vasli;
  X4p = vaprt+vasli+vtrans+vrnaa-vmat-vampd-varna-vadrnr-vada;
  X5p = vmat-vtrans-vpolyam;
  X6p = vpolyam-vaprt-vade;
  X7p = vimpd-vgmps;
  X8p = vgmps+vrnag+vgprt-vgmpr-vgrna-vgdrnr-vgnuc;
  X9p = vadrnr+vdnaa-vadna-vdada;
  X10p= vgdrnr+vdnag-vgdna-vdgnuc;
  X11p= varna+vgrna-vrnaa-vrnag;
  X12p= vgdna+vadna-vdnaa-vdnag;
  X13p= vdada+vada+vinuc-vhprt-vhxd-vhx;
  X14p= vhxd+vgua-vxd-vx;
  X15p= vgnuc+vdgnuc-vgua-vgprt;
  X16p= vxd-vua;
  XP = c(X1p, X2p, X3p, X4p, X5p, X6p, X7p, X8p, X9p, X10p, X11p, X12p, X13p, X14p, X15p, X16p);
  V=c(vada,vade,vadna,vadrnr,vampd,vaprt,varna,vasuc,vasli,vdada,vden,vdgnuc,vdnaa,
      vdnag,vgdna,vgdrnr,vgmpr,vgmps,vgnuc,vgprt,vgrna,vgua,vhprt,vhx,vhxd,vimpd,
      vinuc,vmat,vpolyam, vprpps,vpyr,vrnaa,vrnag,vtrans,vua,vx,vxd)
  names(V)<-c("vada","vade","vadna","vadrnr","vampd","vaprt","varna","vasuc","vasli",
              "vdada","vden","vdgnuc","vdnaa","vdnag","vgdna","vgdrnr","vgmpr","vgmps",
              "vgnuc","vgprt","vgrna","vgua","vhprt","vhx","vhxd","vimpd",
      "vinuc","vmat","vpolyam","vprpps","vpyr","vrnaa","vrnag","vtrans","vua","vx","vxd")
  list(XP,V)
}

(y0=c(PRPP=5,IMP=100,SAMP=.2,Ado=2500,SAM=4,Ade=1,XMP=25,GMP=400,dAdo=6,
      dGMP=3,RNA=28600,DNA=5160,HX=10,Xa=5,Gua=5,UA=100))
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=ode(y=y0,times=seq(-20,70,1),fpur,events = list(data = dPRPP10) ) )
graphics.off()
windows(height=8,width=10,xpos=-100)
plot(out,which=c("PRPP","vden","IMP","HX","Gua","vaprt","XMP","Xa","UA"))
# The PRPP bolus pushed both de novo and salvage purine synthesis. Additional de
# novo in the end means additional uric acid coming out, i.e. UA AUC above 
# baseline is greater than the initial drop below baseline. Note also that lags
# increase with distance from the PRPP perturbation point.

# plot(out)

# show that it really did kick PRPP by a factor of 10
t(ode(y=y0,times=seq(-2,2,.5),fpur, events = list(data = dPRPP10) )[,1:2])
t(ode(y=y0,times=seq(-.5,.5,.1),fpur, events = list(data = dPRPP10) )[,1:2])
t(ode(y=y0,times=seq(-.05,.05,.01),fpur, events = list(data = dPRPP10) )[,1:2])
