# This is the SOD model of Kowald et al, JTB 238 (2006) 828–840. 
ic=c(
		SO		=   0, #    4e-11,
# 		SODI	=    5e-6,
		SODII	=    5e-6,
		H2O2	=  0,#     1e-9,
		OH		= 0,#      1e-23,
		L     =  0,#      1e-11,
		LOO   =       0,
		LOOH	=       0
# 		LOO   =       7e-8,   # roughly their final values in Fig. 2 JTB 2006
# 		LOOH	=       1.5e-6
)

parameters=c( 
		k1  =6.6e-7,
		k2  =1.6e9,
		k3  =1.6e9,
		k4  =1e5,
		k5  =2e4,
		k6  =1,
		k7  =3.4e7,
		k9  =1e6,
		k10 =1e3,
		k11 =2.5e8,
		k12 =0.38,
		k13a=8.7e-3,
		k13b=8.7e-3,
		k17 =3e4,
		k18 =7,
		k19 =8.8e4,
		SOD 	   =   1e-5,
		catalase =   1e-5
)

sod<-function(t, state, parameters) {
  state[state<0]=1e20
	with(as.list(c(state, parameters)),{
 				SODI=SOD-SODII
				etcV1=k1
				sooV2=k2*SO*SODII
				sorV3=k3*SO*SODI
				loorV4=k4*SO*LOO
				habWeiV5=k5*SO*H2O2
				sodFenV6=k6*H2O2*SODII
				catGpxV7=k7*H2O2*catalase
				leakV9=k9*OH
				hsoLipPeroxV10=k10*SO/100 # this is also the flux through the fast reaction
				hoLipPeroxV11=k11*OH
				loohDeOxV12=k12*LOOH
				sodoV13a=k13a*SODI
				sodrV13b=k13b*SODII
				lipRadOxV17=k17*L
				chainRxV18=k18*LOO
				looSelfV19=k19*LOO^2
				dSO	 = etcV1 - sooV2 - sorV3 -loorV4 - habWeiV5 -hsoLipPeroxV10
				dH2O2  =	sorV3 + hsoLipPeroxV10 - habWeiV5 -sodFenV6	-catGpxV7			
				dOH	 =	habWeiV5 + 2*sodFenV6	- leakV9 - hoLipPeroxV11 # I think this 2 should be 1	
				dL		 = hsoLipPeroxV10	+ hoLipPeroxV11 + chainRxV18 - lipRadOxV17			
 				dLOO   =	lipRadOxV17	- chainRxV18 -loorV4 - 2*looSelfV19		
				dLOOH	 = chainRxV18 - loohDeOxV12 + loorV4
# 				dSODI = sooV2 + sodrV13b - sorV3 - sodoV13a # replacing constraint with ODE didn't help
				dSODII = sorV3 + sodoV13a - sooV2 - sodrV13b
				list(c(dSO, dH2O2, dOH, dL, dLOO, dLOOH, dSODII),
						c(etcV1=etcV1, sooV2=sooV2, sorV3=sorV3, loorV4=loorV4, 
								habWeiV5=habWeiV5, sodFenV6=sodFenV6, catGpxV7=catGpxV7, 
								leakV9=leakV9, hsoLipPeroxV10=hsoLipPeroxV10, hoLipPeroxV11=hoLipPeroxV11, 
								loohDeOxV12=loohDeOxV12, sodoV13a=sodoV13a, sodrV13b=sodrV13b, 
								lipRadOxV17=lipRadOxV17, chainRxV18=chainRxV18, looSelfV19=looSelfV19)
				)
			}) # end with.
}

library(deSolve)
out <- ode(y = ic, times = seq(0,100,10), func = sod, parms = parameters)

out <- ode(y = ic, times = seq(0,100,10), func = sod, parms = parameters,
           atol = 1e-1, rtol = 1e-1,method="adams")
head(out)

ode(y = ic, times = seq(0,100,10), func = sod, parms = parameters,method="adams")
