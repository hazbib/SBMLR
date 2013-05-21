library(SBMLR) 
curto=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=simulate(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))
# out is a dataframe with class set to deSolve to overload plot, see ?plot.deSolve

library(SBMLR)
# sod=readSBMLR(system.file("models", "sod.R", package = "SBMLR"))  
sod=readSBMLR("/Users/radivot/Downloads/SBMLR/inst/models/sod.R")  
summary(sod)
out=simulate(sod,times=seq(0,1000,1) ) 
head(out)
plot(out)

setwd("/Users/radivot/Downloads/SB")
saveSBML(sod,"sod.xml")
saveSBMLR(sod,"sod.R")
# like it, so write to package
# saveSBML(sod,"/Users/radivot/Downloads/SBMLR/inst/models/sod.xml")

sodX=readSBML("sod.xml")
sodR=readSBMLR("sod.r")
head((sodX==sodR)$species)
head((sodX==sodR)$reactions)
# Values in these two dataframes are TRUE where the initial concentrations,
# fluxes, and reaction rate laws (as strings) are equal.

summary(sod)
sod

(SOadd <- data.frame(var = "SO", time = 0, value = 1e9,method = "mult"))
out=simulate(sod,times=seq(-1000,1000,1),events = list(data = SOadd) ) 
plot(out)

tail(out)
(ssic=out[dim(out)[1],2:8])
out=simulate(sod,X0=ssic,times=seq(-10,10,.1),events = list(data = SOadd) ) 
plot(out)

out=simulate(sod,X0=ssic,times=seq(-1,1,.01),events = list(data = SOadd) ) 
plot(out)

