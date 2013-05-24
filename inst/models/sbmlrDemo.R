library(SBMLR) 
curto=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=sim(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))
# out is a dataframe with class set to deSolve to overload plot, see ?plot.deSolve

library(SBMLR)
# sod=readSBMLR(system.file("models", "sod.R", package = "SBMLR"))  
sod=readSBMLR("/Users/radivot/Downloads/SBMLR/inst/models/sod.R")  
summary(sod)
out=sim(sod,times=seq(0,1000,1) ) 
tail(out)
(ssic=out[dim(out)[1],2:8])

(SOadd <- data.frame(var = "SO", time = 0, value = 1e9,method = "mult"))
out=sim(sod,X0=ssic,times=seq(-10,10,.1),events = list(data = SOadd) ) 
plot(out)

out=sim(sod,X0=ssic,times=seq(-1,1,.01),events = list(data = SOadd) ) 
plot(out) # now you can actually the spike in of SuperOxide

