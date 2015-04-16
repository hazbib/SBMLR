# Comment: relative to Curto an Morrison, wherein the *.R file was the editable 
# source code, I'm now switching to the use of JDesigner as my source editor
# (model IDE). Thus, for Kowald's SOD models I  now provide only the JDesigner
# xml file as the model source code. This is my new workflow. 
library(SBMLR)
sod=readSBML(file.path(system.file(package="SBMLR"), "models/SOD2006JD.xml"))  
summary(sod)
out=sim(sod,times=seq(0,1000,1) ) 
plot(out)  #reproduces Fig 2 of JTB 2006    
# Note: In RStudio, switch to console and hit return repeatedly to plot

(ssic=out[dim(out)[1],2:9]) # set ic to steady state
(SOadd <- data.frame(var = "SO", time = 0, value = 1e9,method = "mult"))
out=sim(sod,X0=ssic,times=seq(-10,10,.1),events = list(data = SOadd) ) 
plot(out) # note the SO looks like it dropped, but then why the LOOH go up
# Answer, sampling times are too big at 1 sec, need to get to sub millisecs
out=sim(sod,X0=ssic,times=seq(-.001,.001,.00001),events = list(data = SOadd) ) 
plot(out) # now you can actually see the upward spike of SuperOxide, SO

