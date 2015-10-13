library(SBMLR)
setwd(file.path(system.file(package="SBMLR"), "BMCcancerFolates")) #default dump site 
morr=readSBMLR(file.path(system.file(package="SBMLR"), "models/morrison.r"))  
out1=sim(morr,seq(-20,0,1))
morr$species$EMTX$ic=1
out2=sim(morr,0:30)
outs=data.frame(rbind(out1,out2))
attach(outs)
par(mfrow=c(3,3))
plot(time,EMTX,type="l",xlab="Hours",ylab="EMTX")
plot(time,FH2b+FH2f,type="l",xlab="Hours",ylab="DHF")
plot(time,CH2FH4,type="l",xlab="Hours",ylab="CH2THF")
plot(time,CHOFH4,type="l",xlab="Hours",ylab="CHOTHF")
plot(time,FH4,type="l",xlab="Hours",ylab="THF")
plot(time,CH3FH4,type="l",xlab="Hours",ylab="CH3THF")
plot(time,MTHFR,type="l",xlab="Hours",ylab="MTHFR")
plot(time,GARFT,type="l",xlab="Hours",ylab="GARFT")
plot(time,TYMS,type="l",xlab="Hours",ylab="TYMS")
par(mfrow=c(1,1))
detach(outs)
morr$species$EMTX$ic=0
dev.copy(pdf,file="fig2.pdf", width = 7, height = 7)
dev.off()


