library(SBMLR)
# Axel Kowald's model in "A systems biological analysis links ROS metabolism to 
# mitochondrial protein quality control" Mechanisms of Ageing and Development
# 133 (2012) 331–337
# This was redrawn in JDesigner from Axel's CellDesigner original
sod=readSBML(file.path(system.file(package="SBMLR"), "models/sod2012JD.xml"))  
# Note that JD puts all params in global list
summary(sod)$V0
summary(sod)$globalVec
SOD=10^seq(-7,-3,0.2)
out=NULL
for (i in SOD) {
  sod$globalParameters$O2detox_SOD=i
  (out=rbind(out,sim(sod,times=seq(0,5) )[6,c("O2rad","H2O2")] ))
}
(df=as.data.frame(cbind(SOD,out)))
par(mfrow=c(2,1))
with(df,plot(SOD,O2rad,log="x",type="l"))
with(df,plot(SOD,H2O2,log="x",type="l"))
par(mfrow=c(1,1))

####################
library(SBMLR)
setwd("/Users/radivot/Downloads/SB/sod")
sod=readSBML("sod2012JD.xml") 

# these were both zero for Fig. 2, set here now for Fig. 3
sod$globalParameters$PRXoxidation_k9=80  
sod$globalParameters$CLPPoxidation_k13=300

SOD=10^seq(-6,-3,0.2)
Tss=1e5
out=NULL
for (i in SOD) {
  sod$globalParameters$O2detox_SOD=i
  (out=rbind(out,sim(sod,times=seq(0,Tss,by=Tss/10) )[11,c("PRX","CLPP")] ))
}
(df=as.data.frame(cbind(SOD,out)))  # shape is right, but not enough movement?????
par(mfrow=c(2,1))
with(df,plot(SOD,PRX,log="x",type="l"))
with(df,plot(SOD,CLPP,log="x",type="l"))
par(mfrow=c(1,1)) 


####################  Reproduce top curves in Figures 4a and b
library(SBMLR)
setwd("/Users/radivot/Downloads/SB/sod")
sod=readSBML("sod2012JD.xml") 
summary(sod)
paraQ=seq(0,1.3e-6,0.1e-6)
Tss=1e5
out=NULL
for (i in paraQ) {
  sod$species$paraquat$ic=i
  (out=rbind(out,sim(sod,times=seq(0,Tss,by=Tss/10) )[11,c("O2rad","H2O2")]))
}
(df=as.data.frame(cbind(paraQ,out)))  # shape is right, but not enough movement?????
with(df,plot(paraQ,H2O2,type="l",ylim=c(0,7e-7)))

sod$globalParameters$PRXoxidation_k9=80  
sod$globalParameters$CLPPoxidation_k13=300
out=NULL
for (i in paraQ) {
  sod$species$paraquat$ic=i
  (out=rbind(out,sim(sod,times=seq(0,Tss,by=Tss/10) )[11,"H2O2",drop=F]))
}
(df=as.data.frame(cbind(paraQ,out)))  
with(df,lines(paraQ,H2O2)) # note, plot on the bottom if from Fig. 4a, not lower in 4b

out=NULL
sod$globalParameters$O2detox_SOD=1e-6
for (i in paraQ) {
  sod$species$paraquat$ic=i
  (out=rbind(out,sim(sod,times=seq(0,Tss,by=Tss/10) )[11,"H2O2",drop=F]))
}
(df=as.data.frame(cbind(paraQ,out)))  
with(df,lines(paraQ,H2O2,col="red")) # here red is lower curve in 4b

# install.packages("gglot2") # not available for R 3.0.1??????
# library(gglot2)
