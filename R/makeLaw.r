"makeLaw"<-function(r,p,e){
# takes reactant list r, parameter list p and rate law R expression e 
# and makes a reaction rate law function out of them.
lawTempl=function(r,p){ }
if (is.null(p)) lawTempl=function(r){ }  
i=2
if(!is.null(p))
for (j in 1:length(p)){
body(lawTempl)[[i]]<-call("=",as.name(p[j]),call("[",as.name("p"),p[j]))
i=i+1}
#lawTempl
for (j in 1:length(r)){ 
body(lawTempl)[[i]]<-call("=",as.name(r[j]),call("[",as.name("r"),r[j]))
i=i+1}
body(lawTempl)[[i]]<-e
lawTempl
}
