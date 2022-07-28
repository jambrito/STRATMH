RunAlgs<-function(X,cva,L)
{
RKozak<-function(X,LL,cvt=0.1,min_n=2)
{   cpu_time<-proc.time()
    Kozak=strata.LH(X, CV = cvt,Ls = LL,alloc = c(0.5, 0, 0.5), takeall = 0, algo = "Kozak", model="none",
                    algo.control=list(maxiter=500000))
    nh<-Kozak$nh
    Nh<-Kozak$Nh
    Vh<-Kozak$varh
    n<-sum(nh)
    cv<-sqrt(sum(Nh^2*Vh/nh*(1-nh/Nh)))/sum(X)
    cpu_time<-(proc.time()-cpu_time)[3]
    return(list(n=n,bk=Kozak$bh,nh=nh,Nh=Nh,Vh=Vh,cpu_time=cpu_time,cv=cv))
}  
RLH<-function(X,LL,cvt=0.1,min_n=2)
{
cpu_time<-proc.time()  
LH=strata.LH(X, CV = cvt, Ls = LL,alloc = c(0.5, 0, 0.5), takeall = 0, algo = "Sethi", model="none")
nh<-LH$nh
n<-sum(nh)
Nh<-LH$Nh
Vh<-LH$varh
cv<-sqrt(sum(Nh^2*Vh/nh*(1-nh/Nh)))/sum(X)
if (is.nan(cv)==TRUE) {n=10^5;cv=10} else {if (cv>cvt) {n<-100*n}}
cpu_time<-(proc.time()-cpu_time)[3]
return(list(n=n,bk=LH$bh,nh=nh,Nh=Nh,Vh=Vh,cpu_time=cpu_time,cv=cv))
 }  
 
 


library(stratification)
library(stratvns)
source("stratMH.R")
source("Barcarolli.R")
  
Kozak=RKozak(X, cvt=cva,LL = L)
Vns=STRATVNS(X,L=L,cv=cva,maxstart=3,parallelize=TRUE)
LH=RLH(X, cvt = cva, LL = L)
Ba=Barcarolli(X,L=L,cvt=cva,min_n=2)
stratmh=stratMH(X,L,cvt=cva)


algs=c("Kozak","VNS","LH","Barca","STRATMH")
ns=c(Kozak$n,Vns$n,LH$n,Ba$n,stratmh$fobj)
cvs=c(Kozak$cv,Vns$cv,LH$cv,Ba$cv,stratmh$cv)
G=data.frame(algs=algs,ns=ns,cvs=round(cvs,6))
return(G)
}

#Create a directory called experiments
#copy populations to experiments directory
#Copy this file and files stratMH.R and Barcarolli.R to  directory experiments
setwd("\\Experiments") 
X1=scan("BeefFarms.txt")
X2=scan("chi1.txt",dec=".")
s1=RunAlgs(X1,0.1,3)
s2=RunAlgs(X2,0.05,4)








