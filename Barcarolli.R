Barcarolli<-function(X,L=3,cvt=0.1,min_n=2)
{
library(SamplingStrata)
X=sort(X)
cities <- as.data.frame(list(id=c(1:length(X)),
                                 X1=as.numeric(X),
                                 Y1=as.numeric(X),
                                 domainvalue=rep(1,length(X))))
cv <- as.data.frame(list(DOM="DOM1",CV1=cvt,domainvalue=1 ))
cvf=Inf
Tg<-0
rstart<-0
DF<-Inf

while((cvf>cvt) & (rstart<10))
 {rstart<-rstart+1
  tempo<-proc.time()
  solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = cities,
                            minnumstr = min_n,
                            iter = 50,
                            pops = 20,
                            nStrata = L,
                            parallel = FALSE,
                            showPlot=FALSE)
  tempo<-(proc.time()-tempo)[3]    
  dev.off()
  n=sum(round(solution$aggr_strata$SOLUZ))
  Nh=solution$aggr_strata$N
  Vh=solution$aggr_strata$S1^2*(Nh-1)/Nh

  nh=round(solution$aggr_strata$SOLUZ)
  
  cvf=sqrt(sum(Nh^2*Vh/nh*(1-nh/Nh)))/sum(X)
  if (cvf>cvt)
   {if (abs(cvf-cvt)<DF) 
      {DF<-abs(cvf-cvt)
       nbest=n
       nhbest=nh
       Nhbest=Nh
       Vhbest=Vh
       cvfbest<-cvf
      }
   } else
      {nbest=n
       nhbest=nh
       Nhbest=Nh
       Vhbest=Vh
       cvfbest<-cvf
      }
  Tg=Tg+tempo
  LF<-length(nhbest)
}
if (LF<L) {cvfbest=Inf}
return(list(n=nbest,nh=nhbest,Nh=Nhbest,Vh=Vhbest,tempo=Tg,cv=cvfbest,rstart=rstart))
}  