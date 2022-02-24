stratMH<-function(X,L=3,n=NULL,cvt=0.1,p=50,pe=0.3,pm=0.3,maxgen=50,AV=TRUE,npar=1)
{
strathh<-function(X,L=3,n=NULL,cvt=0.1,p=50,pe=0.3,pm=0.3,maxgen=50,AV=TRUE)
{
  
  Build<-function(L,K,Ns)
  {
    lower<-2
    upper<-K-2*(L-1)
    Y<-matrix(0,ncol=L,nrow=Ns)
    for(i in 1:Ns)
    {
      x<-vector(mode="integer",L)
      x[1]<-sample(c(lower:upper),1)
      for(h in 2:(L-1)) 
      {upper<-K-sum(x[1:(h-1)])-2*(L-h)
      x[h]<-sample(c(lower:upper),1)
      if (upper==lower) {x[h]<-lower}
      }  
      x[L]=K-sum(x)
      Y[i,]<-x
      upper<-K-2*(L-1)
    }
    return(Y)
  }  
  
  ########## Calcula a função objetivo ###########################
  Fobj<-function(x,L,n,cvt,tx,X,Xu)
  {
    NhVh<-function(X,Xu,bi,bj)
    { 
      
      z=which(X>=Xu[bi] & X<=Xu[bj])
      Nh<-length(z)
      Vh<-ifelse(AV==FALSE,var(X[z]),var(X[z])*(Nh-1)/Nh)
      return(list(Nh=Nh,Vh=Vh))
    }  
    Nh<-rep(0,L)
    Vh<-rep(0,L)
    bi<-1
    bj<-x[1]
    sa<-cumsum(x)
    for(h in 1:(L-1)) 
    { ca<-NhVh(X,Xu,bi,bj)
    Nh[h]<-ca$Nh
    Vh[h]<-ca$Vh
    bi<-bj+1
    bj<-sa[h+1]
    }
    ca<-NhVh(X,Xu,bi,bj)  
    Nh[L]<-ca$Nh
    Vh[L]<-ca$Vh
    if(is.null(cvt))
    {A=BSSM_FD(Nh,Vh,tx,n)
    nh=A$nh
    cv=A$cvs
    return(c(Nh,Vh,nh,n,cv))
    }
    if(is.null(n))
    { A=BSSM_FC(Nh,Vh,tx,cvt) 
    nh=A$nh
    n=A$n
    cv=A$cvs
    return(c(Nh,Vh,nh,cv,n))
    }
    
  }  
  
  Fobj2<-function(x,L,n,cvt,tx,X,Xu)
  {
    NhVh<-function(X,Xu,bi,bj)
    { 
      z=which(X>=Xu[bi] & X<=Xu[bj])
      Nh<-length(z)
      Vh<-ifelse(AV==FALSE,var(X[z]),var(X[z])*(Nh-1)/Nh)
      return(list(Nh=Nh,Vh=Vh))
    }  
    Nh<-rep(0,L)
    Vh<-rep(0,L)
    bi<-1
    bj<-x[1]
    sa<-cumsum(x)
    for(h in 1:(L-1)) 
    { ca<-NhVh(X,Xu,bi,bj)
    Nh[h]<-ca$Nh  
    Vh[h]<-ca$Vh
    bi<-bj+1
    bj<-sa[h+1]
    }
    ca<-NhVh(X,Xu,bi,bj)  
    Nh[L]<-ca$Nh
    Vh[L]<-ca$Vh
    if(is.null(cvt))
    {A=BSSM_FD(Nh,Vh,tx,n)
    nh=A$nh
    cv=A$cvs
    return(c(Nh,Vh,nh,n,cv))
    }
    if(is.null(n))
    { n=sum(Nh*sqrt(Vh))^2/(cvt^2*tx^2+sum(Nh*Vh))
    nh=n*Nh*sqrt(Vh)/sum(Nh*sqrt(Vh))
    z=which(nh>Nh)
    if (length(z>0)) 
    {A=BSSM_FC(Nh,Vh,tx,cvt)
    nh=A$nh
    n=A$n
    } 
    cv=sqrt(sum(Nh^2*Vh*(1-nh/Nh)))/tx
    return(c(Nh,Vh,nh,cv,n))
    }
    
  }  
  
  
  
  
  
  
  ###Crossover#################################################
  cross<-function(se,sne,L,te,tne,tc,K)
  {
    S<-matrix(0,ncol=L,nrow=tc)
    tcm<-tc%/%2
    j<-0
    while(j<tcm)
    {p1<-sample(te,1)
    p2<-sample(tne,1)
    c1<-se[p1,]
    c2<-sne[p2,]
    ca=c1
    cb=c2
    dc<-which(abs(c1-c2)>0)
    if (length(dc)>0)
    {j<-j+1
    l<-ifelse(length(dc)>1,sample(dc,1),dc)
    es<-setdiff(1:L,l)
    if (c1[l]<c2[l]) {aux<-c1;c1=c2;c2=aux}
    d<-c1[l]-c2[l]
    aux<-c1[l]
    c1[l]<-c2[l]
    c2[l]<-aux

    soma1<-sum(c1[-l])
    prop<-round(d*c1[-l]/soma1)
    for(h in 1:(L-1)) {c1[es[h]]<-c1[es[h]]+prop[h]}
    soma2<-sum(c2[-l]-2)
    prop<-round(d*(c2[-l]-2)/soma2)
    for(h in 1:(L-1)) {c2[es[h]]<-c2[es[h]]-prop[h]}
    soma1<-sum(c1)
    qm<-which.max(c1)
    d<-soma1-K
    c1[qm]<-c1[qm]-d
    
    soma2<-sum(c2)
    qm<-which.max(c2)
    d<-soma2-K
    c2[qm]<-c2[qm]-d
    S[2*j-1,]<-c1
    S[2*j,]<-c2
    }
    }
    return(S)
  }  
  
  
  cross2<-function(se,sne,L,te,tne,tc,K)
  {
    S<-matrix(0,ncol=L,nrow=tc)
    ST<-matrix(0,ncol=L,nrow=2*L)
    FF<-matrix(0,ncol=3*L+2,nrow=tc)
    tcm<-tc%/%2
    j<-0
    while(j<tcm)
    { p1<-sample(te,1)
    p2<-sample(tne,1)
    c1<-se[p1,]
    c2<-sne[p2,]
    ca=c1
    cb=c2
    
    dc<-which(abs(c1-c2)>0)
    if (length(dc)>0)
    {j<-j+1
    for(f in 1:length(dc))
    { l<-dc[f]
    es<-setdiff(1:L,l)
    if (c1[l]<c2[l]) {aux<-c1;c1=c2;c2=aux}
    d<-c1[l]-c2[l]
    aux<-c1[l]
    c1[l]<-c2[l]
    c2[l]<-aux

    soma1<-sum(c1[-l])
    prop<-round(d*c1[-l]/soma1)
    for(h in 1:(L-1)) {c1[es[h]]<-c1[es[h]]+prop[h]}
    soma2<-sum(c2[-l]-2)
    prop<-round(d*(c2[-l]-2)/soma2)
    for(h in 1:(L-1)) {c2[es[h]]<-c2[es[h]]-prop[h]}  
    soma1<-sum(c1)
    qm<-which.max(c1)
    d<-soma1-K
    c1[qm]<-c1[qm]-d
    soma2<-sum(c2)
    qm<-which.max(c2)
    d<-soma2-K
    c2[qm]<-c2[qm]-d
    
    m1=which(c1<2)
    m2=which(c2<2)
    if (length(m1)>0) 
    {dz=2-c1[m1];c1[m1]<-2;qmax=which.max(c1);c1[qmax]=c1[qmax]-sum(dz)}
    if (length(m2)>0) 
    {dz=2-c2[m2];c2[m2]<-2;qmax=which.max(c2);c2[qmax]=c2[qmax]-sum(dz)}
    
    
    ST[2*f-1,]<-c1
    ST[2*f,]<-c2
    c1=ca
    c2=cb
    } 
    
    ST=ST[which(duplicated.matrix(ST)==FALSE),]
    ST=ST[which(apply(ST,1,sum)==K),]
    
    
    SM=t(apply(ST,1,function(x) Fobj2(x,L,n,cvt,tx,X,Xu)))
    qx=order(SM[,3*L+2])
 
    S[2*j-1,]<-ST[qx[1],]
    S[2*j,]<-ST[qx[length(qx)],]
    FF[2*j-1,]<-SM[qx[1],]
    FF[2*j,]<-SM[qx[length(qx)],]  
    ST<-matrix(0,ncol=L,nrow=2*L)
    }
    }
    
    return(list(S=S,FF=FF))
  }  
  
  
  
  mcut<-function(sbest,Xu,L)
  {
    bk<-rep(0,L-1)
    for(h in 1:(L-1))
    {
      bk[h]<-Xu[sbest[h]]
      Xu<-Xu[-c(1:sbest[h])]
    }
    return(bk)
  }
  
  
  
  
  #################Main Program ###############################
  

  library(MultAlloc)
  tx=sum(X)
  Xu<-sort(unique(X))
  K<-length(Xu)
 
  tempo<-proc.time()
  s=Build(L,K,2*p)
  f=t(apply(s,1,function(x) Fobj2(x,L,n,cvt,tx,X,Xu)))
  
  
  idx<-order(f[,3*L+2])
  s<-s[idx[1:p],]
  f<-f[idx[1:p],]
  ne<-round(pe*p)
  nm<-round(pm*p)
  fbest=Inf
  ng<-0
  ibest<-0
  fmg<-Inf
  Message<-ifelse(is.null(n)==FALSE,"coefficient of variation = ","Sample Size = ")
  while((ng<maxgen) & (ng-ibest<=round(0.3*maxgen)))
  {
    ng<-ng+1
    idx<-order(f[,3*L+2],decreasing = FALSE)
    s<-s[idx,]
    f<-f[idx,]
    
    if (f[1,3*L+2]<fbest)
    {  sbest<-s[1,]
    fbest<-f[1,3*L+2]
    fg<-f
    cat("Generation ",ng,Message,fbest,"\n")
    ibest<-ng
    } 
    
  
    smutacao<-Build(L,K,nm)
    fmutacao<-t(apply(smutacao,1,function(x) Fobj2(x,L,n,cvt,tx,X,Xu)))
  
    
    scross<-cross2(s[1:ne,],s[(ne+1):p,],L,ne,p-ne,p-ne-nm,K)
    fcross<-scross$FF
    scross<-scross$S
    
    s<-rbind(s[1:ne,],smutacao,scross)  
    f<-rbind(f[1:ne,],fmutacao,fcross)
    
  }
  
  fg<-Fobj(sbest,L,n,cvt,tx,X,Xu)
  if(is.null(cvt)==TRUE) {r=c("n","cv")} else {r=c("cv","n")}
  cat("Best Solution ",r[1]," = ",fg[3*L+1]," ",r[2]," = ",fg[3*L+2],"\n")
  tempo<-(proc.time()-tempo)[3]
  
  Nh=fg[1:L]
  Sh2=fg[(L+1):(2*L)]
  nh=fg[(2*L+1):(3*L)]
  fobj=fg[3*L+2]
  bk<-mcut(sbest,Xu,L)
  cv=sqrt(sum(Nh^2*Sh2/nh*(1-nh/Nh)))/tx
  return(list(xbest=sbest,bk=bk,Nh=Nh,nh=nh,
              Sh2=Sh2,fobj=fobj,cv=cv,ta=tempo,ibest=ibest))  
  
}
if (npar>1)
 {library(parallel)
  library(MultAlloc)
  nucleos<-detectCores(logical=FALSE)
  clust<-makeCluster(nucleos)
  clusterEvalQ(clust, library(MultAlloc)) 
  clusterEvalQ(clust,"strathh")
  clusterExport(clust,"X")
  Lpar<-rep(L,npar)
  Time_CPU<-proc.time()
  s<-clusterApplyLB(cl=clust,Lpar,function(lx) strathh(X,L=lx,n,cvt,p,pe,pm,maxgen,AV))
  stopCluster(clust) 
  Time_CPU<-(proc.time()-Time_CPU)[3]
  ixmin<-which.min(sapply(s,function(x) x$fobj))
  s<-s[[ixmin]]
  s$ta<-Time_CPU
  
 } else {s<-strathh(X,L,n,cvt,p,pe,pm,maxgen,AV)}  
 
return(s)
  
}