#Author(s):  A.I. McLeod and Y. Zhang
`AR1Est` <-
function(z, MeanValue=0){
stopifnot(length(z)>1)
n=length(z)
m <- MeanValue
mz <- rep(m,n)
a <- sum((z-mz)^2)
b <- sum((z[-1]-mz[-1])*(z[-n]-mz[-n]))
c <- sum((z[c(-1,-n)]-mz[c(-1,-n)])^2)
i <- complex(1,0,1)
x <- ((-16)*b^3+18*a*b*c+24*b^3*n-27*a*b*c*n-9*b*c^2*n-12*b^3*n^2+9*a*b*c*n^2+27*b*c^2*n^2+2*b^3*n^3-18*b*c^2*n^3)
y <- (-b^2*(-2+n)^2+3*c*(-1+n)*(-a-c*n))
f <- complex(1,x^2+4*y^3,0)
z <- (x+sqrt(f))^(1/3)
g <- x^2+4*y^3
z1 <- (x+(-g)^(1/2)*i)^(1/3)
part1 <- (n-2)*b/(3*c*(n-1))
part2 <- (1-sqrt(3)*i)*y/(3*2^(2/3)*c*(n-1)*z)
part3 <- (1+sqrt(3)*i)*z/(6*2^(1/3)*c*(n-1))
Re(part1+part2-part3) 
}


`ARToPacf` <-
function(phi){
phik=phi
L=length(phi)
if(L==0) return(0)
pi=numeric(L)
for (k in 1:L){
    LL=L+1-k
    a <- phik[LL]
    pi[L+1-k] <- a
    phikp1 <- phik[-LL]
    if(is.na(a) || abs(a)==1)
        stop("transformation is not defined, partial correlation = 1")
    phik <- (phikp1+a*rev(phikp1))/(1-a^2)
    }
pi
}

`ChampernowneD` <-
function(z, p, MeanZero=FALSE){
n<-length(z)
if(MeanZero) y<-z
else y<-z-mean(z)
x0<-x<-y
for (i in 1:p)
    x<-c(x,c(rep(0,i),y[1:(n-i)]))
x<-matrix(x, nrow=p+1, ncol=length(x0), byrow=TRUE)
C<-c(x%*%x0)
A<-toeplitz(C)
E<-matrix(0, nrow=p+1, ncol=p+1)
for (j in 1:p)
    for (i in 1:j){
         E[i+1,j+1] <- E[i,j]+y[i]*y[j]+y[n+1-i]*y[n+1-j]
         }
for (j in 1:(p+1))
    for (i in 1:j)
        E[j,i]=E[i,j]
A-E
}

`DetAR` <-
function(phi){
z<-ARToPacf(phi)
1/prod((1-z^2)^(1:length(phi)))
}


`FastLoglikelihoodAR` <-
function(phi, n, CD){
phis<-c(1,-phi)
LL <- -log(DetAR(phi))/2-(n/2)*log(sum(crossprod(phis,CD)*phis)/n)
if (!is.finite(LL)){
   LL<--1e35 }
LL
}


`GetFitAR` <-
function(z, p, ARModel="ARz", ...){
if (ARModel=="ARp")
    GetFitARpLS(z, p, ...)
else #also for ARModel="AR"
    GetFitARz(z, p, ...)
}

`GetFitARz` <-
function(z, pvec, MeanValue=0, ...){
stopifnot(length(z)>0, length(z)>=2*max(pvec), length(pvec)>0,all(pvec>=0))
pVec<-pvec[pvec>0]
y <- z-MeanValue
n<-length(y)
if (length(pVec)==0 || length(pvec)==0) 
    return(list(loglikelihood=-(n/2)*log(sum(y^2)/n),zetaHat=NULL,phiHat=NULL,convergence=0,algorithm="cubic"))
PMAX <- max(pVec)
if (PMAX == 1){
    phiHat <- AR1Est(y)
    LogL <- LoglikelihoodAR(phiHat,z)
    return(list(loglikelihood=LogL, zetaHat=phiHat, phiHat=phiHat, convergence=0))
    }
	
PEFF<-length(pVec)
CD<-ChampernowneD(y,PMAX,MeanZero=TRUE)
xpar<-numeric(PMAX)
EntropyAR<-function(x){
      if (max(abs(x))>0.999)
        out<-1e35
      else {
        xpar[pVec]<-x
        out<- -FastLoglikelihoodAR(PacfToAR(xpar),n,CD)
        }
out
}
xinit<-ARToPacf(ar.burg(y, aic=FALSE, order.max=PMAX, demean=FALSE)$ar)[pVec]
#Sometimes there are problems with "L-BFGS-B" -- it frequently tests the endpoints which is bad news due
#to numerical problems such as ARToPacf(PacfToAR(rep(0.99,20))) is not correct!
#So it is better to use "BFGS" with a penalty instead.
#ans<-optim(xinit,EntropyAR,method="L-BFGS-B", lower=rep(-0.9999,PEFF), upper=rep(0.9999,PEFF),control=list(trace=6),...)
ans<-optim(xinit,EntropyAR,method="BFGS", control="trace", ...)
alg<-1
if(ans$convergence != 0) {
    alg<-2
    warning("Convergence problem. convergence=", ans$convergence)
    warning("Trying Nelder-Mead algorithm ...")
    ans<-optim(xinit,EntropyAR,method="Nelder-Mead", ...)
    if(ans$convergence != 0)
        warning("Still convergence problem, convergence= ", ans$convergence)
}
zetaHat<-ans$par
zetas<-numeric(PMAX)
zetas[pVec]<-zetaHat
list(loglikelihood=-ans$value, zetaHat=ans$par, phiHat=PacfToAR(zetas),convergence=ans$convergence, algorithm=c("BFGS","Nelder-Mead")[alg],pvec=pvec)
}

`LoglikelihoodAR` <-
function(phi, z, MeanValue=0){
if(length(phi)==0) phi=0
phis<-c(1,-phi)
y<-z-MeanValue
n<-length(z)
-log(DetAR(phi))/2 - (n/2)*log(sum(crossprod(phis,ChampernowneD(y,length(phis)-1,MeanZero=TRUE))*phis)/n)
}



SelectModel <-
function(z, lag.max=15, ARModel=c("AR","ARz","ARp"), Criterion="default", Best=3, Candidates=5, t="default"){
#AR - order selection problem
#ARp - AR subset
#ARz - AR subset, partials
stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0, Candidates>0)
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
stopifnot(is.wholenumber(lag.max))
BestCandidates<-Candidates
IsValidCriterionQ <- Criterion %in% c("default", "AIC", "BIC", "UBIC", "EBIC", "BICq", "GIC")
if (!IsValidCriterionQ)
    stop("Criterion = ", Criterion, " not known.")
ARModel <- match.arg(ARModel)
if (Best > BestCandidates)
    BestCandidates<-Best
if (ARModel=="ARp") #subset ARp
    return(GetLeapsAR(z, lag.max=lag.max, Criterion=Criterion, Best=Best, Candidates=Candidates, t=t))
if (ARModel=="ARz")
    SubsetQ <- TRUE
else
    SubsetQ <- FALSE
method<-Criterion
if (Criterion == "default")
    if (SubsetQ)
        method <- "BICq"
    else
        method <- "BIC"
if (!SubsetQ && Criterion=="UBIC")
    method <- "BIC"
#set tuning parameter
P<-0.01
Q<-0.25
G<-1
if (method=="EBIC"  && t!="default")  G <- t
if (method=="QBIC"  && t!="default")  Q <- t
if (method=="GIC"   && t!="default")  P <- t
if (P>=0.25 || P<=0)
    stop("error: GIC tuning parameter invalid")
if (Q<=0 || Q>=1)
    stop("error: BICq tuning parameter invalid")
#approximate likelihood, "AR" or "ARz"
zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
n<-length(z)
LagRange<-1:lag.max
if (method=="UBIC"){
    mColNames<-list(c("p", "UBIC-Exact", "UBIC-Approx"))
    PENALTY1 <- log(n) + 2*lchoose(lag.max, 1)
    penalty<-log(n)*(1+LagRange)+2*lchoose(lag.max, LagRange)
    }
if (method=="EBIC"){
    mColNames<-list(c("p", "UBIC-Exact", "UBIC-Approx"))
    PENALTY1 <- log(n) + 2*G*lchoose(lag.max, 1)
    penalty<-log(n)*(1+LagRange)+2*G*lchoose(lag.max, LagRange)
    }
if (method=="BICq"){
    mColNames<-list(c("p", "BICq-Exact", "BICq-Approx"))
    PENALTY1 <- log(n) - 2*log(Q/(1-Q))
    penalty<-(1+LagRange)*PENALTY1
    }
if (method=="BIC"){
    PENALTY1 <- log(n)
    mColNames<-list(c("p", "BIC-Exact", "BIC-Approx"))
    penalty<-(1+LagRange)*PENALTY1
    }
if (method=="AIC"){
    PENALTY1 <- 2
    mColNames<-list(c("p", "AIC-Exact", "AIC-Approx"))
    penalty<-(1+LagRange)*PENALTY1
    }
if (method=="GIC"){
    mColNames<-list(c("p", "GIC-Exact", "GIC-Approx"))
    PENALTY1 <- qchisq(p=(1+sqrt(1-4*P))/2, df=1)
    penalty<-(1+LagRange)*PENALTY1
    }
if (SubsetQ)
    LagsEntering<-order(abs(zta),decreasing=TRUE)
else  
    LagsEntering<-1:lag.max
LLapprox <- n*log(cumprod(1-zta[LagsEntering]^2))
AnIC <- LLapprox + penalty
#
IndCandidates<-order(AnIC)[1:BestCandidates]
AnICexact<-numeric(BestCandidates+1)
if (SubsetQ){ #subset. AR model subset selection.
    m<-as.list(numeric(BestCandidates+1))
    for (isub in 1:BestCandidates){
        ModelLags<-sort(LagsEntering[1:IndCandidates[isub]])
        LL<-GetFitAR(z-mean(z), ModelLags)$loglikelihood
        k<-length(ModelLags)+1 #mean is included and k>=2 here
        if (method=="UBIC") {
            UBIC <- -2*LL + log(n)*k + 2*lchoose(lag.max+1, k) 
            AnICexact[isub]<-UBIC
            m[[isub]] <- list(p=ModelLags, UBIC=UBIC)
            }
        if (method=="EBIC") {
            EBIC <- -2*LL + log(n)*k + 2*G*lchoose(lag.max+1, k) 
            AnICexact[isub]<-EBIC
            m[[isub]] <- list(p=ModelLags, EBIC=EBIC)
            }
        if (method=="BICq") {
            BICq <- -2*LL + log(n)*k -2*(k*log(Q)+(lag.max+1-k)*log(1-Q))
            AnICexact[isub]<-BICq
            m[[isub]] <- list(p=ModelLags, BICq=BICq)
            }
        if (method=="AIC"){
            AIC <- -2*LL+2*k
            AnICexact[isub]<-AIC
            m[[isub]] <- list(p=ModelLags, AIC=AIC)
            }
        if (method=="BIC") {
            BIC <- -2*LL+log(n)*k
            AnICexact[isub]<-BIC
            m[[isub]] <- list(p=ModelLags, BIC=BIC)
            }
        if (method=="GIC") {
            GIC <- -2*LL+k*qchisq(p=(1+sqrt(1-4*P))/2, df=1)
            AnICexact[isub]<-GIC
            m[[isub]] <- list(p=ModelLags, GIC=GIC)
            }
    }
        #null model, note: k=1
        LL<-GetFitAR(z-mean(z), 0)$loglikelihood
        if (method=="UBIC") {#parameters=1, just mean
            UBIC <- -2*LL + PENALTY1
            AnICexact[BestCandidates+1]<-UBIC
            m[[BestCandidates+1]] <- list(p=0, UBIC=UBIC)
        }
       if (method=="EBIC") {
             EBIC <- -2*LL + PENALTY1
            AnICexact[BestCandidates+1]<-EBIC
            m[[BestCandidates+1]] <- list(p=0, EBIC=EBIC)
        }
        if (method=="BICq") {
            BICq <- -2*LL + PENALTY1 
            AnICexact[BestCandidates+1]<-BICq
            m[[BestCandidates+1]] <- list(p=0, BICq=BICq)
        }
        if (method=="AIC"){
            AIC <- -2*LL+PENALTY1
            AnICexact[BestCandidates+1] <- AIC
            m[[BestCandidates+1]]<-list(p=0,AIC=AIC)
            }
        if (method=="BIC") {
            BIC <- -2*LL+PENALTY1
            AnICexact[BestCandidates+1] <- BIC
            m[[BestCandidates+1]]<-list(p=0,BIC=BIC)
            }
        if (method=="GIC") {
            GIC <- -2*LL+PENALTY1
            AnICexact[isub]<-GIC
            m[[isub]] <- list(p=ModelLags, GIC=GIC)
            }
        #final model select based on exact likelihood
        i<-order(AnICexact)
        m<-m[i]
        m<-m[1:Best] 
        attr(m, "model")<-ARModel              
    }
else  { #non-subset. AR model order selection.
    AnICexact<-numeric(BestCandidates)
    AnICApprox<-numeric(BestCandidates)
    for (i in 1:BestCandidates){
        p<-LagsEntering[IndCandidates[i]]-1
        AnICApprox[i]<-AnIC[p+1]
        ans<-GetFitAR(z-mean(z), 0:p)
        LL<-ans$loglikelihood
#mean included in all models, k=p+1
        if (method=="AIC")
                penalty<-2*(p+1)
        if (method=="BIC")
                penalty<-log(n)*(p+1)
        if (method=="BICq")
                penalty<-log(n)*(1+p)-2*((p+1)*log(Q/(1-Q)))
        if (method=="EBIC") #is equivalent to BIC really!
                penalty<-log(n) + 2*G*lchoose(p+1, 1)
        if (method=="UBIC")
                penalty<-log(n) + 2*lchoose(p+1, 1)
        if (method=="GIC")
                penalty<-(1+p)*PENALTY1
        if (method=="BIC")
                penalty<-log(n)*(p+1)      
        AnICexact[i]<- -2*LL+penalty
    }
    m<-c(LagsEntering[IndCandidates]-1,AnICexact,AnICApprox)
    m<-matrix(m,ncol=3)
    if (Best==1) {
        m<-m[order(AnICexact),drop=FALSE]
        m<-m[1:Best,drop=FALSE]
    }
    else {
        m<-m[order(AnICexact),]
        m<-m[1:Best,]
    }
    if (Best > 1)
        dimnames(m)<-c(list(1:Best), mColNames)
}
if (SubsetQ) class(m)<-"Selectmodel"
if (Best > 1)
    m
else 
    if (is.list(m))
        m[[1]]$p
    else
        as.vector(m[1])
}

`PacfToAR` <-
function(zeta){
L=length(zeta)
if (L==0) return(numeric(0))
if (L==1) return(zeta)
phik=zeta[1]
for (k in 2:L){
    phikm1=phik
    phik=c(phikm1-zeta[k]*rev(phikm1),zeta[k])
    }
phik
}


"plot.Selectmodel" <-
function(x, ...){
ans<-x #for clarity we rename object
par(mar = c(5,4,4,5)+0.1)
on.exit(par(mar=c(5, 4, 4, 2) + 0.1)) #default
if (is.list(ans)){    
    X<-Y<-a<-numeric(0)
    for (i in 1:length(ans)){
        p<-ans[[i]]$p
        X<-c(X,p)
        aic<-ans[[i]][[2]]
        Y<-c(Y,rep(aic,length(p)))
        a<-c(a,aic)
        }
    ylab<-names(ans[[1]])[2]
    plot(X,Y,xlab="lags", ylab=ylab,pch=15, cex=2)
    abline(h=a, lty=2)
    title(main=paste(attr(ans,"model"),"model selection"))
    ytic<-a
    r<-round(exp((min(Y)-ytic)/2),3)
    axis(4,at=ytic, labels=r)
    mtext(side=4, line=2.7, text="Relative Plausibility")
    }
if (is.matrix(ans)){
    X<-ans[,1]
    Y<-ans[,2]
    ylab<-dimnames(ans)[[2]][2]
    plot(X,Y,ylab=ylab,pch=15,cex=2, xlab="p")
    title(main="AR(p) model selection")
    ytic<-Y
    r<-round(exp((min(Y)-ytic)/2),3)
    axis(4,at=Y,labels=r)
    mtext(side=4, line=2.7, text="Relative Plausibility")
    }
invisible()
}

