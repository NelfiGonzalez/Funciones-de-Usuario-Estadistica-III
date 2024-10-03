#Nueva versión de función Descomp.Loess, incorporando I.P aproximados
Descomp.Loessv02=function(serie.ajuste,tipo.descomp,grado,criterio,h,level=0.95){
library(fANCOVA)
library(forecast)
s=tsp(serie.ajuste)[3]
tajuste=1:length(serie.ajuste)
tpron=(length(serie.ajuste)+1):(length(serie.ajuste)+h)
if(end(serie.ajuste)[2]<s){
start.pron=c(end(serie.ajuste)[1],end(serie.ajuste)[2]+1)
}
if(end(serie.ajuste)[2]==s){
start.pron=c(end(serie.ajuste)[1]+1,1)
}

descom=decompose(serie.ajuste,type=tipo.descomp)
St=descom$seasonal
estacionini=cycle(serie.ajuste)[1]
if(estacionini==1){
deltas_i=descom$figure
}
if(estacionini!=1){
j=estacionini;deltas_i=c(descom$figure[(s-j+2):s],descom$figure[1:(s-j+1)])
}
efectos=deltas_i
if(s==12){
names(efectos)=c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
}
if(s==4){
names(efectos)=c("s1","s2","s3","s4")
}
efectos=data.frame(efectos)
names(efectos)=""
cat("Efectos estacionales estimados")
print(efectos)

indi=seasonaldummy(serie.ajuste,h)
indi2=cbind(indi,1-apply(indi,1,sum))
Stnuevo=ts(indi2%*%deltas_i,frequency=s,start=start.pron)

if(tipo.descomp=="additive"){
ytd=serie.ajuste-St
}
if(tipo.descomp=="multiplicative"){
ytd=serie.ajuste/St
}
ajusteLoess1=loess.as(tajuste,ytd,degree=grado,criterion=criterio,family="gaussian",plot=F)
cat("\n")
cat("Resumen Loess sobre serie desestacionalizada:")
cat("\n")
print(summary(ajusteLoess1))
alfa.optim1=ajusteLoess1$pars$span 
nep1=round(ajusteLoess1$enp)
p1=nep1+s-1 #numero aproximado de parametros en ajuste Modelo
Tt1=ts(fitted(ajusteLoess1),frequency=s,start=start(serie.ajuste))
LOESS=predict(loess(ytd~tajuste,span=alfa.optim1,degree=grado,control=loess.control(surface="direct")),data.frame(tajuste=tpron),se=TRUE)
Ttnuevo1=ts(LOESS$fit,freq=s,start=start.pron)
if(tipo.descomp=="additive"){
ythat1=St+Tt1
ytpron1=Ttnuevo1+Stnuevo
et1=serie.ajuste-ythat1
df1=length(serie.ajuste)-nep1+s-1 
MSE1=sum(et1^2)/df1 
varpred=MSE1+(LOESS$se.fit)^2
LIP=ytpron1-qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
LSP=ytpron1+qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
}
if(tipo.descomp=="multiplicative"){
ythat1=St*Tt1
ytpron1=Ttnuevo1*Stnuevo
et1=serie.ajuste-ythat1
df1=length(serie.ajuste)-nep1+s-1 
MSE1=sum(et1^2)/df1 
varpred=MSE1+(Stnuevo^2)*(LOESS$se.fit)^2
LIP=ytpron1-qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
LSP=ytpron1+qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
}
tablapron1=ts(cbind(forecast=ytpron1,lower=LIP,upper=LSP),freq=s,start=start.pron)
result=list(deltasi=deltas_i,alfa.optim=alfa.optim1,nep=nep1,p=p1,St=St,Tt=Tt1,ytd=ytd,fitted=ythat1,forecast=tablapron1,residuals=et1,MSE=MSE1)
result
}

