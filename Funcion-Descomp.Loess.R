#Funci�n de usuario Descomp.Loess() para ajuste y pron�sticos  de validaci�n cruzada, por la combinaci�n del filtro de descomposici�n con loess
#Requiere la librer�a fANCOVA
#La funci�n combina el filtro de la descomposici�n cl�sica disponible en la funci�n R decompose(), para la estimaci�n y 
#pron�stico de la componente estacional, con la  regresi�n loess disponible en la funci�n R loess.as() de la librer�a fANCOVA, 
#para la estimaci�n y pron�stico de la tendencia. Combinando estos resultados mediante un modelo de descomposici�n aditiva o multiplicativa,
#obtiene la estimaci�n y el pron�stico total de una serie de tiempo, en la estrategia de validaci�n cruzada en la cual se ajusta con n obs.
#y se pron�stican los �ltimos m datos observados

##Sus argumentos son
 #1. serie.ajuste: Corresponde al conjunto de datos de la muestra de ajuste (serie recortada a n datos)
 #2. h: Longitud de los pron�sticos ex-post 
 #3. tipo.descomp: Para indicar el tipo de descomposici�n sus valores pueden ser "additive" en el caso aditivo � "multiplicative"
 #   en el caso multiplicativo
 #4. grado: Para especificar el grado del polinomio local, siendo 1 en el loess lineal � 2 en el loess cuadr�tico.
 #5. criterio: Para especificar el criterio de informaci�n a usar en la selecci�n del par�metro del suavizamiento loess, 
 #   puede ser "aicc" si el criterio es el AICC (Criterio de informaci�n de Akaike corregido), � "gcv"
 #   si el criterio es el GCV (criterio de validaci�n cruzada generalizada)

##Sus resultados: La funci�n produce un objeto tipo lista con los siguientes componentes
 #deltasi: Vector con valores estimados de los efectos estacionales en el orden de los periodos de un a�o calendario
 #alfa.optim: Escalar, correspondiendo al valor del par�metro de suavizamiento loess �ptimo seg�n criterio escogido
 #nep: Escalar, es el n�mero equivalente de par�metros de la regresi�n loess realizada sobre la serie desestacionalizada
 #de acuerdo al tipo de descomposici�n
 #p: Escalar, es el n�mero aproximado de p�rametros del ajuste combinando el fitro de la descomposici�n
 #en la estimaci�n de la estacionalidad, con la regresi�n loess en la estimaci�n de la tendencia
 #St: Objeto tipo ts (serie de tiempo) que corresponde a la estimaci�n de la componente estacional en los periodos de ajuste,
 #de acuerdo al tipo de descomposici�n realizada
 #Tt: Objeto tipo ts (serie de tiempo) que corresponde a la estimaci�n de la componente de tendencia en los periodos de ajuste
 #ytd: Objeto tipo ts (serie de tiempo) que corresponde a la versi�n desestacionalizada de la serie en lo periodos de ajuste,
 #de acuerdo al tipo de descomposici�n
 #fitted: Objeto tipo ts (serie de tiempo) que corresponde a la estimaci�n de la serie (juntando tendencia y estacionalidad aditiva o
 #multiplicativamente) en lo periodos de ajuste.
 #tablapron: Objetivo tipo mts (serie de tiempo m�ltiple) con tres columnas: Pron_Tt (pron�stico de la tendencia), 
 #Pron_St (pron�stico de la estacionalidad) y Pron_serie (pron�stico de la serie al combinar aditiva o multiplicativamente los pron�sticos
 #de las componentes).
 #ytpron: Objeto tipo ts (serie de tiempo) que corresponde a los pron�sticos de la serie en lo periodos ex-post.
 #residuals: Objeto tipo ts (serie de tiempo) que corresponde a los residuos del ajuste calculados como la diferencia entre valores observados
 #menos valores estimados, en los periodos de ajuste.
 #MSE: Escalar, es una aproximaci�n al MSE del ajuste total de la serie en los periodos de ajuste.

Descomp.Loess=function(serie.ajuste,h,tipo.descomp,grado,criterio){
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
p1=nep1+s-1 #n�mero aproximado de par�metros en ajuste Modelo
Tt1=ts(fitted(ajusteLoess1),frequency=s,start=start(serie.ajuste))
Ttnuevo1=predict(loess(ytd~tajuste,span=alfa.optim1,degree=grado,control=loess.control(surface="direct")),data.frame(tajuste=tpron),se=FALSE)
Ttnuevo1=ts(Ttnuevo1,freq=s,start=start.pron)
if(tipo.descomp=="additive"){
ythat1=St+Tt1
ytpron1=Ttnuevo1+Stnuevo
}
if(tipo.descomp=="multiplicative"){
ythat1=St*Tt1
ytpron1=Ttnuevo1*Stnuevo
}
tablapron1=cbind(Pron_Tt=Ttnuevo1,Pron_St=Stnuevo,Pron_serie=ytpron1)
et1=serie.ajuste-ythat1
df1=length(serie.ajuste)-nep1+s-1 
MSE1=sum(et1^2)/df1 
result=list(deltasi=deltas_i,alfa.optim=alfa.optim1,nep=nep1,p=p1,St=St,Tt=Tt1,ytd=ytd,fitted=ythat1,tablapron=tablapron1,ytpron=ytpron1,residuals=et1,MSE=MSE1)
result
}

