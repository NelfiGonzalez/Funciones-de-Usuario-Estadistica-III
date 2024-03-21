#La funcion predict_expo() esta orientada a la prediccion de la respuesta futura de forma puntual y por intervalos de prediccion
#Sobre modelos de regresion exponencial donde la funcion de regresion es de la forma g(x,beta)=exp(b0+b1*x1+b2*x2+...+bp*xp)=exp(t(x)%*%beta), 
#con x=c(1,x1,x2,...,xp) vector con predictores 
#beta vector de parametros (b0,b1,...,bp)
#también puede calcular I.C para la respuesta media en lugar de I.P
#Argumentos:
#modelo: Objeto tipo nls(), donde se ha guardado la estimacion del modelo
#new.data:       data.frame con los valores de los predictores en los puntos de prediccion. El nombre de las variables predictoras debe ser igual al usado
#                en el ajuste del modelo
#level:          Un valor en (0,1) para el nivel de confianza de los intervalos. Por defecto es 0.95
#interval:       Una cadena de caracteres correspondiedo a uno de lso siguientes casos: "none", "confidence", "prediction". Su valor por defecto es "none", es decir no se calculan intervalos
#                solo estimaciones o pronosticos puntuales. Para intervalos de confianza en la respuesta media use "confidence" y para intervalos de
#                prediccion use "prediction"

#Resultados
#La funcion produce un objeto data.frame con tres variables: 
#$forecast:      La estimacion o pronostico puntual de la respuesta
#$lower, $upper: Los limites inferior y superior de los intervalos pedidos.

predict_expo=function(modelo,new.data,level=0.95,interval = c("none", "confidence", "prediction")){
stopifnot("`modelo` debe ser un objeto clase `nls` " =class(modelo)=="nls","`interval`debe ser uno de los siguientes argumentos `none`, `confidence`, `prediction`"=interval%in%c("none", "confidence", "prediction"),"`level` debe ser valor en (0,1)"=level>0 & level<1)
#aprovechando normalidad del vector de parametros estimados y considerando g(beta)=exp(t(x0[i])%*%beta) con  x0[i] 
#el vector de predictores en i esima prediccion, incluye al 1 como primera entrada
x0=cbind(1,as.matrix(new.data)) #matriz con filas siendo los puntos de prediccion o valores de variables explicatorias en la prediccion
mse=summary(modelo)$sigma^2 #estimacion varianza del error del modelo de regresion exponencial
yhat=predict(modelo,newdata=new.data)
varbetas=as.matrix(vcov(modelo)) #matriz de varianzas covarianzas del vector de parametros estimados modelo exponencial
Sigmapred=(yhat^2)*diag(x0%*%varbetas%*%t(x0)) #aproximacion varianza de la respuesta predicha o estimada por metodo delta
if(interval=="none"){
res=predict(modelo,newdata=new.data)
}
if(interval=="confidence"){
#intervalos de confianza para la respuesta media usando aproximacion metodo delta en la varianza de estimacion puntual de la respuesta media
LI=yhat-qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(Sigmapred)
LS=yhat+qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(Sigmapred)
res=data.frame(estimate=yhat,lower=LI,upper=LS)
}
if(interval=="prediction"){
#intervalos de prediccion usando aproximacion metodo delta en la varianza de la prediccion puntual de la respuesta
LI=yhat-qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(mse+Sigmapred)
LS=yhat+qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(mse+Sigmapred)
res=data.frame(forecast=yhat,lower=LI,upper=LS)
}
res
}

