#Funci�n de usuario Filtro.lineal() para el ajuste y pron�sticos a trav�s de filtros lineales basados en medias m�viles
#simples unilaterales o bilaterales, usando la estrategia circular
#Sus argumentos son
#serie:     Para especificar la serie a ajustar, debe ser un objeto tipo ts()
#pesos.wi:  para especificar el vector de pesos de las medias m�viles para la estimaci�n del nivel.
#           para medias m�viles simples unilaterales de par�metro m debe ser de longitud m+1, para i=0,1,...,m, y para bilaterales
#           debe ser de longitud 2m+1, para i=-m,...,m. 
#tipo:      Escalar, puede ser 1 � 2, seg�n si se desea medias m�viles unilaterales o bilaterales
#           para la estimaci�n del nivel
#h:         Escalar, para especificar el n�mero de pron�sticos despu�s del �ltimo tiempo de ajuste
#level:     El nivel de confianza para los intervalos de pron�stico, por defecto es 0.95

#Resultados: la funci�n produce un objeto tipo lista con los siguientes componentes
#nivel:     Serie de tiempo de las medias m�viles calculadas con la estrategia circular
#fitted:    objeto ts() con los valores ajustados de la serie
#residuals: objeto ts() con los valores de la serie de tiempo de los residuos del ajuste
#forecast:  objeto mts(), con los valores de los pron�sticos puntuales y por I.P para los h periodos
#           despu�s del ajuste
#RMSE:      estimaci�n de la ra�z cuadrada del MSE del ajuste.

Filtro.lineal=function(serie,pesos.wi,tipo,h,level=0.95){
nivel=filter(serie,filter=pesos.wi,method="convolution",sides=tipo,circular=TRUE)
s=frequency(serie)
ythat=ts(append(NA,nivel[-length(nivel)]),freq=s,start=start(nivel))
residuos=serie-ythat
RMSE=accuracy(ythat,serie)[2]
tpron=(length(serie)+1):(length(serie)+h)
if(end(serie)[2]<s){
start.pron=c(end(serie)[1],end(serie)[2]+1)
}
if(end(serie)[2]==s){
start.pron=c(end(serie)[1]+1,1)
}
ytpron=ts(rep(nivel[length(nivel)],h),frequency=s,start=start.pron)
df=length(serie)-2
lwr=ytpron-qt((1-level)/2,df=df,lower.tail=FALSE)*RMSE*sqrt(1+sum(pesos.wi)^2)
upr=ytpron+qt((1-level)/2,df=df,lower.tail=FALSE)*RMSE*sqrt(1+sum(pesos.wi)^2)
tablapron=ts(data.frame("Point Forecast"=ytpron,lwr,upr),freq=s,start=start.pron)
result=list(nivel=nivel,fitted=ythat,residuals=residuos,forecast=tablapron,RMSE=RMSE)
result
}

