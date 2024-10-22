#--------------------------------------------------------------------------------------------------
#Funcion usuario exp.crit.inf.resid() para calcular AIC y BIC version exp(C_n^*(p))
##Argumentos:
 #residuales: un vector numerico que debe corresponder a la diferencia Yt-Ythat, donde Yt es
 #el valor observado de la serie en su escala original y Ythat su valor estimado en la misma escala
 #NOTA: en modelos con transformacion de la serie, Yt-Ythat son los seudoresiduos y en modelos con 
 #ajuste sin transformar a la serie esta diferencia corresponde a residuos de ajuste.
 #n.par: Es el numero de parametros atribuido al modelo ajustado

##Resultados arrojados por esta funcion: un vector con los valores del numero de parametros (p), el AIC y el BIC.

exp.crit.inf.resid=function(residuales,n.par){
#Calcula AIC
AIC=exp(log(mean(residuales^2))+2*n.par/length(residuales))
#Calcula BIC
BIC=exp(log(mean(residuales^2))+n.par*log(length(residuales))/length(residuales))
crit=list(p=n.par,AIC=AIC,BIC=BIC)
unlist(crit)
}

#--------------------------------------------------------------------------------------------------
#Funcion para calcular la amplitud media y la cobertura (%) de "m">1 Intervalos de Prediccion (I.P)
##Argumentos:
  #real: serie de tiempo o vector numerico con los valores observados en los periodos de pronostico
  #LIP, LSP: Vectores numericos con los valores de los limites inferior y superior, respectivamente,   #de los intervalos de prediccion

##Resultados arrojados por esta funcion: un vector con los valores de la amplitud promedio y la cobertura (en porcentaje) de los "m" I.P  

amplitud.cobertura=function(real,LIP,LSP){
a=LSP-LIP
am=mean(a)
I=ifelse(real>=LIP & real<=LSP,1,0)
p=mean(I)*100
res=list(Amplitud=am,"Cobertura(%)"=p)
unlist(res)
}

#-------------------------------------------------------------------------------------------------------------------------------------
##Interval Scoring,
#------------------------------------------------------------------------------------------------------------------------------------- 
#ver GNEITING, T. and RAFTERY, A. E. (2007), "Strictly Proper Scoring Rules, Prediction and Estimation", JASA, vol. 102(477), pp. 359-378
#It is a Winkler Score which  is defined as the length of the interval plus a penalty if the observation is outside the interval
#For observations that fall within the interval, the Winkler score is simply the length of the interval. 
#Thus, low scores are associated with narrow intervals. 
#However, if the observation falls outside the interval, the penalty applies, 
#with the penalty proportional to how far the observation is outside the interval.
#-------------------------------------------------------------------------------------------------------------------------------------
##Argumentos:
# real: un valor o vector de valores observados de la serie en los periodos de pronostico
# LimInf,LimSup: valores o vectores de valores de los limites inferior y superior de los intervalos de prediccion correspondientes
# alpha: valor correspondiente al I.P (1-alpha)100%

## Valor arrojado: La funcion produce un numero real positivo correspondiente al Score del intervalo (en el caso de un solo IP) 
#  o al promedio de los Scores de los intervalos de prediccion (en el caso de m IPs)
#-------------------------------------------------------------------------------------------------------------------------------------

IntervalScore=function(real,LimInf,LimSup,alpha){
S.alpha=(LimSup-LimInf)+(2/alpha)*(LimInf-real)*ifelse(real<LimInf,1,0)+(2/alpha)*(real-LimSup)*ifelse(real>LimSup,1,0)
S.alpha.mean=mean(S.alpha)
S.alpha.mean
}


