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

