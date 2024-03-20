#La funcion interpdeltas() esta orientada a la intepretacion de los parametros estacionales en la modelacion global con indicadoras. 
#Permite extraer y construir graficos de los efectos estacionales en modelos de regresion global estacionales que usan indicadoras
#con nivel de referencia el "s" esimo periodo calendario. Particularmente en modelos polinomiales estacionales con error RB y con error ARMA, 
#modelos log polinomiales estacionales con error RB y con error ARMA y modelos exponenciales polinomiales con error RB.
#En los modelos polinomiales estacionales bien sea con error RB o ARMA, extrae los valores estimados de los parametros asociados a las 
#variables indicadoras y crea la grafica de estos valores vs. periodo calendario asignando valor de cero en el ultimo periodo
#En los modelos log polinomiales estacionales con error RB o ARMA y en los modelos exponenciales polinomiales estacionales con error RB, 
#extrae los valores estimados de los parametros asociados a las variables indicadoras, pero es necesario indicar que estos modelos son 
#multiplicativos para que exponencie estas estimaciones, las multiplique por 100 y las grafique contra el periodo calendario asignando 
#el valor de 100 al ultimo periodo.

#ARGUMENTOS
#modelo:                       El nombre del objeto R con el modelo ajustado. Objetos tipo lm(), nls() o Arima()
#ordenp, ordenq,ordenP,ordenQ: Argumentos de valor entero, para especificar los valores de los ordenes p, q, P, Q de modelos ARMA(p,q)(P,Q)[s], 
#                              por defecto estan fijos en 0, es decir, asumiendo error RB. Si el error es un AR(p) el usuario debe especificar
#                              el valor de p, si el error es MA(q) debe especificar el valor de q, si es ARMA(p,q) debe especificar p y q, etc.
#gradopoly:                    Se refiere al grado del polinomio global, siempre debe ser especificado.
#aditivo:                      Argumento logico (TRUE, FALSE), por defecto es TRUE, este valor debe usarse en modelos polinomiales estacionales con indicadoras
#                              con error RB o error ARMA. En caso de tratarse de modelos log polinomiales estacionales con error RB o error 
#                              ARMA, o modelos exponenciales polinomiales estacionales con error RB, este argumento debe especificarse igual a FALSE.
#plotit:                       Argumento logico (TRUE, FALSE), por defecto es TRUE indicando que se desea la grafica de los efectos estacionales  
#                              para su interpretacion en la escala de la serie, segun si el modelo es aditivo o multiplicativo                                               


#RESULTADOS
#Mientras argumento plotit=TRUE, despliega la grafica correspondiente y en la consola R el vector de parametros estacionales estimados (modelos aditivos) o de
#los valores exponenciados de estos parametros y multipliados por 100% (modelos multiplicativos). Ademas la funcion crea un objeto lista con los siguientes componentes
#periodo:      vector de valores enteros de 1 a "s" con "s" la longitud anual
#deltas:       vector con las estimaciones de los parametros asociados a las indicadoras y adiciona el valor de cero al final. Este vector aparece cuando el argumento
#              aditivo=TRUE.
#expdeltas100: vector con los valores exponenciados de las estimaciones de los parametros asociados a las indicadoras, multiplicados por 100 y adiciona 
#              el valor de 100 al final. Este vector aparece cuando el argumento aditivo=FALSE.

interpdeltas=function(modelo,ordenp=0L,ordenq=0L,ordenP=0L,ordenQ=0L,gradopoly,aditivo=TRUE,plotit=TRUE){
#stopifnot("`modelo` debe ser un objeto clase `lm`,`nls` o `Arima` " =class(modelo)=="lm"|class(modelo)=="nls"|class(modelo)=="Arima")
aux=ordenp+ordenq+ordenP+ordenQ+gradopoly+1
if(aditivo==TRUE){
deltas=c(coef(modelo)[-c(1:aux)],0)
names(deltas)=paste0("delta",1:length(deltas))
print(deltas[1:(length(deltas)-1)])
if(plotit==TRUE){
plot(1:length(deltas),deltas,type="b",pch=1,lwd=3,xlab="periodo calendario",ylab=expression(widehat(delta)[i]),xaxt="n")
axis(1,at=1:length(deltas),1:length(deltas))
}
periodo=1:length(deltas)
res=list(periodo=periodo,deltas=deltas)
}
if(aditivo==FALSE){
expdeltas100=exp(c(coef(modelo)[-c(1:aux)],0))*100
names(expdeltas100)=paste0("exp","(delta",1:length(expdeltas100),")*100%")
print(expdeltas100[1:(length(expdeltas100)-1)])
if(plotit==TRUE){
plot(1:length(expdeltas100),expdeltas100,type="b",pch=1,lwd=3,xlab="periodo calendario",ylab=expression(exp(widehat(delta)[i])*100),xaxt="n")
axis(1,at=1:length(expdeltas100),labels=1:length(expdeltas100))
}
periodo=1:length(expdeltas100)
res=list(periodo=periodo,expdeltas100=expdeltas100)
}
}

