#Funcion de usuario regexpo.ErrorARMAv02() para el ajuste y pronostico de modelos de regresion exponenciales con error ARMA, hace un ajuste aproximado en el cual estima primero los parametros de las componentes estructurales
#(tendencia y estacionalidad) por minimos cuadrados no lineales, luego, ajusta a los residuos estructurales de este modelo, un modelo ARMA estacionario de media cero, mediante la funcion Arima() de la libreria forecast. 
#Los valores estimados de la serie son aproximados sumando las estimaciones del ajuste de las componentes estructurales y las estimaciones obtenidas del modelo ARMA sobre los 
# residuos estructurales del modelo preliminar. El anterior procedimiento tambien es aplicado para la construccion de pronosticos puntuales (no hay posibilidad de realizar pronosticos por intervalo). 

#Los argumentos son:

# respuesta:     Un objeto ts con los datos a ser ajustados.
#data, newdata:  Objetos tipo data.frame, cuyas columnas corresponden a los valores de las variables predictoras del modelo en el ajuste y en el pronostico, respectivamente.
# control:       una lista opcional de configuraciones de control para la funcion nls(). Consulte nls.control() para conocer los nombres de los valores de control
#                configurables y su efecto.
# order:         Una especificacion de la parte no estacional del modelo ARIMA: los tres componentes (p, d, q) son el orden AR, el grado de diferenciacion y el orden MA, respectivamente.
# seasonal:      Una especificacion de la parte estacional del modelo ARIMA, mas el periodo (que por defecto es frequency(respuesta)). Este argumente debe ser dado como una lista 
#                con el orden de los componentes y el periodo, pero una especificacion de solo un vector numerico de longitud 3 se convertira en una lista adecuada con las especificaciones como el orden
# fixed:         Vector numerico opcional de la misma longitud que la suma de los ordenes de la ecuacion del ARMA estacionario de media cero a ajustar sobre residuos estructurales del ajuste exponencial. 
#                Ver ayuda de funcion arima() sobre detalles de este argumento.
# method:        Metodo de ajuste: maxima verosimilitud ("ML"), Maxima verosilmilitud combinada con minimos cuadrados condicionales ("CSS-ML") o minimos de cuadrados condicionales ("CSS"). 
#                El valor predeterminado (a menos que existan valores faltantes) es "CSS-ML", en el cual se usan minimos cuadrados condicionales para encontrar los valores iniciales y 
#                luego aplica maxima verosimilitud. 
# optim.method:  El valor pasado al argumento  ‘method’ para ‘optim’ en la funcion Arima(). Por defecto es "BFGS".
# optim.control: Una lista de parametros de control para ‘optim’, usandos en la funcion Arima().

# Resultados:    La funcion produce una lista con las siguientes componentes
# coefficients:  Matriz con la tabla de parametros estimados, sus errores estandar, estadistico T0 y valor P asociado. Tenga en cuenta que no hay una estimacion conjunta de los parametros
#                de regresion de la funcion exponencial y de los parametros del modelo ARMA del error estructural, de modo que los errores estandar provienen del ajuste separado de las dos 
#                estructuras, pero los valores p son calculados bajo una distribucion t-student cuyos grados de libertad son df=n - total de parametros. el objeto tabla puede ser obtenido
#                con la funcion coef() aplicada al objeto R donde se guarde el resultado de la funcion regexpo.ErrorARMAv02().
# fitted:        Objeto tipo ts con los valores ajustados de la respuesta. Estos valores pueden ser accesados mediante la funcion fitted() 
#                sobre el objeto R donde se guarde el resultado de la funcion regexpo.ErrorARMAv02().
# residuals:     Objeto tipo ts con los residuos del ajuste de la respuesta. Estos valores pueden ser accesados mediante la funcion residuals() 
#                sobre el objeto R donde se guarde el resultado de la funcion regexpo.ErrorARMAv02(). 
# sigma2:        numerico con la estimacion de la varianza de las innovaciones del modelo ARMA definido para el error estructural del modelo de regresion. Puede ser accesado de la siguiente manera
#                nombre_objeto$sigma2, donde 'nombre_objeto' es el nombre del objeto R donde se guarda el resultado de la funcion regexpo.ErrorARMAv02().
# forecast:      Objeto tipo ts con los pronosticos puntuales de la respuesta para h=ncol(newdata) periodos despues del ajuste.


regexpo.ErrorARMAv02=function(respuesta,data,newdata,control=nls.control(),order= c(0L, 0L, 0L),seasonal=list(order = c(0L, 0L, 0L),period=NA),fixed=NULL,method=c("CSS-ML", "ML", "CSS"),optim.method="BFGS",optim.control = list()){
yt=respuesta
model.aux=lm(log(yt)~.,data)
names.vars=c(1,names(data))
names.param=c("Intercept",paste0("param_",names(data)))
miformula=as.formula(paste("yt~",paste(paste("exp(",paste(names.param,names.vars,sep="*",collapse="+"),sep=""),")",sep="")))
coef0=as.list(coef(model.aux));names(coef0)=as.list(names.param)
modelo=nls(miformula,start=coef0,control=control,data=data)
serie.et=ts(residuals(modelo),freq=frequency(yt),start=start(yt))
ciclos=Arima(serie.et,order=order,seasonal=seasonal,include.mean=F,fixed=fixed,method=method,optim.method=optim.method,optim.control=optim.control)
p1=length(coef(modelo)[coef(modelo)!=0])
p2=p1+length(coef(ciclos)[coef(ciclos)!=0])
df=length(yt)-p2
tabla=rbind(coeftest(ciclos,df=df),coeftest(modelo,df=df))[,1:4]
ythat=ts(fitted(modelo),freq=frequency(yt),start=start(yt))+fitted(ciclos)
sigma2=ciclos$sigma2
s=tsp(yt)[3]
tajuste=1:length(yt)
h=nrow(newdata)
tpron=(length(yt)+1):(length(yt)+h)
if(end(yt)[2]<s){
start.pron=c(end(yt)[1],end(yt)[2]+1)
}
if(end(yt)[2]==s){
start.pron=c(end(yt)[1]+1,1)
}
resid=residuals(ciclos)
ytpron=ts(predict(modelo,newdata=newdata,interval="prediction",level=0.95),freq=frequency(yt),start=start.pron)+ts(forecast(ciclos,h=h)$mean,freq=frequency(yt),start=start.pron)
result=list(coefficients=tabla,fitted=ythat,residuals=resid,sigma2=sigma2,p=p2,forecast=ytpron)
result
}
