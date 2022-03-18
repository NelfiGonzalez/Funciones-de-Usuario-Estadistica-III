#Función de usuario para aplicar el Suavizamiento exponencial Holt-Winters
#aditivo y multiplicativo y obtener los valores suavizados de efectos 
#estacionales en el orden de los periodos del año calendario y no en el 
#orden de los "s" periodos de pronósticos posteriores al ajuste  
#Permite hacer el suavizamiento óptimo pero también predefinir los valores 
#de los parámetros de suavizamiento. 

##Sus argumentos son similares a los de la función HoltWinters() de R:
 #serie: Un objeto serie de tiempo
 #alpha: parametro del suavizamiento del nivel, valor entre 0 y 1
 #beta: parámetro del suavizamiento de la pendiente, valor entre 0 y 1
 #gamma: parámetro del suavizamiento de los efectos estaionales, valor entre 0 y 1
 #seasonal: cadena de caracteres, por defecto "additive" para el caso Holt-Winters aditivo. 
 #Para el caso multiplicativo se define como "multiplicative"
 #optim.start: Vector con componentes con nombre alpha, beta y gamma que contienen los valores iniciales 
 #para el optimizador. Sólo se deben especificar los valores necesarios. Ignorado en el caso de un parámetro.
 #optim.control: Lista opcional con parámetros de control adicionales pasados ​​a optim si se usa. 
 #Ignorado en el caso de un parámetro.

##Sus resultados: La función produce un objeto tipo lista con los siguientes componentes
 #suaviza: Una lista con los objetos que genera ordinariamente la función HoltWinters()
 #coefi: un data.frame con una sola columna de longitud s+2, que guarda exclusivamente los valores al final del
 #suavizamiento del nivel, la pendiente y los efectos estacionales (estos últimos en el orden de los periodos del año
 #calendario).

##Nota: Para obtener valores como residuos, valores ajustados, predicciones, etc. debe operarse sobre la componente de la lista $suaviza

SuavizamientoEstacional=function(serie,alpha=NULL,beta=NULL,gamma=NULL,seasonal="additive",optim.start=c(alpha = 0.3, beta = 0.1, gamma = 0.1),optim.control = list()){
suaviza=HoltWinters(serie,alpha=alpha,beta=beta,gamma=gamma,seasonal=seasonal,optim.start=optim.start,optim.control=optim.control)
s=frequency(suaviza$x)
if(end(suaviza$x)[2]<s){
estacionini=end(suaviza$x)[2]+1
}
if(end(suaviza$x)[2]==s){
estacionini=1
}

if(estacionini==1){
efectossuav=suaviza$coef[-c(1,2)]
}
if(estacionini!=1){
j=estacionini;efectossuav=c(suaviza$coef[-c(1,2)][(s-j+2):s],suaviza$coef[-c(1,2)][1:(s-j+1)])
}
if(s==12){
names(efectossuav)=c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
}
if(s==4){
names(efectossuav)=c("s1","s2","s3","s4")
}
pars=data.frame(c(suaviza$alpha,suaviza$beta,suaviza$gamma))
names(pars)=""
coefi=data.frame(c(suaviza$coef[1:2],efectossuav))
names(coefi)=""
cat("\n")
cat("Call:")
cat("\n")
print(suaviza$call)
cat("\n")
cat("Smoothing parameters:")
print(pars)
cat("\n")
cat("Coefficients:")
print(coefi)
result=list(suaviza=suaviza,coefi=coefi)
result
}



