
#FUNCION DE USUARIO "regexponencialv02" PARA EL AJUSTE DE MODELOS EXPONENCIALES
##Sus argumentos son:
 #respuesta: el vector de valores de la serie en el ajuste
 #data: un marco de datos con los predictores a usar
 #control: una lista opcional de configuraciones de control. Consulte nls.control() para conocer los nombres de los valores de control
 #configurables y su efecto.

##Sus resultados: Son los mismos predeterminados para los objetos tipo nls generados por la funcion nls().

regexponencialv02=function(respuesta,data,control=nls.control()){
yt=respuesta
model.aux=lm(log(yt)~.,data)
names.vars=c(1,names(data))
names.param=c("Intercept",paste0("param_",names(data)))
miformula=as.formula(paste("yt~",paste(paste("exp(",paste(names.param,names.vars,sep="*",collapse="+"),sep=""),")",sep="")))
coef0=as.list(coef(model.aux));names(coef0)=as.list(names.param)
modelo=nls(miformula,start=coef0,control=control,data=data)
modelo
}
