#Funcion para generar las variables trigonometricas a usar en modelos de regresion global
#para cualesquier frecuencias e indexando como se desee cada una de ellas.
#Sus argumentos son:
#tiempo:      Para ingresar el vector de valores con el indice de tiempo para los cuales se debe calcular las funciones trigonometricas. 
#Frecuencias: Para ingresar el vector con las frecuencias Fj correspondientes a las ondas sinusoidales armonicas a considerar.
#indicej:     Para ingresar el objeto vectorial indicando el indice j a usar como identificador para cada frecuencia en el vector 'Frecuencias',
#             de manera que a la componente 'Frecuencias[i]' le es asignado el indice identificador 'indicej[i]'.

#Resultados:  La funcion produce un data.frame cuyas variables son denotadas 'sinFj', 'cosFj', donde el sufijo j es segun el indice 
#             identificador asignado en el argumento 'indicej' a las frecuencias especificadas en el argumento 'Frecuencias'.   
#             Si en el conjunto de frecuencias se encuentra 1/2, para esta solo resulta la componente coseno, que corresponde a cos(pi*tiempo).
#             Los pares de variables 'sinj', 'cosj', aparecen en el data.frame en el orden de menor a mayor frecuencia ingresada en el argumento
#             'Frecuencias', pero conservando el sufijo j que se le haya asignado en el argumento 'indicej'.

Mytrigon=function(tiempo,Frecuencias,indicej){
F=sort(Frecuencias)
if(max(F)!=0.5){
K=length(F)
trig=matrix(,ncol=2*K,nrow=length(tiempo))
for(i in 1:K){
trig[,2*i-1]=sin(2*F[i]*pi*tiempo)
trig[,2*i]=cos(2*F[i]*pi*tiempo)
}
trig=data.frame(trig)
J=indicej[order(Frecuencias)]
names(trig)=paste0(rep(c("sinF","cosF"),times=K),rep(J,each=2))
}
if(max(F)==0.5){
K=length(F)-1
trig=matrix(,ncol=2*K+1,nrow=length(tiempo))
for(i in 1:K){
trig[,2*i-1]=sin(2*F[i]*pi*tiempo)
trig[,2*i]=cos(2*F[i]*pi*tiempo)
}
trig[,2*K+1]=cos(2*F[K+1]*pi*tiempo)
trig=data.frame(trig)
J=indicej[order(Frecuencias)]
names(trig)=c(paste0(rep(c("sinF","cosF"),times=K),rep(J[1:K],each=2)),paste0("cosF",J[K+1]))
}
trig
}

