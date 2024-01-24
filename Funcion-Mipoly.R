#Funcion para generar los predictores de los polinomios de modelos de regresion global.
#Sus argumentos son:
#tiempo: Se refiere al vector de valores con el indice de tiempo t.
#grado: Se refiere al grado p del polinomio deseado.
#Resultados: La funcion produce un data.frame cuyas variables son t y sus potencias t^j denominadas tj, para j=2,..., p, segun el grado polinomial p.

Mipoly=function(tiempo,grado){
if(grado>1){
poli=data.frame(poly(tiempo,degree=grado,raw=TRUE))
names(poli)=c("t",paste0("t",2:ncol(poli)))
}
if(grado==1){
poli=data.frame(t=tiempo)
}
poli
}
