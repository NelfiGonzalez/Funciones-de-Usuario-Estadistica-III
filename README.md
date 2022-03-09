# Funciones de Usuario Estadística-III 
## [Universidad Nacional de Colombia Sede Medellín](https://medellin.unal.edu.co/)
Funciones **creadas** para facilitar algunas rutinas en modelos del curso Estadística III de la Escuela de Estadística de la Facultad de Ciencias.
### Lista de funciones disponibles
* En archivo Funciones-Criterios.Informacion-Calidad.Intervalos.R se encuentra **exp.crit.inf.resid()** para calcular AIC y BIC, y **amplitud.cobertura()** para calcular la amplitud media y la cobertura (%) de "m">1 Intervalos de Predicción (I.P).
* En archivo Funcion-regexponencial.R se encuentra la función **regexponencial()** para el ajuste de modelos exponenciales polinomiales y exponenciales polinomiales estacionales
* En archivo Funcion-Filtro.lineal.R se encuentra la función **Filtro.lineal()** ara el ajuste y pronósticos a través de filtros lineales basados en medias móviles simples unilaterales o bilaterales, usando la estrategia circular.
* En archivo Funcion-Loess.Optimo.R se encuentra la función **Loess.Optimo()**, diseñada para ajustes y pronósticos de una serie no estacional, mediante loess óptimo escogiendo con criterio AICC o GCV el valor del parámetro de suavizamiento loess. Esta función opera con la librería fANCOVA de la cual hace uso de la función loess.as(), como también de la función loess() y predict().
* En archivo Funcion-SuavizamientoEstacional.R se encuentra la función **SuavizamientoEstacional()** Suavizamiento exponencial Holt-Winters aditivo y multiplicativo y obtener los valores suavizados de efectos estacionales en el orden de los periodos del año calendario y no en el orden de los "s" periodos de pronósticos posteriores al ajuste. Permite hacer el suavizamiento óptimo pero también predefinir los valores de los parámetros de suavizamiento. $\mu$
