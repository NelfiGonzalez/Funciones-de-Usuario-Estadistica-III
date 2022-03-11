# Funciones de Usuario Estadística-III 
## [Universidad Nacional de Colombia Sede Medellín](https://medellin.unal.edu.co/)
Funciones **creadas** para facilitar algunas rutinas R aplicadas en varios modelos vistos en el curso **Estadística III** de la Escuela de Estadística de la Facultad de Ciencias.
### Lista de funciones disponibles
* En archivo Funciones-Criterios.Informacion-Calidad.Intervalos.R se encuentra la función **exp.crit.inf.resid()** para calcular AIC y BIC, y **amplitud.cobertura()** para calcular la amplitud media y la cobertura (%) de "m">1 Intervalos de Predicción (I.P).
* En archivo Funcion-regexponencial.R se encuentra la función **regexponencial()** para el ajuste de modelos exponenciales polinomiales y exponenciales polinomiales estacionales
* En archivo Funcion-Filtro.lineal.R se encuentra la función **Filtro.lineal()** ara el ajuste y pronósticos a través de filtros lineales basados en medias móviles simples unilaterales o bilaterales, usando la estrategia circular.
* En archivo Funcion-Loess.Optimo.R se encuentra la función **Loess.Optimo()**, diseñada para ajustes y pronósticos de una serie no estacional, mediante loess óptimo escogiendo con criterio AICC o GCV el valor del parámetro de suavizamiento loess. Esta función opera con la librería fANCOVA de la cual hace uso de la función loess.as(), como también de la función loess() y predict().
* En archivo Funcion-SuavizamientoEstacional.R se encuentra la función **SuavizamientoEstacional()** Suavizamiento exponencial Holt-Winters aditivo y multiplicativo y obtener los valores suavizados de efectos estacionales en el orden de los periodos del año calendario y no en el orden de los "s" periodos de pronósticos posteriores al ajuste. Permite hacer el suavizamiento óptimo pero también predefinir los valores de los parámetros de suavizamiento. 
* En archivo Funcion-Descomp.Loess.R se encuentra la función **Descomp.Loess()** para el ajuste y pronósticos, por la combinación del filtro de descomposición con loess; Requiere la librería fANCOVA. Combina el filtro de la descomposición clásica disponible en la función R decompose(), para la estimación y pronóstico de la componente estacional, con la  regresión loess disponible en la función R loess.as() de la librería fANCOVA, para la estimación y pronóstico de la tendencia. 
* En archivo Funciones-BP.LB.test-pruebaDW1.R se encuentran las funciones de usuario **BP.LB.test()** para para realizar las pruebas Ljung-Box o Box-Pierce, para m=6, 12, 18,..., [maxlag/6] donde maxlag es el máximo orden de rezagos hasta el que se desea aplicar el test. Esta función hace uso de la función R Box.test(), y la función de usuario **pruebaDW1()** para evaluar el test Durbin-Watson para autocorrelación de orden 1 en un MRLM, la cual usa la función R durbinWatsonTest() de la librería car, y proporciona simultáneamente los valores P correspondientes a los tests con hipótesis alternativa de autocorrelación de orden 1 positiva y negativa.
