# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO: El objetivo de este script es filtrar la cantidad de genes con las que estamos trabajando ya que no todos tienen relevancia. 
#Vamos a liminar aquellos que tengan mas de un 70% de los pacientes con mas de 10 cuentas. Luego aplicaremos un logaritmo a la muestras
#para estabilizar la varianza. 

############################################################################################################################################

#PRE-PROCESAMIENTO LUAD: 

##########################################################################################################################################

#Cargamos la matriz generada en el paso anterior y establecemos el directorio de trabajo; 
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/2_ScriptsPreProcesamiento/")
matrizLuad <- readRDS("MatrizCuentasLuad.rds")

#cargamos las librerias necesarias
library(tidyverse)

#Establecemos los elementos necesarios para generar un bucle que me va a poner un 0 donde sea menor de 10 y un 1 donde sea mayor: 
matrizCuentas <- matrix(0,nrow=nrow(matrizLuad),ncol=ncol(matrizLuad)) #Matriz de 1 y 0
i <- 1
j <- 1
filtro <- 10
start <- 1
end <- length(matrizLuad)
numeroColumnas <- seq(from=start,to=end,by=1) #Necesario para que me recorra toda la matriz
numeroFilas <- seq(from=start,to=nrow(matrizLuad),by=1) #Necesario para que me recorra toda la matriz

#Vamos a asginar un 0 o un 1 dependiendo de si ese valor de cuenta es mayor o menor que 10 
for (i in seq_along(numeroFilas)){ #Es la fila
  for (j in seq_along(numeroColumnas)){ #Es la columna
    cuenta <- matrizLuad[i,j] #Para indexar en el interior de un data.frame
    if (cuenta>=filtro){
      matrizCuentas[i,j] <- 1
    }else{
      matrizCuentas[i,j] <- 0
    }
  }
}
#Sumamos los valores de 0 y 1
sumaFilas <- rowSums(matrizCuentas)

#Filtramos por aquellos que tengan menos de un 70% de cuentas: 
indices <- which(sumaFilas>=49)

Luad_Filtrado <- matrizLuad[indices,]
saveRDS(Luad_Filtrado,"Luad_Filtrado.rds")

#Separamos los genes en un vector aparte

Luad_Filtrado_Logaritmo <- log(Luad_Filtrado+0.5)
Luad_Filtrados_Logaritmo_Escalado <- data.frame(scale(Luad_Filtrado_Logaritmo))
saveRDS(Luad_Filtrados_Logaritmo_Escalado, 'Luad_Escalado.rds')

############################################################################################################################################

#PRE-PROCESAMIENTO LUAD: 

##########################################################################################################################################

matrizLusc <- readRDS("MatrizCuentasLusc.rds")

#cargamos las librerias necesarias
library(tidyverse)

#Establecemos los elementos necesarios para generar un bucle que me va a poner un 0 donde sea menor de 10 y un 1 donde sea mayor: 
matrizCuentas <- matrix(0,nrow=nrow(matrizLusc),ncol=ncol(matrizLusc)) #Matriz de 1 y 0
i <- 1
j <- 1
filtro <- 10
start <- 1
end <- length(matrizLusc)
numeroColumnas <- seq(from=start,to=end,by=1) #Necesario para que me recorra toda la matriz
numeroFilas <- seq(from=start,to=nrow(matrizLusc),by=1) #Necesario para que me recorra toda la matriz

#Vamos a asginar un 0 o un 1 dependiendo de si ese valor de cuenta es mayor o menor que 10 
for (i in seq_along(numeroFilas)){ #Es la fila
  for (j in seq_along(numeroColumnas)){ #Es la columna
    cuenta <- matrizLusc[i,j] #Para indexar en el interior de un data.frame
    if (cuenta>=filtro){
      matrizCuentas[i,j] <- 1
    }else{
      matrizCuentas[i,j] <- 0
    }
  }
}
#Sumamos los valores de 0 y 1
sumaFilas <- rowSums(matrizCuentas)

#Filtramos por aquellos que tengan menos de un 70% de cuentas: 
indices <- which(sumaFilas>=385)

Lusc_Filtrado <- matrizLusc[indices,]
saveRDS(Lusc_Filtrado,"Lusc_Filtrado.rds")

#Separamos los genes en un vector aparte

Lusc_Filtrado_Logaritmo <- log(Lusc_Filtrado+0.5)
Lusc_Filtrados_Logaritmo_Escalado <- data.frame(scale(Lusc_Filtrado_Logaritmo))
saveRDS(Lusc_Filtrados_Logaritmo_Escalado, 'Lusc_Escalado.rds')
