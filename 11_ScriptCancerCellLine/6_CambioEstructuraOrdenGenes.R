# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 27-05-2021

#OBJETIVO:Script para cambiar el orden del data.frame con el orden de los genes 

################################################################################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)

Genes <- readRDS("OrdenGenes_Lusc.rds")
################################################################################################################################################

#PRUEBA DE COMO QUERIA PLANTEARLO: Empezamos solo con uno haciendo poco a poco los pasos para luego unirlo. 
G1 <- (Genes[[1]])
G1 <- as.data.frame(G1)

#Necesito dos columnas, en una la linea celular y en la otra los genes. 
i <- 1
Union <- list()
for(i in seq(1,nrow(G1),1)){
  Linea <- t(G1[i,])
  NombreLinea <- rep(colnames(Linea),nrow(Linea))
  
  Union[[i]] <- as.data.frame(cbind(NombreLinea,Linea))
  colnames(Union[[i]])[2] <- "Gen"
  names(Union)[[i]] <- colnames(Linea)
}

i <- 1
Linea1 <- Union[[i]]
Linea2 <- Union[[i+1]]
UnionLineas <- rbind(Linea1,Linea2)

for(i in seq(3,length(Union),1)){
  Linea <- Union[[i]]
  
  UnionLineas <- rbind(UnionLineas,Linea)

}
################################################################################################################################################

# V e r s i o n  D e f i n i t i va: solo hacemos la union en un unico data.frame sin anotaciones. 

i <- 1 #Indice que va cambiando de grupo
UnionFinal <- list()
for (i in seq_along(Genes)){
  Gen <- as.data.frame(Genes[[i]])
  
  j <- 1 #Indice que va cambiando de linea celular
  Union <- list() #Vacia porque en cada grupo se empieza de nuevo. 
  
  #En este bucle generamos un data.frame por cada linea celular en el que: 
    #1) Repetimo el nombre de la linea el numero de genes que tiene ese grupo
    #2) Juntamos la repeticion del nombre de la linea con el gen
    #3) Lo guardamos en una lista que tendra 23 objetos
  
  for (j in seq(1,nrow(Gen),1)){
    Linea <- t(Gen[j,]) #Selecionamos la fila correspondiente y al hacer la traspuesta nos quedamos con la columna
    NombreLinea <- rep(colnames(Linea),nrow(Linea)) #Repetimos el nombre de la linea el numero de genes que esta contiene. 
    #El numero de genes viene especificamos por el numero de filas por haber hecho la traspuesta. 
    
    Union[[j]] <- as.data.frame(cbind(NombreLinea,Linea))
    colnames(Union[[j]])[2] <- "Gen"
    names(Union)[[j]] <- colnames(Linea)
  }
  
  #En esta parte del bucle unimos todo en el mismo data.frame.
  #Tenemos que unir por bloques, es decir los X genes de cada una de las lineas uno detras del otro 
  
  #Para ello: 
    #1) Primero unimos los dos primeros bloques fuera del bucle 
    #2) Dentro del bucle unimos apartir del tercer elemento hasta el final 
  
  h <- 1 #Indice para ir uniendo las lineas
  Linea1 <- Union[[h]]
  Linea2 <- Union[[h+1]]
  UnionLineas <- rbind(Linea1,Linea2)
  Longitud <- length(Union)
  
  Secuencia <- seq(3,length(Union),1) #Los valores que tiene que coger H
  
  for(h in seq_along(Secuencia)){
    Index <- Secuencia[h] #Lo especifico aqui primero para asegurarme que realmente empieza en el numero 3
    LineaNueva <- Union[[Index]] #Cogemos la linea a partir de las 3
    
    UnionLineas <- rbind(UnionLineas,LineaNueva) #Unimos lo que ya teniamos con lo nuevo
  }
  
  #Guardamos cada uno de los 3 data.frames ya que habra uno por grupo
  UnionFinal[[i]] <- UnionLineas
  
}

saveRDS(UnionFinal,"EstructuraLonger_OrdenGenes.rds")

################################################################################################################################################
#Es la lista donde tengo los data.frame ordenados con menos columnas
Union <- readRDS("EstructuraLonger_OrdenGenes.rds")
#Es donde tengo el data.frame donde cada linea es una lineaCelular
Genes <- readRDS("OrdenGenes_Lusc.rds")

# Ahora vamos ha hacer una comprobacion de que realmente se me han unido de forma correcta:
  #1) En cada grupo tendria que haber 23 lineas celulares unicas
  #2) Compobar que el orden del los genes es el mismo 

i <- 1 #Indice para seleccionar grupo
#El numero de lineas en todos es el mismo: 
NumeroLineasCelulares <- nrow(Genes[[1]])

for ( i in seq_along(Union)){
  #Aqui seleccionamos el data.frame donde tenemos dos columnas con todos los genes de ese grupo
  GenesColumna <- Union[[i]] 
  #Aqui selecionamos el data.frame donde cada linea es una linea celular, por lo tanto los genes estan en filas 
  GenesFila <- Genes[[i]]
  #Sacamos el numero de genes que tiene ese grupo en concreto
  NumeroGenes <- ncol(GenesFila)
  
  # P R I M E R A  C O M P R O B A C I O N: Comprobar que las los genes son iguales
  j <- 1 #Para las columnas
  h <- 1 #Para las filas
  #Para ir yendo por bloques en el data.frame de las columnas
  SecuenciasColumna <- seq(1,nrow(GenesColumna),NumeroGenes)
  #Para ir yendo por filas en el data.frame de las filas
  SecuenciasFila <- seq(1,NumeroLineasCelulares,1)
  
  for(j in seq_along(GenesColumna)){
    #Seleccionamos los indices inciales y finales para luego poder indexar
    IndiceInicial <- SecuenciasColumna[j]
    IndiceFinal <- SecuenciasColumna[j+1]-1
    
    #Indexamos y me quedo solo con la segunda columna que es la que tiene los genes para que me haga un vector
    Genes_Columna <- GenesColumna[(IndiceInicial:IndiceFinal),]
    Genes_Columna <- Genes_Columna[,2]
    
    #Indexamos en el data.frame de fila para obtener el vector con los genes
    IndiceFila <- SecuenciasFila[h]
    Genes_Fila <- t(GenesFila[h,])
    
    #Sumamos para poder seguir indexando ya que es un bucle con dos indices
    h <- h+1
    
    #Comprobamos que los genes sean iguales 
    Comprobacion <- all(Genes_Columna==Genes_Fila)
    
    #Ponemos la comprobacion en el caso de que no sea igual. 
    if(Comprobacion==FALSE){print(paste("La linea",colnames(Genes_Fila),"no coinciden los genes en el grupo",i))}
  }
  
  # S E G U N D A  C O M P R O B A C I O N: que las lineas se repiten X veces. 
  k <- 1
  for( k in seq(SecuenciasFila)){
    #Seleccionamos la linea, el nombre
    IndiceFila <- SecuenciasFila[k]
    Linea <- rownames(GenesFila)[IndiceFila]
    
    #Agrupamos por ese nombre
    Linea_Columna <- GenesColumna %>% filter(NombreLinea==Linea)
    #Vemos que tenga coincida con el numero de genes
    NumerodeLineasColumna <- nrow(Linea_Columna)
    
    #Ponemos la comprobacion para que nos la imprima
    if(NumerodeLineasColumna!=NumeroGenes){print(paste("La linea",Linea,"tiene mas valores que los que le corresponden"))}
  }
  
}




