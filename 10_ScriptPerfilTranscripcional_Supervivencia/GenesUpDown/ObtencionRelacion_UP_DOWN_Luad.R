# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:Script para obtener los data.frame de 0 y 1 para saber si son superiores a la medio o a la mediana de LUAD

##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)
#########################################################################################################

#······································ L U A D ···············································

Matriz_1 <- readRDS("GenesUpDown/Matrices_UP_DOWN_Luad_1.rds")
Matriz_2 <- readRDS("GenesUpDown/Matrices_UP_DOWN_Luad_2.rds")
Matriz_3 <- readRDS("GenesUpDown/Matrices_UP_DOWN_Luad_3.rds")

Matriz_1_UP <- t(Matriz_1[[1]])
Matriz_1_DOWN <- t(Matriz_1[[2]])
Matriz_2_UP <- t(Matriz_2[[1]])
Matriz_2_DOWN <- t(Matriz_2[[2]])
Matriz_3_UP <- t(Matriz_3[[1]])
Matriz_3_DOWN <- t(Matriz_3[[2]])

MatrizTraspuesta_1 <- list(Matriz_1_UP,Matriz_1_DOWN)
MatrizTraspuesta_2 <- list(Matriz_2_UP,Matriz_2_DOWN)
MatrizTraspuesta_3 <- list(Matriz_3_UP,Matriz_3_DOWN)

saveRDS(MatrizTraspuesta_1,"MatrizTranscripcion_UP_DOWN_Luad_1.rds")
saveRDS(MatrizTraspuesta_2,"MatrizTranscripcion_UP_DOWN_Luad_2.rds")
saveRDS(MatrizTraspuesta_3,"MatrizTranscripcion_UP_DOWN_Luad_3.rds")

Matrices <- list(Matriz_1_UP,Matriz_1_DOWN,Matriz_2_UP,Matriz_2_DOWN,Matriz_3_UP,Matriz_3_DOWN)
SumatorioMatrices <- list()

i <- 1
for(i in seq_along(Matrices)){
  Matriz <- Matrices[[i]]
  SumatorioMatriz <- t(rowSums(Matriz))
  SumatorioMatrices[[i]] <- SumatorioMatriz
}

names(SumatorioMatrices) <- c("1_UP","1_DOWN","2_UP","2_DOWN","3_UP","3_DOWN")
condiciones <- c("1_UP","1_DOWN","2_UP","2_DOWN","3_UP","3_DOWN")
saveRDS(SumatorioMatrices,"SumatorioMatrices_UP_DOWN_Luad.rds")

Media_Mediana <- data.frame(Media=rep(0,length(SumatorioMatrices)),
                            Mediana=rep(0,length(SumatorioMatrices)),
                            Condicion=condiciones)
i <- 1

for (i in seq_along(SumatorioMatrices)){
  Matriz <- SumatorioMatrices[[i]]
  Media <- mean(Matriz)
  Mediana <- median(Matriz)
  
  Media_Mediana[[i,1]] <- Media
  Media_Mediana[[i,2]] <- Mediana
}
Matriz <- Matriz 
Informacion_1_UP <- data.frame(Media=rep(0,length(Matriz)),
                            Mediana=rep(0,length(Matriz)))
Informacion_1_DOWN <- data.frame(Media=rep(0,length(Matriz)),
                               Mediana=rep(0,length(Matriz)))

rownames(Informacion_1_UP) <- rownames(Matriz_1_DOWN)
rownames(Informacion_1_DOWN) <- rownames(Matriz_1_DOWN)

Informacion_2_UP <- data.frame(Media=rep(0,length(Matriz)),
                            Mediana=rep(0,length(Matriz)))
Informacion_2_DOWN <- data.frame(Media=rep(0,length(Matriz)),
                               Mediana=rep(0,length(Matriz)))

rownames(Informacion_2_UP) <- rownames(Matriz_1_DOWN)
rownames(Informacion_2_DOWN) <- rownames(Matriz_1_DOWN)

Informacion_3_UP <- data.frame(Media=rep(0,length(Matriz)),
                            Mediana=rep(0,length(Matriz)))
Informacion_3_DOWN <- data.frame(Media=rep(0,length(Matriz)),
                               Mediana=rep(0,length(Matriz)))

rownames(Informacion_3_UP) <- rownames(Matriz_1_DOWN)
rownames(Informacion_3_DOWN) <- rownames(Matriz_1_DOWN)


InformacionTotal <- list(Informacion_1_UP,Informacion_1_DOWN,Informacion_2_UP,Informacion_2_DOWN,Informacion_3_UP,Informacion_3_DOWN)
names(InformacionTotal) <- condiciones

i <- 1
for (i in seq_along(SumatorioMatrices)){
  Sumatorio <- SumatorioMatrices[[i]]
  Media <- Media_Mediana[i,1]
  Mediana <- Media_Mediana[i,2]
  Informacion <- InformacionTotal[[i]]
  
  j <- 1
  
  for (j in seq_along(Sumatorio)){
    Valor <- Sumatorio[j]
    
    if (Valor > Media){Informacion[j,1] <- 1}
    else{Informacion[j,1] <- 0}
    
    if(Valor > Mediana){Informacion[j,2] <- 1}
    else(Informacion[j,2] <- 0)
  }
  
  InformacionTotal[[i]] <- Informacion
}

saveRDS(InformacionTotal,"Informacion_UP_DOWN_SinClinica_Luad.rds")
