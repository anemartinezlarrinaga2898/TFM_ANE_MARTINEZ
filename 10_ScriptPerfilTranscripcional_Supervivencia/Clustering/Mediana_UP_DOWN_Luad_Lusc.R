# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-05-2021

#OBJETIVO:Script para obtener si la expresion de dichos pacientes esta por encima o por debajo de la MEDIANA. 

##############################################################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)

#########################################################################################################

#······································ L U A D ···············································

# CALCULO DE LA MEDIANA: 

Matrices <- readRDS("Clustering/Grupos_UP_DOWN_Luad.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 2 DOWN","Grupo 3 UP","Grupo 3 DOWN")
SumatorioMatrices <- list()
MedianaMatrices <- vector()
i <- 1

for (i in seq_along(Matrices)){
  Matriz <- Matrices[[i]]
  Sumatorio <- rowSums(Matriz)
  Mediana <- median(Sumatorio)
  SumatorioMatrices[[i]] <- Sumatorio
  MedianaMatrices[i] <- Mediana
}

names(MedianaMatrices) <- Condiciones
names(SumatorioMatrices) <- Condiciones

saveRDS(SumatorioMatrices,"SumatorioMatrices_UP_DOWN_Luad.rds")
saveRDS(MedianaMatrices,"MedianaMatrices_UP_DOWN_Luad.rds")

#···············································································································

# ESTIMACION DE SI SE ENCUENTRA POR ENCIMA O POR DEBAJO: 

Matrices <- readRDS("Clustering/Grupos_UP_DOWN_Luad.rds")
SumatorioMatrices <- readRDS("Clustering/SumatorioMatrices_UP_DOWN_Luad.rds")
MedianaMatrices <- readRDS("Clustering/MedianaMatrices_UP_DOWN_Luad.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 2 DOWN","Grupo 3 UP","Grupo 3 DOWN")


Informacion_1_UP <- data.frame(Mediana=rep(0,nrow(Matrices[[1]])),
                               Grupo=rep(0,nrow(Matrices[[1]])))
                               
Informacion_1_DOWN <- data.frame(Mediana=rep(0,nrow(Matrices[[1]])),
                               Grupo=rep(0,nrow(Matrices[[1]])))                               

rownames(Informacion_1_UP) <- rownames(Matrices[[1]])
rownames(Informacion_1_DOWN) <- rownames(Matrices[[1]])

Informacion_2_UP <- data.frame(Mediana=rep(0,nrow(Matrices[[3]])),
                               Grupo=rep(0,nrow(Matrices[[3]])))

Informacion_2_DOWN <- data.frame(Mediana=rep(0,nrow(Matrices[[3]])),
                                 Grupo=rep(0,nrow(Matrices[[3]])))                               

rownames(Informacion_2_UP) <- rownames(Matrices[[3]])
rownames(Informacion_2_DOWN) <- rownames(Matrices[[3]])

Informacion_3_UP <- data.frame(Mediana=rep(0,nrow(Matrices[[5]])),
                               Grupo=rep(0,nrow(Matrices[[5]])))

Informacion_3_DOWN <- data.frame(Mediana=rep(0,nrow(Matrices[[5]])),
                                 Grupo=rep(0,nrow(Matrices[[5]])))                               

rownames(Informacion_3_UP) <- rownames(Matrices[[5]])
rownames(Informacion_3_DOWN) <- rownames(Matrices[[5]])


InformacionTotal <- list(Informacion_1_UP,Informacion_1_DOWN,Informacion_2_UP,Informacion_2_UP,Informacion_3_UP,Informacion_3_DOWN)
names(InformacionTotal) <- Condiciones

i <- 1

for (i in seq_along(Matrices)){
  Sumatorio <- SumatorioMatrices[[i]]
  Mediana <- MedianaMatrices[i]
  Informacion <- InformacionTotal[[i]]
  
  j <- 1
  
  for (j in seq_along(Sumatorio)){
    Valor <- Sumatorio[j]
    
    if ( Valor > Mediana ){
      Informacion[j,1] <- 1
      Informacion[j,2] <- "UP"}
    
    if ( Valor < Mediana ){
      Informacion[j,1] <- 0
      Informacion[j,2] <- "DOWN"}
  }
  
  InformacionTotal[[i]] <- Informacion
  
}

saveRDS(InformacionTotal,"InformacionMediana_UP_DOWN_Luad.rds")

#########################################################################################################

#······································ L U S C ···············································

# CALCULO DE LA MEDIANA: 

Matrices <- readRDS("Clustering/Grupos_UP_DOWN_Lusc.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 3 DOWN")
SumatorioMatrices <- list()
MedianaMatrices <- vector()
i <- 1

for (i in seq_along(Matrices)){
  Matriz <- Matrices[[i]]
  Sumatorio <- rowSums(Matriz)
  Mediana <- median(Sumatorio)
  SumatorioMatrices[[i]] <- Sumatorio
  MedianaMatrices[i] <- Mediana
}

names(MedianaMatrices) <- Condiciones
names(SumatorioMatrices) <- Condiciones

saveRDS(SumatorioMatrices,"SumatorioMatrices_UP_DOWN_Lusc.rds")
saveRDS(MedianaMatrices,"MedianaMatrices_UP_DOWN_Lusc.rds")

#···············································································································

# ESTIMACION DE SI SE ENCUENTRA POR ENCIMA O POR DEBAJO: 

Matrices <- readRDS("Clustering/Grupos_UP_DOWN_Lusc.rds")
SumatorioMatrices <- readRDS("Clustering/SumatorioMatrices_UP_DOWN_Lusc.rds")
MedianaMatrices <- readRDS("Clustering/MedianaMatrices_UP_DOWN_Lusc.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 3 DOWN")


Informacion_1_UP <- data.frame(Mediana=rep(0,nrow(Matrices[[1]])),
                               Grupo=rep(0,nrow(Matrices[[1]])))

Informacion_1_DOWN <- data.frame(Mediana=rep(0,nrow(Matrices[[1]])),
                                 Grupo=rep(0,nrow(Matrices[[1]])))                               

rownames(Informacion_1_UP) <- rownames(Matrices[[1]])
rownames(Informacion_1_DOWN) <- rownames(Matrices[[1]])

Informacion_2_UP <- data.frame(Mediana=rep(0,nrow(Matrices[[3]])),
                               Grupo=rep(0,nrow(Matrices[[3]])))

rownames(Informacion_2_UP) <- rownames(Matrices[[3]])

Informacion_3_DOWN <- data.frame(Mediana=rep(0,nrow(Matrices[[4]])),
                                 Grupo=rep(0,nrow(Matrices[[4]])))                               

rownames(Informacion_3_DOWN) <- rownames(Matrices[[4]])


InformacionTotal <- list(Informacion_1_UP,Informacion_1_DOWN,Informacion_2_UP,Informacion_3_DOWN)
names(InformacionTotal) <- Condiciones

i <- 1

for (i in seq_along(Matrices)){
  Sumatorio <- SumatorioMatrices[[i]]
  Mediana <- MedianaMatrices[i]
  Informacion <- InformacionTotal[[i]]
  
  j <- 1
  
  for (j in seq_along(Sumatorio)){
    Valor <- Sumatorio[j]
    
    if ( Valor > Mediana ){
      Informacion[j,1] <- 1
      Informacion[j,2] <- "UP"}
    
    if ( Valor < Mediana ){
      Informacion[j,1] <- 0
      Informacion[j,2] <- "DOWN"}
  }
  
  InformacionTotal[[i]] <- Informacion
  
}

saveRDS(InformacionTotal,"InformacionMediana_UP_DOWN_Lusc.rds")












