# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-05-2021

#OBJETIVO:Script para obtener que genes estan UP y que genes estan DOWN. 

##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)

#########################################################################################################

#······································ L U A D ···············································

Clinica <- readRDS("DatosPartida/ClinicalData_PacientesLuad_cBio.rds")
Clinica <- as.data.frame(Clinica[,c(2,67,68)])
colnames(Clinica) <- gsub(" ","",colnames(Clinica))
Clinica$OverallSurvivalStatus <- gsub("\\:.*","",Clinica$OverallSurvivalStatus)

Informacion <- readRDS("Clustering/InformacionMediana_UP_DOWN_Luad.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 2 DOWN","Grupo 3 UP","Grupo 3 DOWN")

Pacientes <- list()

i <- 1
for (i in seq_along(Informacion)){
  DataFrame <- Informacion[[i]]
  Pacientes[[i]] <- rownames(DataFrame)
}

names(Pacientes) <- Condiciones

i <- 1
ListaMedianaPacientes <- list()

for(i in seq_along(Pacientes)){
  PatientID <- Pacientes[[i]]
  Matriz <- Informacion[[i]]
  ListaMedianaPacientes[[i]] <- cbind(PatientID,Matriz)
}

saveRDS(ListaMedianaPacientes,"ListaMediana_Pacientes_UP_DOWN_Luad.rds")

i <- 1
Clinica_Total <- list()

for (i in seq_along(ListaMedianaPacientes)){
  Matriz <- ListaMedianaPacientes[[i]]
  Union <- inner_join(Matriz,Clinica,by="PatientID")
  Clinica_Total[[i]] <- Union
}

saveRDS(Clinica_Total,"Clinica_Total_Luad.rds")


#########################################################################################################

#······································ L U S C ···············································

Clinica <- readRDS("DatosPartida/ClinicalData_PacientesLusc_cBio.rds")
Clinica <- as.data.frame(Clinica[,c(2,67,68)])
colnames(Clinica) <- gsub(" ","",colnames(Clinica))
Clinica$OverallSurvivalStatus <- gsub("\\:.*","",Clinica$OverallSurvivalStatus)

Informacion <- readRDS("Clustering/InformacionMediana_UP_DOWN_Lusc.rds")
Condiciones <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 3 DOWN")

Pacientes <- list()

i <- 1
for (i in seq_along(Informacion)){
  DataFrame <- Informacion[[i]]
  Pacientes[[i]] <- rownames(DataFrame)
}

names(Pacientes) <- Condiciones

i <- 1
ListaMedianaPacientes <- list()

for(i in seq_along(Pacientes)){
  PatientID <- Pacientes[[i]]
  Matriz <- Informacion[[i]]
  ListaMedianaPacientes[[i]] <- cbind(PatientID,Matriz)
}

saveRDS(ListaMedianaPacientes,"ListaMediana_Pacientes_UP_DOWN_Lusc.rds")

i <- 1
Clinica_Total <- list()

for (i in seq_along(ListaMedianaPacientes)){
  Matriz <- ListaMedianaPacientes[[i]]
  Union <- inner_join(Matriz,Clinica,by="PatientID")
  Clinica_Total[[i]] <- Union
}

names(Clinica_Total) <- Condiciones
saveRDS(Clinica_Total,"Clinica_Total_Lusc.rds")


