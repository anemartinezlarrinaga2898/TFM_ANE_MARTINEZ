# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:Script para obtener la informacion clinica
##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)
#########################################################################################################

#······································ L U A D ···············································

Informacion <- readRDS("GenesUpDown/Luad/Informacion_UP_DOWN_SinClinica_Luad.rds")

Clinica <- readRDS("DatosPartida/ClinicalData_PacientesLuad_cBio.rds")
Clinica <- as.data.frame(Clinica[,c(2,67,68)])
colnames(Clinica) <- gsub(" ","",colnames(Clinica))
Clinica$OverallSurvivalStatus <- gsub("\\:.*","",Clinica$OverallSurvivalStatus)

UP_1 <- Informacion[[1]]
PatientID <- rownames(UP_1)

Listas <- list()
i <- 1
for (i in seq_along(Informacion)){
  Matriz <- Informacion[[i]]
  Matriz_Paciente <- cbind(PatientID,Matriz)
  Listas[[i]] <- Matriz_Paciente
}
saveRDS(Listas,"Media_Mediana_Pacientes_Luad.rds")
Clinica_Total <- list()

for (i in seq_along(Listas)){
  Matriz <- Listas[[i]]
  Union <- full_join(Matriz,Clinica,by="PatientID")
  Clinica_Total[[i]] <- Union
}

saveRDS(Clinica_Total,"InformacionPacientes_UP_DOWN_Luad.rds")

#······································ L U S C ···············································

Informacion <- readRDS("GenesUpDown/Lusc/Informacion_UP_DOWN_SinClinica_Lusc.rds")

Clinica <- readRDS("DatosPartida/ClinicalData_PacientesLusc_cBio.rds")
Clinica <- as.data.frame(Clinica[,c(2,67,68)])
colnames(Clinica) <- gsub(" ","",colnames(Clinica))
Clinica$OverallSurvivalStatus <- gsub("\\:.*","",Clinica$OverallSurvivalStatus)

UP_1 <- Informacion[[1]]
PatientID <- rownames(UP_1)

Listas <- list()
i <- 1
for (i in seq_along(Informacion)){
  Matriz <- Informacion[[i]]
  Matriz_Paciente <- cbind(PatientID,Matriz)
  Listas[[i]] <- Matriz_Paciente
}
saveRDS(Listas,"Media_Mediana_Pacientes_Lusc.rds")
Clinica_Total <- list()

i <- 1
for (i in seq_along(Listas)){
  Matriz <- Listas[[i]]
  Union <- full_join(Matriz,Clinica,by="PatientID")
  Clinica_Total[[i]] <- Union
}

saveRDS(Clinica_Total,"InformacionPacientes_UP_DOWN_Lusc.rds")
