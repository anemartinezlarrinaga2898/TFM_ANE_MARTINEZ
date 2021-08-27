# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-05-2021

#OBJETIVO:Script ppara obtener las matrices con los pacientes de cada grupo. 

##############################################################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)

#########################################################################################################

#······································ L U A D ···············································

#Matrices de expresion de los genes diferencialmente expresados en todos los pacientes: 

Matriz_1 <- readRDS("GenesUpDown/Luad/Matrices_UP_DOWN_Luad_1.rds")
Matriz_2 <- readRDS("GenesUpDown/Luad/Matrices_UP_DOWN_Luad_2.rds")
Matriz_3 <- readRDS("GenesUpDown/Luad/Matrices_UP_DOWN_Luad_3.rds")

UP_1 <- Matriz_1[[1]]
DOWN_1 <- Matriz_1[[2]]

UP_2 <- Matriz_2[[1]]
DOWN_2 <- Matriz_2[[2]]

UP_3 <- Matriz_3[[1]]
DOWN_3 <- Matriz_3[[2]]

MatricesTotales <- list(UP_1,DOWN_1,UP_2,DOWN_2,UP_3,DOWN_3)
names(MatricesTotales) <- c("1_UP","1_DOWN","2_UP","2_DOWN","3_UP","3_DOWN")

saveRDS(MatricesTotales,"MatricesTotales_UP_DOWN_Luad.rds")

#Pertenencia a los grupos: 

K3.1 <- readRDS("PertenenciaGrupos/K3.1_NMF_Luad.rds")
K3.2 <- readRDS("PertenenciaGrupos/K3.2_NMF_Luad.rds")
K3.3 <- readRDS("PertenenciaGrupos/K3.3_NMF_Luad.rds")

Pacientes_3.1 <- K3.1$Sample_ID
Pacientes_3.2 <- K3.2$Sample_ID
Pacientes_3.3 <- K3.3$Sample_ID

saveRDS(Pacientes_3.1,"Pacientes_Grupo_3.1_Luad.rds")
saveRDS(Pacientes_3.2,"Pacientes_Grupo_3.2_Luad.rds")
saveRDS(Pacientes_3.3,"Pacientes_Grupo_3.3_Luad.rds")

PacientesTotales <- list(Pacientes_3.1,Pacientes_3.1,Pacientes_3.2,Pacientes_3.2,Pacientes_3.3,Pacientes_3.3)
names(PacientesTotales) <- c("Grupo 1","Grupo 1","Grupo 2","Grupo 2","Grupo 3","Grupo 3")

i <- 1
Genes_Pacientes_Totales <- list()

for (i in seq_along(MatricesTotales)){
  Matriz <- MatricesTotales[[i]] #Los genes estan en FILAS y los pacientes estan en las COLUMNAS
  Pacientes <- PacientesTotales[[i]]
  
  j <- 1
  Gen_Paciente <- matrix(0,nrow=length(Pacientes),ncol=nrow(Matriz))
  Gen_Paciente <- as.data.frame(Gen_Paciente)
  colnames(Gen_Paciente) <- rownames(Matriz) #Las columnas son los genes 
  rownames(Gen_Paciente) <- Pacientes # Las dilas son los pacientes
  for (j in seq_along(Pacientes)){
    Paciente <- Pacientes[j]
    Paciente_Gen_Expresion <- t(Matriz[,which(colnames(Matriz)==Paciente)])
    Gen_Paciente[j,] <- Paciente_Gen_Expresion
  }
  
  Genes_Pacientes_Totales[[i]] <- Gen_Paciente
}

names(Genes_Pacientes_Totales) <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 2 DOWN","Grupo 3 UP","Grupo 3 DOWN")

saveRDS(Genes_Pacientes_Totales,"Grupos_UP_DOWN_Luad.rds")
  
#########################################################################################################

#······································ L U S C ···············································  
  
#Matrices de expresion de los genes diferencialmente expresados en todos los pacientes: 

Matriz_1 <- readRDS("GenesUpDown/Lusc/Matrices_UP_DOWN_Lusc_1.rds")
Matriz_2 <- readRDS("GenesUpDown/Lusc/Matrices_UP_DOWN_Lusc_2.rds")
Matriz_3 <- readRDS("GenesUpDown/Lusc/Matrices_UP_DOWN_Lusc_3.rds")

UP_1 <- Matriz_1[[1]]
DOWN_1 <- Matriz_1[[2]]

UP_2 <- Matriz_2[[1]]

DOWN_3 <- Matriz_3[[1]]

MatricesTotales <- list(UP_1,DOWN_1,UP_2,DOWN_3)
names(MatricesTotales) <- c("1_UP","1_DOWN","2_UP","3_DOWN")

saveRDS(MatricesTotales,"MatricesTotales_UP_DOWN_Lusc.rds")

#Pertenencia a los grupos: 

K3.1 <- readRDS("PertenenciaGrupos/K3.1_NMF_Lusc.rds")
K3.2 <- readRDS("PertenenciaGrupos/K3.2_NMF_Lusc.rds")
K3.3 <- readRDS("PertenenciaGrupos/K3.3_NMF_Lusc.rds")

Pacientes_3.1 <- K3.1$Sample_ID
Pacientes_3.2 <- K3.2$Sample_ID
Pacientes_3.3 <- K3.3$Sample_ID

saveRDS(Pacientes_3.1,"Pacientes_Grupo_3.1_Lusc.rds")
saveRDS(Pacientes_3.2,"Pacientes_Grupo_3.2_Lusc.rds")
saveRDS(Pacientes_3.3,"Pacientes_Grupo_3.3_Lusc.rds")

PacientesTotales <- list(Pacientes_3.1,Pacientes_3.1,Pacientes_3.2,Pacientes_3.3)
names(PacientesTotales) <- c("Grupo 1","Grupo 1","Grupo 2","Grupo 3")

i <- 1
Genes_Pacientes_Totales <- list()

for (i in seq_along(MatricesTotales)){
  Matriz <- MatricesTotales[[i]] #Los genes estan en FILAS y los pacientes estan en las COLUMNAS
  Pacientes <- PacientesTotales[[i]]
  
  j <- 1
  Gen_Paciente <- matrix(0,nrow=length(Pacientes),ncol=nrow(Matriz))
  Gen_Paciente <- as.data.frame(Gen_Paciente)
  colnames(Gen_Paciente) <- rownames(Matriz) #Las columnas son los genes 
  rownames(Gen_Paciente) <- Pacientes # Las dilas son los pacientes
  for (j in seq_along(Pacientes)){
    Paciente <- Pacientes[j]
    Paciente_Gen_Expresion <- t(Matriz[,which(colnames(Matriz)==Paciente)])
    Gen_Paciente[j,] <- Paciente_Gen_Expresion
  }
  
  Genes_Pacientes_Totales[[i]] <- Gen_Paciente
}

names(Genes_Pacientes_Totales) <- c("Grupo 1 UP","Grupo 1 DOWN","Grupo 2 UP","Grupo 3 DOWN")

saveRDS(Genes_Pacientes_Totales,"Grupos_UP_DOWN_Lusc.rds")  
  
  