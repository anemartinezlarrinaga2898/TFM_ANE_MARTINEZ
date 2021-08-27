# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 09-04-2021


#OBJETIVO: analizar el consensus clustering de Kmeans, y sacar la asignacion a cada grupo. 

####################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(tidyverse)
NumeroCluster <- seq(2,7,1)
###################################################################################

# ······························L U A D······································

Grupos_Kmeans_Luad <- readRDS("PertenciaGruposKmeans_Luad.rds")
ID_Pacientes_Luad <- unique(Grupos_Kmeans_Luad$item)


# Prueba <- Grupos_Kmeans_Luad[grep(ID_Pacientes_Luad[1],Grupos_Kmeans_Luad$item),]
# Prueba_1 <- Prueba %>% filter(k==NumeroCluster[1])
# Prueba_2 <- Prueba_1[which.max(Prueba_1$itemConsensus),]
# 
# Prueba3 <- Grupos_Kmeans_Luad[grep(ID_Pacientes_Luad[2],Grupos_Kmeans_Luad$item),]
# Prueba_4 <- Prueba3 %>% filter(k==NumeroCluster[1])
# Prueba_5 <- Prueba_4[which.max(Prueba_4$itemConsensus),]
# 
# Total <- rbind(Prueba_2,Prueba_5)

i <- 1
j <- 1
AsignacionCluster <- data.frame()

for (i in seq_along(NumeroCluster)){
  Grupo <- NumeroCluster[i] #Selecionamos el K
  AgrupacionCluster <- Grupos_Kmeans_Luad %>% filter(k==Grupo) #Asignamos grupo
  for (j in seq_along(ID_Pacientes_Luad)){
    Paciente <- ID_Pacientes_Luad[j] #Asignamos paciente
    InfoPaciente <- AgrupacionCluster[grep(Paciente,AgrupacionCluster$item),] #Lo buscamos
    InfoPacienteMaximo <- InfoPaciente[which.max(InfoPaciente$itemConsensus),] #Selecionamos el valor maximo
    AsignacionCluster_Kmeans <- rbind(AsignacionCluster,InfoPacienteMaximo) #Generamos el data.frame final
    AsignacionCluster <- AsignacionCluster_Kmeans
  }
}

#Comprobacion:

i <- 1
i <- i+1
for (i in seq_along(NumeroCluster)){
  Grupo <- NumeroCluster[i]
  Comprobacion <- AsignacionCluster_Kmeans %>% filter(k==Grupo)
  NumeroFilas_Comprobacion <- nrow(Comprobacion)
  
  if (NumeroFilas_Comprobacion==length(ID_Pacientes_Luad)){next}
  else{print(paste(Grupo,"no tiene 61 pacientes"))}
} 

TotalObservaciones <- 61*6 #61 pacientes por 6 clusters cantidad de filas que tendrian que salir en el data.frame final

saveRDS(AsignacionCluster_Kmeans,"AsignacionCluster_Kmeans_Luad.rds")

###################################################################################

# ······························L U S C······································

Grupos_Kmeans_Lusc <- readRDS("PertenciaGruposKmeans_Lusc.rds")
ID_Pacientes_Lusc <- unique(Grupos_Kmeans_Lusc$item)

i <- 1
j <- 1
AsignacionCluster <- data.frame()
AsignacionCluster_Kmeans<- data.frame()

for (i in seq_along(NumeroCluster)){
  Grupo <- NumeroCluster[i] #Selecionamos el K
  AgrupacionCluster <- Grupos_Kmeans_Lusc %>% filter(k==Grupo) #Asignamos grupo
  for (j in seq_along(ID_Pacientes_Lusc)){
    Paciente <- ID_Pacientes_Lusc[j] #Asignamos paciente
    InfoPaciente <- AgrupacionCluster[grep(Paciente,AgrupacionCluster$item),] #Lo buscamos
    InfoPacienteMaximo <- InfoPaciente[which.max(InfoPaciente$itemConsensus),] #Selecionamos el valor maximo
    AsignacionCluster_Kmeans <- rbind(AsignacionCluster,InfoPacienteMaximo) #Generamos el data.frame final
    AsignacionCluster <- AsignacionCluster_Kmeans
  }
}

#Comprobacion:

i <- 1
i <- i+1
for (i in seq_along(NumeroCluster)){
  Grupo <- NumeroCluster[i]
  Comprobacion <- AsignacionCluster_Kmeans %>% filter(k==Grupo)
  NumeroFilas_Comprobacion <- nrow(Comprobacion)
  
  if (NumeroFilas_Comprobacion==length(ID_Pacientes_Lusc)){next}
  else{print(paste(Grupo,"no tiene 495 pacientes"))}
} 

TotalObservaciones <- 495*6 #61 pacientes por 6 clusters cantidad de filas que tendrian que salir en el data.frame final

saveRDS(AsignacionCluster_Kmeans,"AsignacionCluster_Kmeans_Lusc.rds")
