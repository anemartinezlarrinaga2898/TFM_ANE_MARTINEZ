# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:Script para obtener las matrices de cuentas de mis genes UP y de mis genes DOWN en LUSC

##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)
#########################################################################################################
#······································ L U A D ···············································
Lusc <- readRDS("DatosPartida/Lusc_Filtrado.rds")
rownames(Lusc) <- gsub("\\..*","",rownames(Lusc))

Lusc_1_UP <- readRDS("GenesUpDown/Lusc_1_UP.rds")
Lusc_1_DOWN <- readRDS("GenesUpDown/Lusc_1_DOWN.rds")
Genes_1_UP <- rownames(Lusc_1_UP)
Genes_1_DOWN <- rownames(Lusc_1_DOWN)

Lusc_2_UP <- readRDS("GenesUpDown/Lusc_2_UP.rds")
Genes_2_UP <- rownames(Lusc_2_UP)

Lusc_3_DOWN <- readRDS("GenesUpDown/Lusc_3_DOWN.rds")
Genes_3_DOWN <- rownames(Lusc_3_DOWN)

GenesTotales <- list(Genes_1_UP,Genes_1_DOWN,Genes_2_UP,Genes_3_DOWN)

# 1) Luad_1_UP
# 2) Luad_1_DOWN
# 3) Luad_2_UP
# 4) Luad_3_DOWN

#En las matrices primero guardo el UP y luego guardo el DOWN
Matrices_1 <- list()
Matrices_2 <- list()
Matrices_3 <- list()

i <- 1
for (i in seq_along(GenesTotales)){
  Genes <- GenesTotales[[i]]
  j <- 1
  MatrizTranscripcional <- data.frame()
  
  for (j in seq_along(Genes)){
    Gen <- Genes[j]
    Gen_Paciente <- Lusc[which(rownames(Lusc)==Gen),]
    MatrizTranscripcional <- rbind(MatrizTranscripcional,Gen_Paciente)
  }
  
  if(i==1){Matrices_1[[1]] <- MatrizTranscripcional}
  if(i==2){Matrices_1[[2]] <- MatrizTranscripcional}
  if(i==3){Matrices_2[[1]] <- MatrizTranscripcional}
  if(i==4){Matrices_3[[1]] <- MatrizTranscripcional}
}


saveRDS(Matrices_1,"Matrices_UP_DOWN_Lusc_1.rds")
saveRDS(Matrices_2,"Matrices_UP_DOWN_Lusc_2.rds")
saveRDS(Matrices_3,"Matrices_UP_DOWN_Lusc_3.rds")
