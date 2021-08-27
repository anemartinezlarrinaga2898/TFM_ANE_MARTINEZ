# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:Script para obtener que genes estan UP y que genes estan DOWN. 

##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)
#########################################################################################################

#······································ L U A D ···············································

Comparacion1 <- readRDS("DatosPartida/ResultadosFiltrados_Luad_Comparacion1.rds")
Luad_1_UP <- Comparacion1 %>% dplyr::filter(log2FoldChange>2)
Luad_1_DOWN <- Comparacion1 %>% dplyr::filter(log2FoldChange<2)

saveRDS(Luad_1_UP,"Luad_1_UP.rds")
saveRDS(Luad_1_DOWN,"Luad_1_DOWN.rds")

Comparacion2 <- readRDS("DatosPartida/ResultadosFiltrados_Luad_Comparacion2.rds")
Luad_2_UP <- Comparacion2 %>% dplyr::filter(log2FoldChange>2)
Luad_2_DOWN <- Comparacion2 %>% dplyr::filter(log2FoldChange<2)

saveRDS(Luad_2_UP,"Luad_2_UP.rds")
saveRDS(Luad_2_DOWN,"Luad_2_DOWN.rds")

Comparacion3 <- readRDS("DatosPartida/ResultadosFiltrados_Luad_Comparacion3.rds")
Luad_3_UP <- Comparacion3 %>% dplyr::filter(log2FoldChange>2)
Luad_3_DOWN <- Comparacion3 %>% dplyr::filter(log2FoldChange<2)

saveRDS(Luad_3_UP,"Luad_3_UP.rds")
saveRDS(Luad_3_DOWN,"Luad_3_DOWN.rds")

#······································ L U S C ···············································

Comparacion1 <- readRDS("DatosPartida/ResultadosFiltrados_Lusc_Comparacion1.rds")
Lusc_1_UP <- Comparacion1 %>% dplyr::filter(log2FoldChange>2)
Lusc_1_DOWN <- Comparacion1 %>% dplyr::filter(log2FoldChange<2)

saveRDS(Lusc_1_UP,"Lusc_1_UP.rds")
saveRDS(Lusc_1_DOWN,"Lusc_1_DOWN.rds")

Comparacion2 <- readRDS("DatosPartida/ResultadosFiltrados_Lusc_Comparacion2.rds")
Lusc_2_UP <- Comparacion2 %>% dplyr::filter(log2FoldChange>2)
saveRDS(Lusc_2_UP,"Lusc_2_UP.rds")

Comparacion3 <- readRDS("DatosPartida/ResultadosFiltrados_Lusc_Comparacion3.rds")
Lusc_3_DOWN <- Comparacion3 %>% dplyr::filter(log2FoldChange<2)
saveRDS(Lusc_3_DOWN,"Lusc_3_DOWN.rds")

