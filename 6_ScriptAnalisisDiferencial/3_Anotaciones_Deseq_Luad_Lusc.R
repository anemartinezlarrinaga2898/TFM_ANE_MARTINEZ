# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-04-2021

#OBJETIVO: Volver a anotar mis genes en las 6 comparaciones de los 2 grupos con la nueva base de datos. 

#########################################################################################################
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(tidyverse)

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")

#########################################################################################################

# ········································L U A D ··············································

# C O M P A R A C I O N: 

DDS <- readRDS("AnalisisExpresion_1_Luad.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_1_Luad.rds")


DDS <- readRDS("AnalisisExpresion_2_Luad.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_2_Luad.rds")

DDS <- readRDS("AnalisisExpresion_3_Luad.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_3_Luad.rds")

# ········································L U S C ··············································

# C O M P A R A C I O N: 

DDS <- readRDS("AnalisisExpresion_1_Lusc.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_1_Lusc.rds")


DDS <- readRDS("AnalisisExpresion_2_Lusc.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_2_Lusc.rds")

DDS <- readRDS("AnalisisExpresion_3_Lusc.rds")
Resultados <- results(DDS)

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

saveRDS(Resultados,"ResultadosAnalasisExpresion_3_Lusc.rds")
