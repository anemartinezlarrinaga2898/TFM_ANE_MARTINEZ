# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-06-2021

#OBJETIVO:
   
    #1) Obtencion de la matriz

#######################################################################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(biomaRt)

ensembel <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#######################################################################################################################################

# Obtencion de la matriz de las lineas con todos los genes: 

LINES <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/Cell_lines_annotations_20181226.txt")
Datos_LineasCelulares <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/TMP.txt")
Datos <- Datos_LineasCelulares[,-2]
Datos$gene_id <- gsub("\\..*","",Datos$gene_id)
rownames(Datos) <- Datos$gene_id
colnames(Datos)[1] <- "GeneID"
Datos <- Datos[,-1]

AnotacionesLineas <- unique(LINES$PATHOLOGIST_ANNOTATION)
LineasLusc <- LINES[which(LINES$PATHOLOGIST_ANNOTATION=="Lung:NSCLC_Squamous"),]

Matriz <- Datos[,LineasLusc$CCLE_ID]

saveRDS(Matriz,"MatrizLineasCelulares_TotalGenes.rds")