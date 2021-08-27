# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 01-07-2021

#OBJETIVO:Encontrar los niveles de expresion de los factores de transcripcion en todos los genes que se encuentran diferencialmente expresados. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/17_FactoresTranscripcion/")
################################################################################################################################

Expresion1 <- readRDS("GeneList_Comparacion_1_Luad_Total.rds")
Expresion2 <- readRDS("GeneList_Comparacion_2_Luad_Total.rds")
Expresion3 <- readRDS("GeneList_Comparacion_3_Luad_Total.rds")

NivelesExpresion <- list(Expresion1,Expresion2,Expresion3)

ListaFactoresTranscripcion <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
colnames(ListaFactoresTranscripcion) <- c("Ensembl_version","EntrezID","Symbol","ensembl_gene_id")

FactoresTranscripcion_Huellas_NivelesExpresion <- list()

i <- 1

for(i in seq_along(NivelesExpresion)){
  
  Expresion <- NivelesExpresion[[i]]
  Expresion <- data.frame(Expresion)
  Expresion <- data.frame(Expresion,ensembl_gene_id=rownames(Expresion))
  
  X <- inner_join(ListaFactoresTranscripcion,Expresion,by="ensembl_gene_id")
  X <- X[,c(4,2,3,5)]
  
  FactoresTranscripcion_Huellas_NivelesExpresion[[i]] <- X
 
}

saveRDS(FactoresTranscripcion_Huellas_NivelesExpresion,"FT_Luad_NE_Anotados_Totales.rds")

################################################################################################################################

Expresion1 <- readRDS("GeneList_Comparacion_1_Lusc_Total.rds")
Expresion2 <- readRDS("GeneList_Comparacion_2_Lusc_Total.rds")
Expresion3 <- readRDS("GeneList_Comparacion_3_Lusc_Total.rds")

NivelesExpresion <- list(Expresion1,Expresion2,Expresion3)

ListaFactoresTranscripcion <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
colnames(ListaFactoresTranscripcion) <- c("Ensembl_version","EntrezID","Symbol","ensembl_gene_id")

FactoresTranscripcion_Huellas_NivelesExpresion <- list()

i <- 1

for(i in seq_along(NivelesExpresion)){
  
  Expresion <- NivelesExpresion[[i]]
  Expresion <- data.frame(Expresion)
  Expresion <- data.frame(Expresion,ensembl_gene_id=rownames(Expresion))
  
  X <- inner_join(ListaFactoresTranscripcion,Expresion,by="ensembl_gene_id")
  X <- X[,c(4,2,3,5)]
  
  FactoresTranscripcion_Huellas_NivelesExpresion[[i]] <- X
  
}

saveRDS(FactoresTranscripcion_Huellas_NivelesExpresion,"FT_Lusc_NE_Anotados_Totales.rds")

