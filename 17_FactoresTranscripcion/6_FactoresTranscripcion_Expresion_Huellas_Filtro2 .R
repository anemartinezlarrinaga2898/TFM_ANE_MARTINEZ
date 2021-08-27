# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 01-07-2021

#OBJETIVO:Encontrar los niveles de expresion de los factores de transcripcion. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/17_FactoresTranscripcion/")

################################################################################################################################

#······························· L U A D ·····································

Huella1 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Luad_1.rds")
Huella2 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Luad_2.rds")
Huella3 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Luad_3.rds")

Huellas <- list(Huella1,Huella2,Huella3)

Expresion1 <- readRDS("GeneList_Comparacion_1_Luad_Total.rds")
Expresion2 <- readRDS("GeneList_Comparacion_2_Luad_Total.rds")
Expresion3 <- readRDS("GeneList_Comparacion_3_Luad_Total.rds")

NivelesExpresion <- list(Expresion1,Expresion2,Expresion3)

ListaFactoresTranscripcion <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
colnames(ListaFactoresTranscripcion) <- c("Ensembl_version","EntrezID","Symbol","ensembl_gene_id")

FactoresTranscripcion_Huellas_NivelesExpresion <- list()
HuellasNivelesExpresion <- list()

i <- 1

for(i in seq_along(Huellas)){
  
  Huella <- Huellas[[i]]
  
  Expresion <- NivelesExpresion[[i]]
  Expresion <- data.frame(Expresion)
  Expresion <- data.frame(Expresion,ensembl_gene_id=rownames(Expresion))
  
  X <- inner_join(Huella,Expresion,by="ensembl_gene_id")
  
  HuellasNivelesExpresion[[i]] <- X
  
  Y <- inner_join(X,ListaFactoresTranscripcion,by="ensembl_gene_id")
  Y <- Y[,c(3,4,6,7)]
  rownames(Y) <- Y[,1]
  Y <- Y[,-1]
  
  FactoresTranscripcion_Huellas_NivelesExpresion[[i]] <- Y
}

saveRDS(HuellasNivelesExpresion,"Huellas_Luad_Anotados_NivelesExpresion.rds")
saveRDS(FactoresTranscripcion_Huellas_NivelesExpresion,"FT_Luad_NE_Anotados.rds")
################################################################################################################################

#······························· L U S C ·····································

Huella1 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Lusc_1.rds")
Huella2 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Lusc_2.rds")
Huella3 <- readRDS("Huella_Ensembel_HUGO_EntrezGene_Lusc_3.rds")

Huellas <- list(Huella1,Huella2,Huella3)

Expresion1 <- readRDS("GeneList_Comparacion_1_Lusc_Total.rds")
Expresion2 <- readRDS("GeneList_Comparacion_2_Lusc_Total.rds")
Expresion3 <- readRDS("GeneList_Comparacion_3_Lusc_Total.rds")

NivelesExpresion <- list(Expresion1,Expresion2,Expresion3)

ListaFactoresTranscripcion <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
colnames(ListaFactoresTranscripcion) <- c("Ensembl_version","EntrezID","Symbol","ensembl_gene_id")

FactoresTranscripcion_Huellas_NivelesExpresion <- list()
HuellasNivelesExpresion <- list()

i <- 1

for(i in seq_along(Huellas)){
  
  Huella <- Huellas[[i]]
  
  Expresion <- NivelesExpresion[[i]]
  Expresion <- data.frame(Expresion)
  Expresion <- data.frame(Expresion,ensembl_gene_id=rownames(Expresion))
  
  X <- inner_join(Huella,Expresion,by="ensembl_gene_id")
  
  HuellasNivelesExpresion[[i]] <- X
  
  Y <- inner_join(X,ListaFactoresTranscripcion,by="ensembl_gene_id")
  Y <- Y[,c(3,4,6,7)]
  rownames(Y) <- Y[,1]
  Y <- Y[,-1]
  
  FactoresTranscripcion_Huellas_NivelesExpresion[[i]] <- Y
}

saveRDS(HuellasNivelesExpresion,"Huellas_Lusc_Anotados_NivelesExpresion.rds")
saveRDS(FactoresTranscripcion_Huellas_NivelesExpresion,"FT_Lusc_NE_Anotados.rds")
