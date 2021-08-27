# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-06-2021

#OBJETIVO:Script para analizar los top 10 genes que tenemos en cada una de las huellas tanto UP regulado como DOWN regulado 

#################################################################################################################
library(tidyverse)

#######################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/12_Huellas/")

GenesGrupo1 <- readRDS("ResultadosFiltrados_Luad_Comparacion1.rds")
GenesGrupo2 <- readRDS("ResultadosFiltrados_Luad_Comparacion2.rds")
GenesGrupo3 <- readRDS("ResultadosFiltrados_Luad_Comparacion3.rds")

Genes <- list(GenesGrupo1,GenesGrupo2,GenesGrupo3)
names(Genes) <- c("Genes Grupo 1","Genes Grupo 2","Genes Grupo 3")

ListaGenesOrganizada_DOWN <- list()
ListaGenesOrganizada_UP <- list()
i <- 1

for( i in seq_along(Genes)){
  ListaGenes <- Genes[[i]]
  
  ListaGenesOrganizada_DOWN[[i]] <- ListaGenes[order(ListaGenes$log2FoldChange),]
  ListaGenesOrganizada_UP[[i]] <- ListaGenes[order(ListaGenes$log2FoldChange),]
}

names(ListaGenesOrganizada_DOWN) <- c("Genes Grupo 1","Genes Grupo 2","Genes Grupo 3")
names(ListaGenesOrganizada_UP) <- c("Genes Grupo 1","Genes Grupo 2","Genes Grupo 3")




