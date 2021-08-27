# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 24-05-2021

#OBJETIVO:Script para obtener las descripciones io las matrices de mis huellas en los datos de las lineas celulares. 

#############################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)

Datos_LineasCelulares <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/TMP.txt")

#############################################################################################################################################

Huellas1 <- readRDS("Huella_Ensembel_Lusc_1.rds")
Huellas2 <- readRDS("Huella_Ensembel_Lusc_2.rds")
Huellas3 <- readRDS("Huella_Ensembel_Lusc_3.rds")

Huellas <- list(Huellas1,Huellas2,Huellas3)
names(Huellas) <- c("Huella Grupo 1","Huella Grupo 2","Huella Grupo 3")
saveRDS(Huellas,"Huellas_Totales_Ensembel_Lusc.rds")

Datos <- Datos_LineasCelulares[,-2]
Datos$gene_id <- gsub("\\..*","",Datos$gene_id)
rownames(Datos) <- Datos$gene_id
colnames(Datos)[1] <- "GeneID"
Datos <- Datos[,-1]

Huellas_LineasCelulares <- list()
GenesNoEncontrados_LineasCelulares <- list()
i <- 1

for ( i in seq_along(Huellas)){
  Huella <- Huellas[[i]]
  
  ObtencionHuella <- data.frame()
  GenesNoEncontrados <- rep(0,length(Huella))
  j <- 1
  
  for ( j in seq_along(Huella)){
    Gen <- Huella[j]
    
    X <- Datos[which(rownames(Datos)==Gen),]
    ObtencionHuella <- rbind(ObtencionHuella,X)
    
    NumeroFilas <- nrow(X)
    
    if (NumeroFilas==0){GenesNoEncontrados[j] <- Gen}
    
  }
  
  Huellas_LineasCelulares[[i]] <- ObtencionHuella
  GenesNoEncontrados_LineasCelulares[[i]] <- GenesNoEncontrados
}

remove <- "0"
i <- 1
for (i in seq_along(GenesNoEncontrados_LineasCelulares)){
  Genes <- GenesNoEncontrados_LineasCelulares[[i]]
  EliminacionCeros <- Genes [! Genes %in% remove]
  
  GenesNoEncontrados_LineasCelulares[[i]] <- EliminacionCeros
}

names(GenesNoEncontrados_LineasCelulares) <- c("Grupo 1","Grupo 2","Grupo 3")
names(Huellas_LineasCelulares) <- c("Grupo 1","Grupo 2","Grupo 3")

saveRDS(GenesNoEncontrados_LineasCelulares,"GenesNoEncontrados_LineasCelulares_Lusc.rds")
saveRDS(Huellas_LineasCelulares,"Huellas_LineasCelulares_Lusc.rds")



