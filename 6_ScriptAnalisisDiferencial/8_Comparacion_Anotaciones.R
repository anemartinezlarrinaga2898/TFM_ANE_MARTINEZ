# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-04-2021

#OBJETIVO: Comprobar cual de las dos bases de datos es capaz de anotar mas genes. 

#########################################################################################################
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")

# ·········································· L U A D ················································

DDS_1 <- readRDS("AnalisisExpresion_1_Luad.rds")
Resultados_1 <- results(DDS_1)
DDS_2 <- readRDS("AnalisisExpresion_2_Luad.rds")
Resultados_2 <- results(DDS_2)
DDS_3 <- readRDS("AnalisisExpresion_3_Luad.rds")
Resultados_3 <- results(DDS_3)

# ANOTACIONES CON ORG.HS.EG.DB

ResultadosORG.HS.DB_1 <- readRDS("ResultadosComparacion1_Luad_DF.rds")
ResultadosORG.HS.DB_2 <- readRDS("ResultadosComparacion2_Luad_DF.rds")
ResultadosORG.HS.DB_3 <- readRDS("ResultadosComparacion3_Luad_DF.rds")

Nan_1 <- ResultadosORG.HS.DB_1[which(is.na(ResultadosORG.HS.DB_1$SYMBOL)==TRUE),]
Nan_2 <- ResultadosORG.HS.DB_2[which(is.na(ResultadosORG.HS.DB_2$SYMBOL)==TRUE),]
Nan_3 <- ResultadosORG.HS.DB_3[which(is.na(ResultadosORG.HS.DB_3$SYMBOL)==TRUE),]

Nan_1 <- ResultadosORG.HS.DB_1[which(is.na(ResultadosORG.HS.DB_1$ENTREZID)==TRUE),]
Nan_2 <- ResultadosORG.HS.DB_2[which(is.na(ResultadosORG.HS.DB_2$ENTREZID)==TRUE),]
Nan_3 <- ResultadosORG.HS.DB_3[which(is.na(ResultadosORG.HS.DB_3$ENTREZID)==TRUE),]

#ANOTACIONES CON ENSDB.HSAPIENS.V75

Resultados_1$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(Resultados_1), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

ResultadosEnsDb_1 <- as.data.frame(mutate(as.data.frame(Resultados_1)),
                                     row.names=rownames(Resultados_1))

Nan_1 <- ResultadosEnsDb_1[which(is.na(ResultadosEnsDb_1$SYMBOL)==TRUE),]

Resultados_2$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados_2), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

ResultadosEnsDb_2 <- as.data.frame(mutate(as.data.frame(Resultados_2)),
                                   row.names=rownames(Resultados_2))

Nan_2 <- ResultadosEnsDb_2[which(is.na(ResultadosEnsDb_2$SYMBOL)==TRUE),]

Resultados_3$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados_3), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

ResultadosEnsDb_3 <- as.data.frame(mutate(as.data.frame(Resultados_3)),
                                   row.names=rownames(Resultados_3))

Nan_3 <- ResultadosEnsDb_3[which(is.na(ResultadosEnsDb_3$SYMBOL)==TRUE),]


# ·········································· L U S C ················································

ResultadosORG.HS.DB_1 <- readRDS("ResultadosComparacion1_Lusc_DF.rds")
ResultadosORG.HS.DB_2 <- readRDS("ResultadosComparacion2_Lusc_DF.rds")
ResultadosORG.HS.DB_3 <- readRDS("ResultadosComparacion3_Lusc_DF.rds")

Nan_1 <- ResultadosORG.HS.DB_1[which(is.na(ResultadosORG.HS.DB_1$SYMBOL)==TRUE),]
Nan_2 <- ResultadosORG.HS.DB_2[which(is.na(ResultadosORG.HS.DB_2$SYMBOL)==TRUE),]
Nan_3 <- ResultadosORG.HS.DB_3[which(is.na(ResultadosORG.HS.DB_3$SYMBOL)==TRUE),]

#ANOTACIONES CON ENSDB.HSAPIENS.V75

DDS_1 <- readRDS("AnalisisExpresion_1_Lusc.rds")
Resultados_1 <- results(DDS_1)
DDS_2 <- readRDS("AnalisisExpresion_2_Lusc.rds")
Resultados_2 <- results(DDS_2)
DDS_3 <- readRDS("AnalisisExpresion_3_Lusc.rds")
Resultados_3 <- results(DDS_3)

Resultados_1$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados_1), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

ResultadosEnsDb_1 <- as.data.frame(mutate(as.data.frame(Resultados_1)),
                                   row.names=rownames(Resultados_1))

Nan_1 <- ResultadosEnsDb_1[which(is.na(ResultadosEnsDb_1$SYMBOL)==TRUE),]

Resultados_2$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados_2), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

ResultadosEnsDb_2 <- as.data.frame(mutate(as.data.frame(Resultados_2)),
                                   row.names=rownames(Resultados_2))

Nan_2 <- ResultadosEnsDb_2[which(is.na(ResultadosEnsDb_2$SYMBOL)==TRUE),]

Resultados_3$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                               keys=rownames(Resultados_3), 
                               column="SYMBOL", 
                               keytype="GENEID", 
                               multiVals="first")

ResultadosEnsDb_3 <- as.data.frame(mutate(as.data.frame(Resultados_3)),
                                   row.names=rownames(Resultados_3))

Nan_3 <- ResultadosEnsDb_3[which(is.na(ResultadosEnsDb_3$SYMBOL)==TRUE),]
