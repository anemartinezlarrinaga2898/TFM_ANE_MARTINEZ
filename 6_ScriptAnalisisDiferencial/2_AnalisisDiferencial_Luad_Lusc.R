# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-04-2021

#OBJETIVO: Analisis diferencial de LUAD y de LUSC

#################################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)


#################################################################################################

# ··················· LUAD ······································

# 1) COMPARACION 1: 

Matriz <- readRDS("Comparacion1_Luad.rds")
ColData <- readRDS("Comparacion1ColData_Luad.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

    # OBJETO DESE1 COMPARACION 1: 
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                                     colData= ColData,
                                     design = ~ Grupo)

DDS_1_Luad <- DESeq(dds)
DDS_1_Luad_Resultados <- results(DDS_1_Luad)
DDS_1_Luad_Resultados

saveRDS(DDS_1_Luad,"AnalisisExpresion_1_Luad.rds")

#2) COMPARACION 2: 

Matriz <- readRDS("Comparacion2_Luad.rds")
ColData <- readRDS("Comparacion2ColData_Luad.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

# OBJETO DESE1 COMPARACION 2: 
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS_2_Luad <- DESeq(dds)
DDS_2_Luad_Resultados <- results(DDS_2_Luad)
DDS_2_Luad_Resultados


saveRDS(DDS_2_Luad,"AnalisisExpresion_2_Luad.rds")

#3) COMPARACION 3: 

Matriz <- readRDS("Comparacion3_Luad.rds")
ColData <- readRDS("Comparacion3ColData_Luad.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

# OBJETO DESE1 COMPARACION 3: 
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS_3_Luad <- DESeq(dds)
DDS_3_Luad_Resultados <- results(DDS_3_Luad)
DDS_3_Luad_Resultados


saveRDS(DDS_3_Luad,"AnalisisExpresion_3_Luad.rds")

#################################################################################################

# ··················· LUSC ······································

# 1) COMPARACION 1: 

Matriz <- readRDS("Comparacion1_Lusc.rds")
ColData <- readRDS("Comparacion1ColData_Lusc.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

# OBJETO DESE1 COMPARACION 1: 
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS <- DESeq(dds)
Resultados <- results(DDS)
saveRDS(DDS,"AnalisisExpresion_1_Lusc.rds")

#2) COMPARACION 2: 

Matriz <- readRDS("Comparacion2_Lusc.rds")
ColData <- readRDS("Comparacion2ColData_Lusc.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

# OBJETO DESE1 COMPARACION 2: 
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS <- DESeq(dds)
Resultados <- results(DDS)

saveRDS(DDS,"AnalisisExpresion_2_Lusc.rds")

#3) COMPARACION 3: 

Matriz <- readRDS("Comparacion3_Lusc.rds")
ColData <- readRDS("Comparacion3ColData_Lusc.rds")

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

# OBJETO DESE1 COMPARACION 3: 

dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS <- DESeq(dds)
Resultados <- results(DDS)
saveRDS(DDS,"AnalisisExpresion_3_Lusc.rds")

