# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 12-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
# Todo se hara con el BIOLOGICAL PROCESS
# 1) analisis gsea con todos los genes y los filtrados

#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
#########################################################################################################

# ···································· L U A D················································

# COMPARACION 1: 

geneList <- readRDS("GeneList/GeneList_Comparacion_1_Luad_Filtrado.rds")
gseaGO_F <- gseGO(geneList = geneList, 
                  OrgDb = org.Hs.eg.db, 
                  ont = 'BP',  
                  minGSSize = 10, 
                  keyType = "ENSEMBL",
                  pvalueCutoff = 0.05,
                  verbose = FALSE,
                  eps=0) 
saveRDS(gseaGO_F,"AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_1_Filtrados.rds")

geneList <- readRDS("GeneList/GeneList_Comparacion_1_Luad_Total.rds")
gseaGO_T <- gseGO(geneList = geneList, 
                  OrgDb = org.Hs.eg.db, 
                  ont = 'BP',  
                  minGSSize = 10, 
                  keyType = "ENSEMBL",
                  pvalueCutoff = 0.05,
                  verbose = FALSE,
                  eps = 0) 
saveRDS(gseaGO_T,"AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_1_Total.rds")

# COMPARACION 2: 

geneList <- readRDS("GeneList/GeneList_Comparacion_2_Luad_Filtrado.rds")

gseaGO_F <- gseGO(geneList = geneList, 
                  OrgDb = org.Hs.eg.db, 
                  ont = 'BP',  
                  minGSSize = 10, 
                  keyType = "ENSEMBL",
                  pvalueCutoff = 0.05,
                  verbose = FALSE,
                  eps=0) 

saveRDS(gseaGO_F,"AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_2_Filtrados.rds")

geneList <- readRDS("GeneList/GeneList_Comparacion_2_Luad_Total.rds")

gseaGO_T <- gseGO(geneList = geneList, 
                   OrgDb = org.Hs.eg.db, 
                   ont = 'BP',  
                   minGSSize = 10, 
                   keyType = "ENSEMBL",
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   eps = 0) 
saveRDS(gseaGO_T,"AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_2_Total.rds")

# COMPARACION 3: 

geneList <- readRDS("GeneList/GeneList_Comparacion_3_Luad_Filtrado.rds")
gseaGO_F <- gseGO(geneList = geneList, 
                  OrgDb = org.Hs.eg.db, 
                  ont = 'BP',  
                  minGSSize = 10, 
                  keyType = "ENSEMBL",
                  pvalueCutoff = 0.05,
                  verbose = FALSE,
                  eps=0) 
saveRDS(gseaGO_F,"AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_3_Filtrados.rds")


geneList <- readRDS("GeneList/GeneList_Comparacion_3_Luad_Total.rds")
gseaGO1_T <- gseGO(geneList = geneList, 
                   OrgDb = org.Hs.eg.db, 
                   ont = 'BP',  
                   minGSSize = 10, 
                   keyType = "ENSEMBL",
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   eps = 0) 
saveRDS(gseaGO1_T,"AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_3_Total.rds")


