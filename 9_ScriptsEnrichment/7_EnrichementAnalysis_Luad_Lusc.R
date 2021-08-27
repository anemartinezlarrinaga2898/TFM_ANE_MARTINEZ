# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 06-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS de mis TERMINOS. 

#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

#########################################################################################################

# ···································· L U A D················································

# C O M P A R A C I O N  1: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion1_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion1.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_1 <- enrichGO(gene=GenesFiltrados,
                        universe = GenesTotales,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE,
                        pAdjustMethod = "BH",
                        ont="BP")
saveRDS(ENRICHEMENT_1,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion1_Luad.rds")
ENRICHMENT_1 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion1_Luad.rds")
Comparacion1_Luad<- barplot(ENRICHMENT_1,title="Enrichment Luad comparacion 1")
Comparacion1_Luad

# C O M P A R A C I O N  2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion2.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_2 <- enrichGO(gene=GenesFiltrados,
                        universe = GenesTotales,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE,
                        pAdjustMethod = "BH",
                        ont="BP")

saveRDS(ENRICHEMENT_2,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion2_Luad.rds")
ENRICHEMENT_2 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion2_Luad.rds")
Comparacion2_Luad<- barplot(ENRICHEMENT_2,title="Enrichment Luad comparacion 2")
Comparacion2_Luad

# C O M P A R A C I O N  3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion3.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_3 <- enrichGO(gene=GenesFiltrados,
                        universe = GenesTotales,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE,
                        pAdjustMethod = "BH",
                        ont="BP")
saveRDS(ENRICHEMENT_3,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion3_Luad.rds")
ENRICHEMENT_3 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion3_Luad.rds")
Comparacion3_Luad<- barplot(ENRICHEMENT_3,title="Enrichment Luad comparacion 3")
Comparacion3_Luad

# ···································· L U S C················································

# C O M P A R A C I O N  1: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion1_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion1.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_1 <- enrichGO(gene=GenesFiltrados,
                          universe = GenesTotales,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE,
                          pAdjustMethod = "BH",
                          ont="BP")

saveRDS(ENRICHEMENT_1,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion1_Lusc.rds")
ENRICHEMENT_1 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion1_Lusc.rds")
Comparacion1_Lusc<- barplot(ENRICHEMENT_1,title="Enrichment Lusc comparacion 1")
Comparacion1_Lusc

# C O M P A R A C I O N  2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion2.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_2 <- enrichGO(gene=GenesFiltrados,
                          universe = GenesTotales,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE,
                          pAdjustMethod = "BH",
                          ont="BP")
saveRDS(ENRICHEMENT_2,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion2_Lusc.rds")
ENRICHEMENT_2 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion2_Lusc.rds")
Comparacion2_Luad<- barplot(ENRICHEMENT_2,title="Enrichment Lusc comparacion 2")
Comparacion2_Luad

# C O M P A R A C I O N  3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion3.rds")

GenesTotales <- rownames(ResultadosTotales)
GenesFiltrados <- rownames(ResultadosFiltrados)

ENRICHEMENT_3 <- enrichGO(gene=GenesFiltrados,
                          universe = GenesTotales,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE,
                          pAdjustMethod = "BH",
                          ont="BP")
saveRDS(ENRICHEMENT_3,"Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion3_Lusc.rds")
ENRICHEMENT_3 <- readRDS("Analisis_ORA_BP/GenesFiltrados/Enrichment_Comparacion3_Lusc.rds")
Comparacion3_Luad<- barplot(ENRICHEMENT_3,title="Enrichment Lusc comparacion 3")
Comparacion3_Luad
