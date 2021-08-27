# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 12-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
                    # Todo se hara con el BIOLOGICAL PROCESS
  # 1) Realizar el analisis de enrichment con los genes que son UP y con los genes que son DOWN 
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

UP <- ResultadosFiltrados %>% dplyr::filter(log2FoldChange>2)
DOWN<- ResultadosFiltrados %>% dplyr::filter(log2FoldChange< -2)

UP <- UP[order(UP$log2FoldChange,decreasing = TRUE),]
DOWN <- DOWN[order(DOWN$log2FoldChange,decreasing=TRUE),]

GenesUP <- rownames(UP)
GenesDOWN <- rownames(DOWN)
GenesTotales <- rownames(ResultadosTotales)

set.seed(123456789)
ENRICHEMENT_UP_LUAD_1 <- enrichGO(gene=GenesUP,
                          universe = GenesTotales,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE,
                          pAdjustMethod = "BH",
                          ont="BP")

 saveRDS(ENRICHEMENT_UP_LUAD_1,"Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_1_BP.rds")

ENRICHEMENT_UP_LUAD_1 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_1_BP.rds")
Comparacion_UP_LUAD_1<- barplot(ENRICHEMENT_UP_LUAD_1 ,title="Enrichment Analysis LUAD (UP)")
Comparacion_UP_LUAD_1

ENRICHEMENT_DOWN_LUAD_1 <- enrichGO(gene=GenesDOWN,
                                  universe = GenesTotales,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

saveRDS(ENRICHEMENT_DOWN_LUAD_1,"Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_1_BP.rds")

ENRICHEMENT_DOWN_LUAD_1 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_1_BP.rds")
Comparacion_DOWN_LUAD_1<- barplot(ENRICHEMENT_DOWN_LUAD_1 ,title="Enrichment Analysis LUAD (DOWN)")
Comparacion_DOWN_LUAD_1


# C O M P A R A C I O N  2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion2.rds")

UP <- ResultadosFiltrados %>% dplyr::filter(log2FoldChange>2)
DOWN<- ResultadosFiltrados %>% dplyr::filter(log2FoldChange< -2)

UP <- UP[order(UP$log2FoldChange,decreasing = TRUE),]
DOWN <- DOWN[order(DOWN$log2FoldChange,decreasing=TRUE),]

GenesUP <- rownames(UP)
GenesDOWN <- rownames(DOWN)
GenesTotales <- rownames(ResultadosTotales)

set.seed(123456789)
ENRICHEMENT_UP_LUAD_2 <- enrichGO(gene=GenesUP,
                                  universe = GenesTotales,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

saveRDS(ENRICHEMENT_UP_LUAD_2,"Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_2_BP.rds")

ENRICHEMENT_UP_LUAD_2 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_2_BP.rds")
Comparacion_UP_LUAD_2<- barplot(ENRICHEMENT_UP_LUAD_2 ,title="Enrichment Analysis LUAD (UP)")
Comparacion_UP_LUAD_2

ENRICHEMENT_DOWN_LUAD_2 <- enrichGO(gene=GenesDOWN,
                                    universe = GenesTotales,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE,
                                    pAdjustMethod = "BH",
                                    ont="BP")

saveRDS(ENRICHEMENT_DOWN_LUAD_2,"Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_2_BP.rds")

ENRICHEMENT_DOWN_LUAD_2 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_2_BP.rds")
Comparacion_DOWN_LUAD_2<- barplot(ENRICHEMENT_DOWN_LUAD_2,title="Enrichment Analysis LUAD (DOWN)")
Comparacion_DOWN_LUAD_2

# C O M P A R A C I O N  3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion3.rds")

UP <- ResultadosFiltrados %>% dplyr::filter(log2FoldChange>2)
DOWN<- ResultadosFiltrados %>% dplyr::filter(log2FoldChange< -2)

UP <- UP[order(UP$log2FoldChange,decreasing = TRUE),]
DOWN <- DOWN[order(DOWN$log2FoldChange,decreasing=TRUE),]

GenesUP <- rownames(UP)
GenesDOWN <- rownames(DOWN)
GenesTotales <- rownames(ResultadosTotales)

set.seed(123456789)
ENRICHEMENT_UP_LUAD_3 <- enrichGO(gene=GenesUP,
                                  universe = GenesTotales,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

saveRDS(ENRICHEMENT_UP_LUAD_3,"Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_3_BP.rds")

ENRICHEMENT_UP_LUAD_3 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Luad_3_BP.rds")
Comparacion_UP_LUAD_3<- barplot(ENRICHEMENT_UP_LUAD_3 ,title="Enrichment Analysis LUAD (UP)")
Comparacion_UP_LUAD_3

ENRICHEMENT_DOWN_LUAD_3 <- enrichGO(gene=GenesDOWN,
                                    universe = GenesTotales,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE,
                                    pAdjustMethod = "BH",
                                    ont="BP")

saveRDS(ENRICHEMENT_DOWN_LUAD_3,"Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_3_BP.rds")

ENRICHEMENT_DOWN_LUAD_3 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Luad_3_BP.rds")
Comparacion_DOWN_LUAD_3<- barplot(ENRICHEMENT_DOWN_LUAD_3,title="Enrichment Analysis LUAD (DOWN)")
Comparacion_DOWN_LUAD_3

#########################################################################################################

# ···································· L U S C················································

# C O M P A R A C I O N  1: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion1_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion1.rds")

UP <- ResultadosFiltrados %>% dplyr::filter(log2FoldChange>2)
DOWN<- ResultadosFiltrados %>% dplyr::filter(log2FoldChange< -2)

UP <- UP[order(UP$log2FoldChange,decreasing = TRUE),]
DOWN <- DOWN[order(DOWN$log2FoldChange,decreasing=TRUE),]

GenesUP <- rownames(UP)
GenesDOWN <- rownames(DOWN)
GenesTotales <- rownames(ResultadosTotales)

set.seed(123456789)
ENRICHEMENT_UP_LUAD_1 <- enrichGO(gene=GenesUP,
                                  universe = GenesTotales,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

saveRDS(ENRICHEMENT_UP_LUAD_1,"Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Lusc_1_BP.rds")

ENRICHEMENT_UP_LUAD_1 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Lusc_1_BP.rds")
Comparacion_UP_LUAD_1<- barplot(ENRICHEMENT_UP_LUAD_1 ,title="Enrichment Analysis LUSC (UP)")
Comparacion_UP_LUAD_1

ENRICHEMENT_DOWN_LUAD_1 <- enrichGO(gene=GenesDOWN,
                                    universe = GenesTotales,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE,
                                    pAdjustMethod = "BH",
                                    ont="BP")

saveRDS(ENRICHEMENT_DOWN_LUAD_1,"Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Lusc_1_BP.rds")

ENRICHEMENT_DOWN_LUAD_1 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Lusc_1_BP.rds")
Comparacion_DOWN_LUAD_1<- barplot(ENRICHEMENT_DOWN_LUAD_1 ,title="Enrichment Analysis LUSC (DOWN)")
Comparacion_DOWN_LUAD_1


# C O M P A R A C I O N  2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion2.rds")

UP <- ResultadosFiltrados %>% dplyr::filter(log2FoldChange>2)

UP <- UP[order(UP$log2FoldChange,decreasing = TRUE),]

GenesUP <- rownames(UP)
GenesTotales <- rownames(ResultadosTotales)

set.seed(123456789)
ENRICHEMENT_UP_LUAD_2 <- enrichGO(gene=GenesUP,
                                  universe = GenesTotales,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

saveRDS(ENRICHEMENT_UP_LUAD_2,"Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Lusc_2_BP.rds")

ENRICHEMENT_UP_LUAD_2 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_UP_Lusc_2_BP.rds")
Comparacion_UP_LUAD_2<- barplot(ENRICHEMENT_UP_LUAD_2 ,title="Enrichment Analysis LUSC (UP)")
Comparacion_UP_LUAD_2



# C O M P A R A C I O N  3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion3.rds")

DOWN<- ResultadosFiltrados %>% dplyr::filter(log2FoldChange< -2)

DOWN <- DOWN[order(DOWN$log2FoldChange,decreasing=TRUE),]

GenesDOWN <- rownames(DOWN)
GenesTotales <- rownames(ResultadosTotales)

ENRICHEMENT_DOWN_LUAD_3 <- enrichGO(gene=GenesDOWN,
                                    universe = GenesTotales,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE,
                                    pAdjustMethod = "BH",
                                    ont="BP")

saveRDS(ENRICHEMENT_DOWN_LUAD_3,"Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Lusc_3_BP.rds")

ENRICHEMENT_DOWN_LUAD_3 <- readRDS("Analisis_ORA_BP/GenesUpDown/Enrichment_DOWN_Lusc_3_BP.rds")
Comparacion_DOWN_LUAD_3<- barplot(ENRICHEMENT_DOWN_LUAD_3,title="Enrichment Analysis LUSC (DOWN)")
Comparacion_DOWN_LUAD_3


