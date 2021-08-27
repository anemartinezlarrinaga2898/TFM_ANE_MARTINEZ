# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 12-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
#Generacion de las geneList para realizar el enrichment

#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
#########################################################################################################

# ···································· L U A D················································

                                    # COMPARACION 1: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion1_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion1.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_1_Luad_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  T O T A L E S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)

GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_1_Luad_Total.rds")

# ············································································································

                                          # COMPARACION 2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion2.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_2_Luad_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  T O T A L E S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)
GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_2_Luad_Total.rds")

# ············································································································

                                        # COMPARACION 3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Luad_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Luad_Comparacion3.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_3_Luad_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  T O T A L E S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)

GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_3_Luad_Total.rds")

# ···········································································································

#########################################################################################################

# ···································· L U S C················································

# COMPARACION 1: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion1_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion1.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_1_Lusc_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  T O D O S  L O S  G E N E S  F I L T R A D O S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)

GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_1_Lusc_Total.rds")

# ···························································································································

# COMPARACION 2: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion2_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion2.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_2_Lusc_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  T O D O S  L O S  G E N E S  F I L T R A D O S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)

GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_2_Lusc_Total.rds")

# ···························································································································

# COMPARACION 3: 

ResultadosTotales <- readRDS("DiferentialAnalysisResults/ResultadosComparacion3_Lusc_DF.rds")
ResultadosFiltrados <- readRDS("DiferentialAnalysisResults/ResultadosFiltrados_Lusc_Comparacion3.rds")

# A N A L I S I S  G S E A  C O N  L O S  G E N E S  F I L T R A D O S: 

GenesFiltrados <- ResultadosFiltrados
ValoresLogFoldChange <- GenesFiltrados$log2FoldChange
NombresGenesFiltrados <- rownames(GenesFiltrados)

GeneList <- data.frame(GeneId=NombresGenesFiltrados,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_3_Lusc_Filtrado.rds")

# A N A L I S I S  G S E A  C O N  T O D O S  L O S  G E N E S  F I L T R A D O S: 

GenesTotales<- ResultadosTotales
GenesTotales <- GenesTotales %>% filter(padj<0.05)
ValoresLogFoldChange <- GenesTotales$log2FoldChange
NombresGenesTotales <- rownames(GenesTotales)

GeneList <- data.frame(GeneId=NombresGenesTotales,LogFoldChange=ValoresLogFoldChange)
geneList <-  GeneList[,2]
names(geneList) <-  as.character(GeneList[,1])
geneList <-  sort(geneList, decreasing = TRUE)

saveRDS(geneList,"GeneList/GeneList_Comparacion_3_Lusc_Total.rds")
