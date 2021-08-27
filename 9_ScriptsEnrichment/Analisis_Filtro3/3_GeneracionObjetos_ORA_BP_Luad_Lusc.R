# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 01-06-2021

#OBJETIVO: Generacion de las DF con el filtro 3 para poder hacer el BP con el algoritmo ORA

################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)

################################################################################################################################
#··································  L U A D ··········································
#RC1: Resultados comparacion 1

RC1 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion1_Luad_DF.rds")
RC2 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion2_Luad_DF.rds")
RC3 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion3_Luad_DF.rds")

#RFC1: Resultados comparacion filtrados 1

RCF1 <- RC1 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
RCF2 <- RC2 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
RCF3 <- RC3 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

saveRDS(RCF1,"Analisis_Filtro3/ResultadosFiltrados/RFC1_Luad.rds")
saveRDS(RCF2,"Analisis_Filtro3/ResultadosFiltrados/RFC2_Luad.rds")
saveRDS(RCF3,"Analisis_Filtro3/ResultadosFiltrados/RFC3_Luad.rds")

# GT : Genes Totales
GT <- rownames(RC1)

# GF1: Genes filtrados 1
GF1 <- rownames(RCF1)
GF2 <- rownames(RCF2)
GF3 <- rownames(RCF3)

saveRDS(GT,"Analisis_Filtro3/GFT_Luad.rds")
saveRDS(GF1,"Analisis_Filtro3/GF1_Luad.rds")
saveRDS(GF2,"Analisis_Filtro3/GF2_Luad.rds")
saveRDS(GF3,"Analisis_Filtro3/GF3_Luad.rds")

# RCF1U: Resultados filtrados 1 UP 
# RCF1D: Resultados filtrados 1 DOWN 

RCF1U <-RCF1 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF1D <-RCF1 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

RCF2U <-RCF2 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF2D <-RCF2 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

RCF3U <-RCF3 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF3D <-RCF3 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

saveRDS(RCF1U,"Analisis_Filtro3/ResultadosFiltrados/RCF1U_Luad.rds")
saveRDS(RCF1D,"Analisis_Filtro3/ResultadosFiltrados/RCF1D_Luad.rds")

saveRDS(RCF2U,"Analisis_Filtro3/ResultadosFiltrados/RCF2U_Luad.rds")
saveRDS(RCF2D,"Analisis_Filtro3/ResultadosFiltrados/RCF2D_Luad.rds")

saveRDS(RCF3U,"Analisis_Filtro3/ResultadosFiltrados/RCF3U_Luad.rds")
saveRDS(RCF3D,"Analisis_Filtro3/ResultadosFiltrados/RCF3D_Luad.rds")

GF1U <- rownames(RCF1U)
GF1D <- rownames(RCF1D)

GF2U <- rownames(RCF2U)
GF2D <- rownames(RCF2D)

GF3U <- rownames(RCF3U)
GF3D <- rownames(RCF3D)

saveRDS(GF1U,"Analisis_Filtro3/ResultadosFiltrados/GF1U_Luad.rds")
saveRDS(GF1D,"Analisis_Filtro3/ResultadosFiltrados/GF1D_Luad.rds")

saveRDS(GF2U,"Analisis_Filtro3/ResultadosFiltrados/GF2U_Luad.rds")
saveRDS(GF2D,"Analisis_Filtro3/ResultadosFiltrados/GF2D_Luad.rds")

saveRDS(GF3U,"Analisis_Filtro3/ResultadosFiltrados/GF3U_Luad.rds")
saveRDS(GF3D,"Analisis_Filtro3/ResultadosFiltrados/GF3D_Luad.rds")

################################################################################################################################
#··································  L U S C ··········································
#RC1: Resultados comparacion 1

RC1 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion1_Lusc_DF.rds")
RC2 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion2_Lusc_DF.rds")
RC3 <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion3_Lusc_DF.rds")

#RFC1: Resultados comparacion filtrados 1

RCF1 <- RC1 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
RCF2 <- RC2 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
RCF3 <- RC3 %>% dplyr::filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

saveRDS(RCF1,"Analisis_Filtro3/ResultadosFiltrados/RFC1_Lusc.rds")
saveRDS(RCF2,"Analisis_Filtro3/ResultadosFiltrados/RFC2_Lusc.rds")
saveRDS(RCF3,"Analisis_Filtro3/ResultadosFiltrados/RFC3_Lusc.rds")

# GT : Genes Totales
GT <- rownames(RC1)

# GF1: Genes filtrados 1
GF1 <- rownames(RCF1)
GF2 <- rownames(RCF2)
GF3 <- rownames(RCF3)

saveRDS(GT,"Analisis_Filtro3/GFT_Lusc.rds")
saveRDS(GF1,"Analisis_Filtro3/GF1_Lusc.rds")
saveRDS(GF2,"Analisis_Filtro3/GF2_Lusc.rds")
saveRDS(GF3,"Analisis_Filtro3/GF3_Lusc.rds")

# RCF1U: Resultados filtrados 1 UP 
# RCF1D: Resultados filtrados 1 DOWN 

RCF1U <-RCF1 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF1D <-RCF1 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

RCF2U <-RCF2 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF2D <-RCF2 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

RCF3U <-RCF3 %>% dplyr::filter (padj<0.05 & (log2FoldChange >1.5))
RCF3D <-RCF3 %>% dplyr::filter (padj<0.05 & (log2FoldChange <1.5))

saveRDS(RCF1U,"Analisis_Filtro3/ResultadosFiltrados/RCF1U_Lusc.rds")
saveRDS(RCF1D,"Analisis_Filtro3/ResultadosFiltrados/RCF1D_Lusc.rds")

saveRDS(RCF2U,"Analisis_Filtro3/ResultadosFiltrados/RCF2U_Lusc.rds")
saveRDS(RCF2D,"Analisis_Filtro3/ResultadosFiltrados/RCF2D_Lusc.rds")

saveRDS(RCF3U,"Analisis_Filtro3/ResultadosFiltrados/RCF3U_Lusc.rds")
saveRDS(RCF3D,"Analisis_Filtro3/ResultadosFiltrados/RCF3D_Lusc.rds")

GF1U <- rownames(RCF1U)
GF1D <- rownames(RCF1D)

GF2U <- rownames(RCF2U)
GF2D <- rownames(RCF2D)

GF3U <- rownames(RCF3U)
GF3D <- rownames(RCF3D)

saveRDS(GF1U,"Analisis_Filtro3/ResultadosFiltrados/GF1U_Lusc.rds")
saveRDS(GF1D,"Analisis_Filtro3/ResultadosFiltrados/GF1D_Lusc.rds")

saveRDS(GF2U,"Analisis_Filtro3/ResultadosFiltrados/GF2U_Lusc.rds")
saveRDS(GF2D,"Analisis_Filtro3/ResultadosFiltrados/GF2D_Lusc.rds")

saveRDS(GF3U,"Analisis_Filtro3/ResultadosFiltrados/GF3U_Lusc.rds")
saveRDS(GF3D,"Analisis_Filtro3/ResultadosFiltrados/GF3D_Lusc.rds")
