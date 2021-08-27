# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 17-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
# Generacion de graficos para el analisis de GSEA

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)

library(DOSE)
##############################################################################################################

# ···································· L U A D················································

#Trabajaremos primero con el BP de los genes filtrados y totales: 

# C o m p a r a c i o n  1 : 

Resultados_1_Filtrados <- readRDS("AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_1_Filtrados.rds")
Resultados_1_Total <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_1_Total.rds")

TablaResultados_GSEA_1_Filtrado <- data.frame(Resultados_1_Filtrados@result)

GraficoResultados_1_F <- ggplot(TablaResultados_GSEA_1_Filtrado, aes(x=NES, y=Gene_set))+
                         geom_point(aes(colour=PAdjust),size=3.5)+
                         scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))
GraficoResultados_1_F

DotPlot_1_F <- dotplot(Resultados_1_Filtrados,showCategory=26)
DotPlot_1_F

TablaResultados_GSEA_1_Total <- data.frame(Gene_set=Resultados_1_Total@result$Description,
                                              NES=Resultados_1_Total@result$NES,
                                              PAdjust=Resultados_1_Total@result$p.adjust)

GraficoResultados_1_T <- ggplot(TablaResultados_GSEA_1_Total, aes(x=NES, y=Gene_set))+
  geom_point(aes(colour=PAdjust),size=3.5)+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))
GraficoResultados_1_T

DotPlot_1_T <- dotplot(Resultados_1_Total,showCategory=10)
DotPlot_1_T

prube <- Resultados_1_Filtrados@result

enrichplot::gseaplot2(Resultados_1_Filtrados,geneSetID = "GO:0009966")


