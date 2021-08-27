# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: generar los VOLCANO PLOTS de las tres comparaciones de la expresion diferencial. 

###########################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)
library(clusterProfiler)
library(ggrepel)
library(pheatmap)
library(enrichplot)
library(ggnewscale)
library(pathview)
library(ggpubr)

###############################################################################################

#······························· L U A D ····································

# 1) C O M P A R A C I O N  1:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_1_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

  #Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)


ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

  #Diferentes tipo de Filtracion:  

# FILTRACION 1: ······························································

ResultadosVolcano_Filtracion1 <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<1,
                                                                      "Gene DownRegulated",
                                                                       "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))
#Contaje Filtracion 1: 

GenesUp <- ResultadosVolcano_Filtracion1[which(ResultadosVolcano_Filtracion1$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano_Filtracion1[which(ResultadosVolcano_Filtracion1$FoldChange=="Gene DownRegulated",),]

# FILTRACION 2: ······························································

ResultadosVolcano_Filtracion2 <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                                FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<2,
                                                                                  "Gene DownRegulated",
                                                                                  "Gene UpRegulated"),
                                                                row.names=rownames(ResultadosFiltrados))
#Contaje Filtracion 2: 

GenesUp <- ResultadosVolcano_Filtracion2[which(ResultadosVolcano_Filtracion2$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano_Filtracion2[which(ResultadosVolcano_Filtracion2$FoldChange=="Gene DownRegulated",),]


VolcanoPlot_Filtrado1 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#63CAA7", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 1: Group 1 vs Group 2- Group 3")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado

# 2) C O M P A R A C I O N  2:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_2_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

#Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)

ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

#Grafico Volcano Plot Comparacion 1: 

ResultadosVolcano <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<0,
                                                                      "Gene DownRegulated",
                                                                      "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))

GenesUp <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene DownRegulated",),]


VolcanoPlot_Filtrado2 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#63CAA7", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 2: Group 2 vs Group 1- Group 3")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado

# 3) C O M P A R A C I O N  3:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_3_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

#Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)

ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

#Grafico Volcano Plot Comparacion 1: 

ResultadosVolcano <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<0,
                                                                      "Gene DownRegulated",
                                                                      "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))

GenesUp <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene DownRegulated",),]


VolcanoPlot_Filtrado3 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#63CAA7", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 3: Group 3 vs Group 1- Group 2")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado

VolcanoPlotsLuad <- ggarrange(VolcanoPlot_Filtrado1,VolcanoPlot_Filtrado2,VolcanoPlot_Filtrado3,
                        labels = c("A","B","C"), ncol = 3, nrow = 1)
annotate_figure(VolcanoPlotsLuad,
                top = text_grob("VolcanoPlots Luad", color = "Black", size = 14))
############################################################################################

#······························· L U S C ····································

# 1) C O M P A R A C I O N  1:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_1_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

#Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)

ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

#Grafico Volcano Plot Comparacion 1: 

ResultadosVolcano <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<0,
                                                                      "Gene DownRegulated",
                                                                      "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))

GenesUp <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene DownRegulated",),]

VolcanoPlot_Filtrado1 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#D5AAFF", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 1: Group 1 vs Group 2- Group 3")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado1

# 2) C O M P A R A C I O N  2:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

#Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)

ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

#Grafico Volcano Plot Comparacion 1: 

ResultadosVolcano <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<0,
                                                                      "Gene DownRegulated",
                                                                      "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))

GenesUp <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene DownRegulated",),]

VolcanoPlot_Filtrado2 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#D5AAFF", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 2: Group 2 vs Group 1- Group 3")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado2

# 3) C O M P A R A C I O N  3:(F I L T R A D O) ····························

Resultados <- readRDS("ResultadosAnalasisExpresion_3_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))

#Filtrado por aquellos: 

ResultadosFiltrados <- ResultadosDataFrame %>%filter(padj<0.05)

ResultadosFiltrados <- ResultadosFiltrados[order(ResultadosFiltrados$padj), ]

#Grafico Volcano Plot Comparacion 1: 

ResultadosVolcano <- ResultadosFiltrados %>% mutate(sig=ifelse(padj<0.05, "FDR<0.05", "Not Sig"),
                                                    FoldChange=ifelse(ResultadosFiltrados$log2FoldChange<0,
                                                                      "Gene DownRegulated",
                                                                      "Gene UpRegulated"),
                                                    row.names=rownames(ResultadosFiltrados))

GenesUp <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene UpRegulated",),]
GenesDown <- ResultadosVolcano[which(ResultadosVolcano$FoldChange=="Gene DownRegulated",),]


VolcanoPlot_Filtrado3 <- ggplot2::ggplot(ResultadosVolcano, ggplot2::aes(log2FoldChange, -log10(padj)))+
  ggplot2::geom_point(ggplot2::aes(col = sig, shape=FoldChange),size=2,alpha=0.8)+
  ggplot2::scale_color_manual(values = c("#D5AAFF", "#FFC172"))+ 
  ggplot2::scale_shape_manual(values = c(15,23))+
  ggplot2::ggtitle("Comparation 3: Group 3 vs Group 1- Group 2")+ 
  ggplot2::theme(axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 10,color="Grey"),
                 axis.text.y = element_text(size = 10),
                 axis.title.y = element_text(size = 10,color="Grey"),
                 plot.title = element_text(size = 14, face = "bold", color = "Black"))

VolcanoPlot_Filtrado3

VolcanoPlotsLusc <- ggarrange(VolcanoPlot_Filtrado1,VolcanoPlot_Filtrado2,VolcanoPlot_Filtrado3,
                              labels = c("A","B","C"), ncol = 3, nrow = 1)
annotate_figure(VolcanoPlotsLusc,
                top = text_grob("VolcanoPlots Lusc", color = "Black", size = 14))


