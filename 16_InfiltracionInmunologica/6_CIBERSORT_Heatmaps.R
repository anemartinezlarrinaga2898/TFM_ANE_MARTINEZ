# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-07-2021

# OBJETIVO:Heatmaps de la infiltracion inmunologica con CIBERSORT

##########################################################################################################################################
directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(MatrixGenerics)
##########################################################################################################################################

# 1) Cargamos los datos: 

LUAD <- read.csv("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/CIBERSORT/ResultadosLuad/TOTALES_1000.csv")
rownames(LUAD) <- LUAD$Input.Sample
LUAD <- LUAD %>% filter(P.value<0.05) #Aquellos en los que la correlacion no es signifivativa. 
LUAD <- LUAD[,-c(1,24,25,26)]

# 2) Convertimos el csv en una matrix: 

LUAD_matriz <- matrix(LUAD,ncol=ncol(LUAD),nrow=nrow(LUAD))
LUAD_matriz <- as.numeric(unlist(LUAD_matriz))
LUAD_matriz<- matrix(LUAD_matriz,ncol=ncol(LUAD),nrow=nrow(LUAD))
rownames(LUAD_matriz) <- rownames(LUAD)
colnames(LUAD_matriz) <- colnames(LUAD)

# 3) Obtenemos las anotaciones para la grafica: 

K3 <- readRDS("K3_NMF_Luad.rds")
K3 <- K3[rownames(LUAD_matriz),]

AnotacionesLuad <- data.frame(Cluster=K3$nmf_subtypes)
rownames(AnotacionesLuad) <- rownames(K3)


# 4) Representamos tal cual, el csv sin hacer ningun tipo de cambio. 

Colors <- c("green","black","red")
breaks <- seq(0,0.2,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(LUAD_matriz, 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = F,
         show_colnames = T,
         show_rownames =F,
         annotation_row = AnotacionesLuad)

########################################################################################################################################################################

#··························· L U S C ·········································

# 1) Cargamos los datos: 

LUAD <- read.csv("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/CIBERSORT/ResultadosLusc/TOTAL.csv")
rownames(LUAD) <- LUAD$Input.Sample
LUAD <- LUAD %>% filter(P.value<0.05) #Aquellos en los que la correlacion no es signifivativa. 
LUAD <- LUAD[,-c(1,24,25,26)]

# 2) Convertimos el csv en una matrix: 

LUAD_matriz <- matrix(LUAD,ncol=ncol(LUAD),nrow=nrow(LUAD))
LUAD_matriz <- as.numeric(unlist(LUAD_matriz))
LUAD_matriz<- matrix(LUAD_matriz,ncol=ncol(LUAD),nrow=nrow(LUAD))
rownames(LUAD_matriz) <- rownames(LUAD)
colnames(LUAD_matriz) <- colnames(LUAD)

# 3) Obtenemos las anotaciones para la grafica: 

K3 <- readRDS("K3_NMF_Lusc.rds")
K3 <- K3[rownames(LUAD_matriz),]

AnotacionesLuad <- data.frame(Cluster=K3$nmf_subtypes)
rownames(AnotacionesLuad) <- rownames(K3)


# 4) Representamos tal cual, el csv sin hacer ningun tipo de cambio. 

Colors <- c("green","black","red")
breaks <- seq(0,0.2,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(LUAD_matriz, 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = F,
         show_colnames = T,
         show_rownames =F,
         annotation_row = AnotacionesLuad)