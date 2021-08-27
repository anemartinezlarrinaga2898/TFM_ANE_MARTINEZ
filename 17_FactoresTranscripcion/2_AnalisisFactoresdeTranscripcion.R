# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-07-2021

#OBJETIVO:Trabajo con los factores de transcripcion

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/17_FactoresTranscripcion/")
library(RColorBrewer)
library(pheatmap)
library(MatrixGenerics)
################################################################################################################################

# 1) GENERACION HEATMAP MATRIZ ORIGINAL: 

MatrizLuad <- readRDS("MatrizLuad_Anotada.rds")
MatrizLusc <- readRDS("MatrizLusc_Anotada.rds")

FT <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")

FT_Luad <- readRDS("FT_TOTALES_LUAD.rds")
FT_Lusc <- readRDS("FT_TOTALES_LUSC.rds")

FT_LUAD_SinDuplicados <- FT_Luad$SYMBOL[!duplicated(FT_Luad$SYMBOL)]
FT_LUSC_SinDuplicados <- FT_Lusc$SYMBOL[!duplicated(FT_Lusc$SYMBOL)]

G1_Luad <- readRDS("K3.1_NMF_Luad.rds")
G2_Luad <- readRDS("K3.2_NMF_Luad.rds")
G3_Luad <- readRDS("K3.3_NMF_Luad.rds")

G1_Lusc <- readRDS("K3.1_NMF_Lusc.rds")
G2_Lusc <- readRDS("K3.2_NMF_Lusc.rds")
G3_Lusc <- readRDS("K3.3_NMF_Lusc.rds")

#·············································································

MLuad_FT <- MatrizLuad[FT_LUAD_SinDuplicados,] #La matriz de LUAD con todos los factores de transcripcion, con las cuentas 
MLusc_FT <- MatrizLusc[FT_LUSC_SinDuplicados,]

# A.1) Organizacion de los pacientes de LUAD

MG1_FT_Luad <- MLuad_FT[,G1_Luad$Sample_ID]
MG2_FT_Luad <- MLuad_FT[,G2_Luad$Sample_ID]
MG3_FT_Luad <- MLuad_FT[,G3_Luad$Sample_ID]

MLuad_FT <- cbind(MG1_FT_Luad,MG2_FT_Luad)
MLuad_FT <- cbind(MLuad_FT,MG3_FT_Luad) #Matriz de los FT de LUAD
MLuad_FT <- na.omit(MLuad_FT)
# A.2) Organizacion de los pacientes de LUSC

MG1_FT_Lusc <- MLusc_FT[,G1_Lusc$Sample_ID]
MG2_FT_Lusc <- MLusc_FT[,G2_Lusc$Sample_ID]
MG3_FT_Lusc <- MLusc_FT[,G3_Lusc$Sample_ID]

MLusc_FT <- cbind(MG1_FT_Lusc,MG2_FT_Lusc)
MLusc_FT <- cbind(MLusc_FT,MG3_FT_Lusc) #Matriz de los FT de LUSC

# Los objetos MLusc_FT son los objetos FINALES 

# A.3) Vamos a generar las anotaciones para los heatmaps: 

AnotacionesPacientesLuad <- data.frame(Pacientes=c(G1_Luad$Sample_ID,G2_Luad$Sample_ID,G3_Luad$Sample_ID),
                                       Cluster=c(G1_Luad$nmf_subtypes,G2_Luad$nmf_subtypes,G3_Luad$nmf_subtypes))
AnotacionesPacientesLuad$Cluster <- paste("Group",AnotacionesPacientesLuad$Cluster)
rownames(AnotacionesPacientesLuad) <- AnotacionesPacientesLuad$Pacientes

AnotacionesPacientesLusc <- data.frame(Pacientes=c(G1_Lusc$Sample_ID,G2_Lusc$Sample_ID,G3_Lusc$Sample_ID),
                                       Cluster=c(G1_Lusc$nmf_subtypes,G2_Lusc$nmf_subtypes,G3_Lusc$nmf_subtypes))
AnotacionesPacientesLusc$Cluster <- paste("Group",AnotacionesPacientesLusc$Cluster)
rownames(AnotacionesPacientesLusc) <- AnotacionesPacientesLusc$Pacientes

# A.4) Ahora vamos a generar los heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

MLuad_FT_Matrix <- as.numeric(unlist(MLuad_FT))
MLuad_FT_Matrix<- matrix(MLuad_FT_Matrix,ncol=ncol(MLuad_FT),nrow=nrow(MLuad_FT))
rownames(MLuad_FT_Matrix) <- rownames(MLuad_FT)
colnames(MLuad_FT_Matrix) <- colnames(MLuad_FT)

Varianzas <- rowVars(MLuad_FT_Matrix)
names(Varianzas) <- rownames(MLuad_FT_Matrix)
hist(Varianzas)

VarianzasSinCeros <- Varianzas[which(Varianzas>0.1)]

MLuad_FT_Matrix_DF<- data.frame(MLuad_FT_Matrix)

MatrizSinCeros <- MLuad_FT_Matrix_DF[names(VarianzasSinCeros),]

pheatmap(t(MLuad_FT_Matrix), 
                         color = Pallete1, 
                         cluster_rows =F , 
                         breaks = breaks,
                         cluster_cols = T, 
                         show_colnames = T,
                         show_rownames =F)
