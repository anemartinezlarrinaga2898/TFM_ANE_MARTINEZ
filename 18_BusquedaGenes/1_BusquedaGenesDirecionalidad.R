# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-07-2021

#OBJETIVO:Encontrar los genes del paper y la direccionalidad de los mismos. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/18_BusquedaGenes/")
################################################################################################################################

GenesPaper <- readRDS("GenesPaper_Ampliados.rds")
GenesPaper <- GenesPaper[!duplicated(GenesPaper$SYMBOL),]
GenesAAñadir <- data.frame(SYMBOL=c("FOXE1","NKX2-1","MUC1","SERPINB3","MMP13"),
                           TipoCelular=c("Clasico","Secretorio","Secretorio","Basal","Basal"))

GenesPaper <- rbind(GenesPaper,GenesAAñadir)
GenesPaper <- GenesPaper[order(GenesPaper$TipoCelular),]

H1 <- readRDS("ResultadosComparacion1_Lusc_DF.rds")
H1 <- H1 %>% filter(padj < 0.05)
H1$ensembl <- rownames(H1)
H1 <- H1[!duplicated(H1$SYMBOL),]
H1 <- na.omit(H1)
rownames(H1) <- H1$SYMBOL

H2 <- readRDS("ResultadosComparacion2_Lusc_DF.rds")
H2 <- H2 %>% filter(padj<0.05)
H2$ensembl <- rownames(H2)
H2 <- H2[!duplicated(H2$SYMBOL),]
H2 <- na.omit(H2)
rownames(H2) <- H2$SYMBOL

H3 <- readRDS("ResultadosComparacion3_Lusc_DF.rds")
H3 <- H3 %>% filter(padj<0.05)
H3$ensembl <- rownames(H3)
H3 <- H3[!duplicated(H3$SYMBOL),]
H3 <- na.omit(H3)
rownames(H3) <- H3$SYMBOL

#Primero trabajamos con los genes que estan UPREGULATED: 

H1_P <- H1[GenesPaper$SYMBOL,]
H1_P <- na.omit(H1_P)
H1_P <- H1_P %>% filter(log2FoldChange >=0)

H2_P <- H2[GenesPaper$SYMBOL,]
H2_P <- na.omit(H2_P)
H2_P <- H2_P %>% filter(log2FoldChange >=0)

H3_P <- H3[GenesPaper$SYMBOL,]
H3_P <- na.omit(H3_P)
H3_P <- H3_P %>% filter(log2FoldChange >=0)

Clasico <- GenesPaper %>% filter(TipoCelular=="Clasico")
rownames(Clasico) <- Clasico$SYMBOL

Secretorio <- GenesPaper %>% filter(TipoCelular=="Secretorio")
rownames(Secretorio) <- Secretorio$SYMBOL

Basal <- GenesPaper %>% filter(TipoCelular=="Basal")
rownames(Basal) <- Basal$SYMBOL

Primitivo <- GenesPaper %>% filter(TipoCelular=="Primitivo")
rownames(Primitivo) <- Primitivo$SYMBOL

GenesClasico <- c("TP63","AKR1C3","FOXE1","TXN")
GenesSecretorio <- c("NKX2-1","MUC1","ARHGDIB","TNFRSF14")
GenesBasal <- c("LAMB3","S100A7","MMP13","SERPINB3")
GenesPrimitivo <- c("MCM10","E2F3","TYMS","POLA1")

Clasico <- Clasico[GenesClasico,]
Secretorio <- Secretorio[GenesSecretorio,]
Basal <- Basal[GenesBasal,]
Primitivo <- Primitivo[GenesPrimitivo,]

H1_P_C <- na.omit(H1_P[Clasico$SYMBOL,])
H2_P_C <- na.omit(H2_P[Clasico$SYMBOL,])
H3_P_C <- na.omit(H3_P[Clasico$SYMBOL,])

H1_P_S <- na.omit(H1_P[Secretorio$SYMBOL,])
H2_P_S <- na.omit(H2_P[Secretorio$SYMBOL,])
H3_P_S <- na.omit(H3_P[Secretorio$SYMBOL,])

H1_P_B <- na.omit(H1_P[Basal$SYMBOL,])
H2_P_B <- na.omit(H2_P[Basal$SYMBOL,])
H3_P_B <- na.omit(H3_P[Basal$SYMBOL,])

H1_P_P <- na.omit(H1_P[Primitivo$SYMBOL,])
H2_P_P <- na.omit(H2_P[Primitivo$SYMBOL,])
H3_P_P <- na.omit(H3_P[Primitivo$SYMBOL,])


# Busqueda en la huella 1 un perfil que me permita asegurar que estoy con el grupo SECRETORIO 

H1_U <- H1[GenesPaper$SYMBOL,]
H1_U <- na.omit(H1_U)
H1_U <- H1_U %>% filter(log2FoldChange <=0)

H1_U_C <- na.omit(H1_U[Clasico$SYMBOL,])
H1_U_B <- na.omit(H1_U[Basal$SYMBOL,])

#··················

H2_U <- H2[GenesPaper$SYMBOL,]
H2_U <- na.omit(H2_U)
H2_U <- H2_U %>% filter(log2FoldChange <=0)

H2_U_C <- na.omit(H2_U[Clasico$SYMBOL,])
H2_U_B <- na.omit(H2_U[Basal$SYMBOL,])

# busqueda de la huella 3: 

H3_U <- H3[GenesPaper$SYMBOL,]
H3_U <- na.omit(H3_U)
H3_U <- H3_U %>% filter(log2FoldChange <=0)

H3_U_S <- na.omit(H3_U[Secretorio$SYMBOL,])
H3_U_B <- na.omit(H3_U[Basal$SYMBOL,])

#Busqueda de la huella 2: 

H2_U <- H2[GenesPaper$SYMBOL,]
H2_U <- na.omit(H2_U)
H2_U <- H2_U %>% filter(log2FoldChange <=0)

H2_U_C <- na.omit(H2_U[Clasico$SYMBOL,])
H2_U_S <- na.omit(H2_U[Secretorio$SYMBOL,])



