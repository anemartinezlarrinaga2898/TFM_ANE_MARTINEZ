# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-07-2021

#OBJETIVO:Encontrar los genes del paper y la direccionalidad de los mismos. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/18_BusquedaGenes/")
library(pheatmap)
library(RColorBrewer)
################################################################################################################################

GenesPaper <- readRDS("GenesPaper_Ampliados.rds")
GenesPaper <- GenesPaper[!duplicated(GenesPaper$SYMBOL),]
GenesAAñadir <- data.frame(SYMBOL=c("FOXE1","NKX2-1","MUC1","SERPINB3","MMP13"),
                           TipoCelular=c("Clasico","Secretorio","Secretorio","Basal","Basal"))

GenesPaper <- rbind(GenesPaper,GenesAAñadir)
GenesPaper <- GenesPaper[order(GenesPaper$TipoCelular),]

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

GenesCaracteristicos <- rbind(Clasico,Secretorio)
GenesCaracteristicos <- rbind(GenesCaracteristicos,Basal)
GenesCaracteristicos <- rbind(GenesCaracteristicos,Primitivo)


#·······························································································································································

H1 <- readRDS("ResultadosComparacion1_Lusc_DF.rds")

H1$ensembl <- rownames(H1)
H1 <- H1[!duplicated(H1$SYMBOL),]
H1 <- na.omit(H1)
rownames(H1) <- H1$SYMBOL
H1$Grupo <- "H1"
H1 <- H1[,c(2,6,8,9)]

H2 <- readRDS("ResultadosComparacion2_Lusc_DF.rds")

H2$ensembl <- rownames(H2)
H2 <- H2[!duplicated(H2$SYMBOL),]
H2 <- na.omit(H2)
rownames(H2) <- H2$SYMBOL
H2$Grupo <- "H2"
H2 <- H2[,c(2,6,8,9)]

H3 <- readRDS("ResultadosComparacion3_Lusc_DF.rds")

H3$ensembl <- rownames(H3)
H3 <- H3[!duplicated(H3$SYMBOL),]
H3 <- na.omit(H3)
rownames(H3) <- H3$SYMBOL
H3$Grupo <- "H3"
H3 <- H3[,c(2,6,8,9)]

H1_P <- na.omit(H1[GenesPaper$SYMBOL,])
H2_P <- na.omit(H2[GenesPaper$SYMBOL,])
H3_P <- na.omit(H3[GenesPaper$SYMBOL,])

HT <- list(H1,H2,H3)
#·······························································································································································

Grupo <- c("Cluster1","Cluster2","Cluster3")
DataFrameFinal <- data.frame(matrix(0,ncol=length(Grupo),nrow=nrow(GenesCaracteristicos)))
colnames(DataFrameFinal) <- Grupo
rownames(DataFrameFinal) <- rownames(GenesCaracteristicos)

i <- 1

for(i in seq_along(HT)){
  
  Huella <- HT[[i]]
  
  j <- 1
  
  for(j in seq(1,nrow(DataFrameFinal),1)){
    
    GenInteres <- rownames(DataFrameFinal)[j]
    
    HuellaInteres <- Huella[which(rownames(Huella)==GenInteres),]
    
    if(i==1){DataFrameFinal[j,1] <- HuellaInteres[1,1]}
    if(i==2){DataFrameFinal[j,2] <- HuellaInteres[1,1]}
    if(i==3){DataFrameFinal[j,3] <- HuellaInteres[1,1]}
    
  }
  
}

DataFrameFinal[is.na(DataFrameFinal)] <-  0

Matriz <- matrix(DataFrameFinal,ncol=ncol(DataFrameFinal),nrow=nrow(DataFrameFinal))
Matriz<- as.numeric(unlist(Matriz))
Matriz <- matrix(Matriz,ncol=ncol(DataFrameFinal),nrow=nrow(DataFrameFinal))
colnames(Matriz) <- Grupo
rownames(Matriz) <- rownames(GenesCaracteristicos)


AnotacionesCelulas <- GenesCaracteristicos
AnotacionesCelulas <- data.frame(CelularType=AnotacionesCelulas$TipoCelular)
rownames(AnotacionesCelulas) <- GenesCaracteristicos$SYMBOL

AnotacionesGrupo <- data.frame(Cluster=c("Cluster1","Cluster2","Cluster3"))
rownames(AnotacionesGrupo) <- Grupo

Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(Matriz, 
         color = Pallete1,
         cluster_rows =F , 
         breaks = breaks,
         cluster_cols = T,
         show_colnames = T,
         show_rownames =T,
         annotation_row = AnotacionesCelulas,
         annotation_col = AnotacionesGrupo)

