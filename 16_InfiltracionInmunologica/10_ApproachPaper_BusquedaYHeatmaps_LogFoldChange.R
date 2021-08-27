# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 21-07-2021

# OBJETIVO:Infiltracion inmunologica con los valores del logfoldchange
##########################################################################################################################################
directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

###########################################################################################################################################

# 1) Cargamos los genes que forman parte de mi huella: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$GenesInmune),]
rownames(HI) <- HI$GenesInmune

# 2) Cargamos los resultados del analisis y los convertimos en DF 

G1 <- readRDS("ResultadosAnalasisExpresion_1_Luad.rds")
G1 <- as.data.frame(dplyr::mutate(as.data.frame(G1)),
                                     row.names=rownames(G1))
G1 <-  G1[!duplicated(G1$SYMBOL),]
G1 <- na.omit(G1)
rownames(G1) <- G1$SYMBOL

G2 <- readRDS("ResultadosAnalasisExpresion_2_Luad.rds")
G2 <- as.data.frame(dplyr::mutate(as.data.frame(G2)),
                    row.names=rownames(G2))
G2 <- na.omit(G2)
G2 <-  G2[!duplicated(G2$SYMBOL),]
rownames(G2) <- G2$SYMBOL

G3 <- readRDS("ResultadosAnalasisExpresion_3_Luad.rds")
G3 <- as.data.frame(dplyr::mutate(as.data.frame(G3)),
                    row.names=rownames(G3))
G3 <- na.omit(G3)
G3 <-  G3[!duplicated(G3$SYMBOL),]
rownames(G3) <- G3$SYMBOL

# 3) Modificamos los DF de cada grupo para luego poder juntarlos en uno

G1 <- G1[HI$GenesInmune,]
G1 <- na.omit(G1)
HI <- HI[rownames(G1),]
G1 <- G1[,c(2,5)]
G1$Cluster <- "Group 1"
G1$SYMBOL <- rownames(G1)
rownames(G1) <- seq(1,nrow(G1),1)

G2 <- G2[HI$GenesInmune,]
G2 <- na.omit(G2)
G2 <- G2[,c(2,5)]
G2$Cluster <- "Group 2"
G2$SYMBOL <- rownames(G2)
rownames(G2) <- seq(1,nrow(G2),1)

G3 <- G3[HI$GenesInmune,]
G3 <- na.omit(G3)
G3 <- G3[,c(2,5)]
G3$Cluster <- "Group 3"
G3$SYMBOL <- rownames(G3)
rownames(G3) <- seq(1,nrow(G3),1)

# 4) generamos la matriz final para poder agrupar los valores del LogFoldChange

Matriz <- data.frame(matrix(0,ncol=3,nrow=nrow(G1)))
rownames(Matriz) <- rownames(HI)
colnames(Matriz) <- c("Group 1","Group 2","Group 3")

GT <- rbind(G1,G2)
GT <- rbind(GT,G3)

# 5) Hacemos el bucle para generar la matriz correspondiente

i <- 1

for(i in seq(1,nrow(HI),1)){
  
  Gen <- HI$GenesInmune[i]
  
  ValorLogFold <- GT[which(GT$SYMBOL==Gen),]
  ValorLogFold <- t(ValorLogFold)
  
  Matriz[i,] <- ValorLogFold[1,]
  
}

# 6) Generamos el heatmap directamente con todos los genes sin hacer el cambio: 

MatrizInmune_Matrix <- matrix(Matriz,ncol=ncol(Matriz),nrow=nrow(Matriz))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(Matriz),nrow=nrow(Matriz))
rownames(MatrizInmune_Matrix) <- rownames(Matriz)
colnames(MatrizInmune_Matrix) <- colnames(Matriz)

# 7) Generamos las anotaciones del HI para poder añadirlas al heatmap: 

AnotacionesHI <- HI[rownames(Matriz),]
AnotacionesHI <- na.omit(AnotacionesHI)
AnotacionesHI_HM <- data.frame(CelularType=AnotacionesHI$TipoCelular)
rownames(AnotacionesHI_HM) <- AnotacionesHI$GenesInmune

# 8) Generamos anotaciones para las columnas: 

Anotaciones <- data.frame(Grupo=c("Group 1","Group 2","Group 3"))
rownames(Anotaciones) <- Anotaciones$Grupo

# 9) Dibujamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(MatrizInmune_Matrix[182:184,], 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F,
         show_rownames =F,
         annotation_row = AnotacionesHI_HM,
         annotation_col = Anotaciones)


# 10) Calculamos la media de los genes: 

  # 10.1) Generamos un vector con el ID de los tipos celulares 

TipoCelularInmune <- unique(HI$TipoCelular)

  # 10.2) Generamos la matriz final. 

MatrizMedias <- data.frame(matrix(0,ncol=ncol(MatrizInmune_Matrix),nrow=length(TipoCelularInmune)))
rownames(MatrizMedias) <- TipoCelularInmune
colnames(MatrizMedias) <- c("Group 1","Group 2","Group 3")

  # 10.3) Generamos el bucle

i <- 1

for(i in seq(1,length(TipoCelularInmune),1)){
  
  TipoCelular_Busqueda <- HI %>% filter(TipoCelular==TipoCelularInmune[i])
  
  MatrixTipoCelularIndexada <- MatrizInmune_Matrix[TipoCelular_Busqueda$hgnc_symbol,]
  MatrixTipoCelularIndexada <- na.omit(MatrixTipoCelularIndexada)
  
  if(length(MatrixTipoCelularIndexada)==3){MatrizMedias[i,] <- MediaPacientes}
  else{
    MediaPacientes <- colMeans(MatrixTipoCelularIndexada)
    
    MatrizMedias[i,] <- MediaPacientes
  }
}

  # 10.4) Generamos el heatmap: 

pheatmap(MatrizMedias, 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F,
         show_rownames =T,
         annotation_col = Anotaciones)

################################################################################################################

# 1) Cargamos los genes que forman parte de mi huella: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$GenesInmune),]
rownames(HI) <- HI$GenesInmune

# 2) Cargamos los resultados del analisis y los convertimos en DF 

G1 <- readRDS("ResultadosAnalasisExpresion_1_Lusc.rds")
G1 <- as.data.frame(dplyr::mutate(as.data.frame(G1)),
                    row.names=rownames(G1))
G1 <-  G1[!duplicated(G1$SYMBOL),]
G1 <- na.omit(G1)
rownames(G1) <- G1$SYMBOL

G2 <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
G2 <- as.data.frame(dplyr::mutate(as.data.frame(G2)),
                    row.names=rownames(G2))
G2 <- na.omit(G2)
G2 <-  G2[!duplicated(G2$SYMBOL),]
rownames(G2) <- G2$SYMBOL

G3 <- readRDS("ResultadosAnalasisExpresion_3_Lusc.rds")
G3 <- as.data.frame(dplyr::mutate(as.data.frame(G3)),
                    row.names=rownames(G3))
G3 <- na.omit(G3)
G3 <-  G3[!duplicated(G3$SYMBOL),]
rownames(G3) <- G3$SYMBOL

# 3) Modificamos los DF de cada grupo para luego poder juntarlos en uno

G1 <- G1[HI$GenesInmune,]
G1 <- na.omit(G1)
HI <- HI[rownames(G1),]
G1 <- G1[,c(2,5)]
G1$Cluster <- "Group 1"
G1$SYMBOL <- rownames(G1)
rownames(G1) <- seq(1,nrow(G1),1)

G2 <- G2[HI$GenesInmune,]
G2 <- na.omit(G2)
G2 <- G2[,c(2,5)]
G2$Cluster <- "Group 2"
G2$SYMBOL <- rownames(G2)
rownames(G2) <- seq(1,nrow(G2),1)

G3 <- G3[HI$GenesInmune,]
G3 <- na.omit(G3)
G3 <- G3[,c(2,5)]
G3$Cluster <- "Group 3"
G3$SYMBOL <- rownames(G3)
rownames(G3) <- seq(1,nrow(G3),1)

# 4) generamos la matriz final para poder agrupar los valores del LogFoldChange

Matriz <- data.frame(matrix(0,ncol=3,nrow=nrow(G1)))
rownames(Matriz) <- rownames(HI)
colnames(Matriz) <- c("Group 1","Group 2","Group 3")

GT <- rbind(G1,G2)
GT <- rbind(GT,G3)

# 5) Hacemos el bucle para generar la matriz correspondiente

i <- 1

for(i in seq(1,nrow(HI),1)){
  
  Gen <- HI$GenesInmune[i]
  
  ValorLogFold <- GT[which(GT$SYMBOL==Gen),]
  ValorLogFold <- t(ValorLogFold)
  
  Matriz[i,] <- ValorLogFold[1,]
  
}

# 6) Generamos el heatmap directamente con todos los genes sin hacer el cambio: 

MatrizInmune_Matrix <- matrix(Matriz,ncol=ncol(Matriz),nrow=nrow(Matriz))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(Matriz),nrow=nrow(Matriz))
rownames(MatrizInmune_Matrix) <- rownames(Matriz)
colnames(MatrizInmune_Matrix) <- colnames(Matriz)

# 7) Generamos las anotaciones del HI para poder añadirlas al heatmap: 

AnotacionesHI <- HI[rownames(Matriz),]
AnotacionesHI <- na.omit(AnotacionesHI)
AnotacionesHI_HM <- data.frame(CelularType=AnotacionesHI$TipoCelular)
rownames(AnotacionesHI_HM) <- AnotacionesHI$GenesInmune

# 8) Generamos anotaciones para las columnas: 

Anotaciones <- data.frame(Grupo=c("Group 1","Group 2","Group 3"))
rownames(Anotaciones) <- Anotaciones$Grupo

# 9) Dibujamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(MatrizInmune_Matrix[c(178,179,180),], 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F,
         show_rownames =F,
         annotation_row = AnotacionesHI_HM,
         annotation_col = Anotaciones)


# 10) Calculamos la media de los genes: 

# 10.1) Generamos un vector con el ID de los tipos celulares 

TipoCelularInmune <- unique(HI$TipoCelular)

# 10.2) Generamos la matriz final. 

MatrizMedias <- data.frame(matrix(0,ncol=ncol(MatrizInmune_Matrix),nrow=length(TipoCelularInmune)))
rownames(MatrizMedias) <- TipoCelularInmune
colnames(MatrizMedias) <- c("Group 1","Group 2","Group 3")

# 10.3) Generamos el bucle

i <- 1

for(i in seq(1,length(TipoCelularInmune),1)){
  
  TipoCelular_Busqueda <- HI %>% filter(TipoCelular==TipoCelularInmune[i])
  
  MatrixTipoCelularIndexada <- MatrizInmune_Matrix[TipoCelular_Busqueda$hgnc_symbol,]
  MatrixTipoCelularIndexada <- na.omit(MatrixTipoCelularIndexada)
  
  if(length(MatrixTipoCelularIndexada)==3){MatrizMedias[i,] <- MediaPacientes}
  else{
    MediaPacientes <- colMeans(MatrixTipoCelularIndexada)
    
    MatrizMedias[i,] <- MediaPacientes
  }
}

# 10.4) Generamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-1,1,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(MatrizMedias, 
         color = Pallete1, 
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F,
         show_rownames =T,
         annotation_col = Anotaciones)
