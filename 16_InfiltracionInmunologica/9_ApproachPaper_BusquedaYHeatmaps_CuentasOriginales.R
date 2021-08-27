# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 21-07-2021

# OBJETIVO:Infiltracion inmunologica del approach del paper. Tengo que encontrar X genes que tengo en un data.frame y generar el heatmap. 

##########################################################################################################################################
directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
##########################################################################################################################################

#·························· L U A D ············································

Matriz <- readRDS("MatrizLuad_Anotada_Escalada.rds")

# 1) Organizamos la matriz en el orden de los pacientes correspondiente. 

G1 <- readRDS("K3.1_NMF_Luad.rds")
G2 <- readRDS("K3.2_NMF_Luad.rds")
G3 <- readRDS("K3.3_NMF_Luad.rds")

Matriz_G1 <- Matriz[,G1$Sample_ID]
Matriz_G2 <- Matriz[,G2$Sample_ID]
Matriz_G3 <- Matriz[,G3$Sample_ID]

Matriz <- cbind(Matriz_G1,Matriz_G2)
Matriz <- cbind(Matriz,Matriz_G3)

# 2) Generamos las anotaciones para el heatmap: 

Anotaciones <- data.frame(Pacientes=c(G1$Sample_ID,G2$Sample_ID,G3$Sample_ID),
                          Cluster=c(G1$nmf_subtypes,G2$nmf_subtypes,G3$nmf_subtypes))
Anotaciones_HM <- data.frame(Cluster=Anotaciones$Cluster)
rownames(Anotaciones_HM) <- Anotaciones$Pacientes

# 3) Cargamos los genes correspondientes que nos interesa buscar: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$hgnc_symbol),]
rownames(HI) <- HI$hgnc_symbol

# 4) Indexamos por los genes que nos interesan del sistema inmune 

MatrizInmune <- Matriz[HI$hgnc_symbol,]
MatrizInmune <- na.omit(MatrizInmune)

# 5) Indexamos en las anotaciones del sistema inmune, en el DF para tener el mismo numero de genes 

AnotacionesHI <- HI[rownames(MatrizInmune),]
AnotacionesHI <- na.omit(AnotacionesHI)
AnotacionesHI_HM <- data.frame(CelularType=AnotacionesHI$TipoCelular)
rownames(AnotacionesHI_HM) <- AnotacionesHI$hgnc_symbol

MatrizInmune <- Matriz[rownames(AnotacionesHI),]

# 6) Convertimos el data.frame a una matriz para poder hacer el heatmap 

MatrizInmune_Matrix <- matrix(MatrizInmune,ncol=ncol(MatrizInmune),nrow=nrow(MatrizInmune))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(MatrizInmune),nrow=nrow(MatrizInmune))
rownames(MatrizInmune_Matrix) <- rownames(MatrizInmune)
colnames(MatrizInmune_Matrix) <- colnames(MatrizInmune)

# 7) Dibujamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-1,1,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(MatrizInmune_Matrix, 
         color = Pallete1, 
         cluster_rows =T ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = T,
         show_rownames =F,
         annotation_col = Anotaciones_HM,
         annotation_row = AnotacionesHI_HM)


###############################################################################################

#························· L U A D ·································

# Tengo que encontrar la media de los genes que definen un tipo celular. 

Matriz <- readRDS("MatrizLuad_Anotada_Escalada.rds")

# 1) Organizamos la matriz en el orden de los pacientes correspondiente. 

G1 <- readRDS("K3.1_NMF_Luad.rds")
G2 <- readRDS("K3.2_NMF_Luad.rds")
G3 <- readRDS("K3.3_NMF_Luad.rds")

Matriz_G1 <- Matriz[,G1$Sample_ID]
Matriz_G2 <- Matriz[,G2$Sample_ID]
Matriz_G3 <- Matriz[,G3$Sample_ID]

Matriz <- cbind(Matriz_G1,Matriz_G2)
Matriz <- cbind(Matriz,Matriz_G3)

# 2) Cargamos los genes correspondientes que nos interesa buscar: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$hgnc_symbol),]
rownames(HI) <- HI$hgnc_symbol

# 3) Indexamos por los genes que nos interesan del sistema inmune 

MatrizInmune <- Matriz[HI$hgnc_symbol,]
MatrizInmune <- na.omit(MatrizInmune)

# 5) Generamos vectores con el nombre de las celulas y con el nombre de los pacientes 

TipoCelularInmune <- HI$TipoCelular
TipoCelularInmune <- TipoCelularInmune[!duplicated(TipoCelularInmune)]

Pacientes <- colnames(MatrizInmune)

# 6) Generamos el objeto final de la matriz donde las lineas celulares seran la columna y los pacientes las filas 

MatrizFinalInmune <- data.frame(matrix(0,ncol=length(TipoCelularInmune),nrow=length(Pacientes)))
rownames(MatrizFinalInmune) <- Pacientes 
colnames(MatrizFinalInmune) <- TipoCelularInmune

# 7) Tenemos que generar un bucle que me busque los genes que definen a un grupo, e ir calculando la media. 

i <- 1

for(i in seq(1,length(TipoCelularInmune),1)){
  
  TipoCelular_Busqueda <- HI %>% filter(TipoCelular==TipoCelularInmune[i])
  
  MatrixTipoCelularIndexada <- MatrizInmune[TipoCelular_Busqueda$hgnc_symbol,]
  MatrixTipoCelularIndexada <- na.omit(MatrixTipoCelularIndexada)
  
  MediaPacientes <- colMeans(MatrixTipoCelularIndexada)
  names(MediaPacientes) <- Pacientes
  
  MatrizFinalInmune[,i] <- MediaPacientes
  
  
}

# 8) Convertimos el data.frame a matrix 

MatrizInmune_Matrix <- matrix(MatrizFinalInmune,ncol=ncol(MatrizFinalInmune),nrow=nrow(MatrizFinalInmune))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(MatrizFinalInmune),nrow=nrow(MatrizFinalInmune))
rownames(MatrizInmune_Matrix) <- rownames(MatrizFinalInmune)
colnames(MatrizInmune_Matrix) <- colnames(MatrizFinalInmune)

# 9) Generamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-1,1,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(t(MatrizInmune_Matrix), 
         color = Pallete1, 
         cluster_rows =T ,
         breaks = breaks,
         cluster_cols = F,
         show_colnames = F,
         show_rownames =T,
         annotation_col = Anotaciones_HM)

############################################################################################################

#·························· L U S C ············································

Matriz <- readRDS("MatrizLusc_Anotada_Escalada.rds")

# 1) Organizamos la matriz en el orden de los pacientes correspondiente. 

G1 <- readRDS("K3.1_NMF_Lusc.rds")
G2 <- readRDS("K3.2_NMF_Lusc.rds")
G3 <- readRDS("K3.3_NMF_Lusc.rds")

Matriz_G1 <- Matriz[,G1$Sample_ID]
Matriz_G2 <- Matriz[,G2$Sample_ID]
Matriz_G3 <- Matriz[,G3$Sample_ID]

Matriz <- cbind(Matriz_G1,Matriz_G2)
Matriz <- cbind(Matriz,Matriz_G3)

# 2) Generamos las anotaciones para el heatmap: 

Anotaciones <- data.frame(Pacientes=c(G1$Sample_ID,G2$Sample_ID,G3$Sample_ID),
                          Cluster=c(G1$nmf_subtypes,G2$nmf_subtypes,G3$nmf_subtypes))
Anotaciones_HM <- data.frame(Cluster=Anotaciones$Cluster)
rownames(Anotaciones_HM) <- Anotaciones$Pacientes

# 3) Cargamos los genes correspondientes que nos interesa buscar: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$hgnc_symbol),]
rownames(HI) <- HI$hgnc_symbol

# 4) Indexamos por los genes que nos interesan del sistema inmune 

MatrizInmune <- Matriz[HI$hgnc_symbol,]
MatrizInmune <- na.omit(MatrizInmune)

# 5) Indexamos en las anotaciones del sistema inmune, en el DF para tener el mismo numero de genes 

AnotacionesHI <- HI[rownames(MatrizInmune),]
AnotacionesHI <- na.omit(AnotacionesHI)
AnotacionesHI_HM <- data.frame(CelularType=AnotacionesHI$TipoCelular)
rownames(AnotacionesHI_HM) <- AnotacionesHI$hgnc_symbol

MatrizInmune <- Matriz[rownames(AnotacionesHI),]

# 6) Convertimos el data.frame a una matriz para poder hacer el heatmap 

MatrizInmune_Matrix <- matrix(MatrizInmune,ncol=ncol(MatrizInmune),nrow=nrow(MatrizInmune))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(MatrizInmune),nrow=nrow(MatrizInmune))
rownames(MatrizInmune_Matrix) <- rownames(MatrizInmune)
colnames(MatrizInmune_Matrix) <- colnames(MatrizInmune)

# 7) Dibujamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-1,1,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(MatrizInmune_Matrix, 
         color = Pallete1, 
         cluster_rows =T ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F,
         show_rownames =F,
         annotation_col = Anotaciones_HM,
         annotation_row = AnotacionesHI_HM)



#························· L U S C ·································

# Tengo que encontrar la media de los genes que definen un tipo celular. 

Matriz <- readRDS("MatrizLusc_Anotada_Escalada.rds")

# 1) Organizamos la matriz en el orden de los pacientes correspondiente. 

G1 <- readRDS("K3.1_NMF_Lusc.rds")
G2 <- readRDS("K3.2_NMF_Lusc.rds")
G3 <- readRDS("K3.3_NMF_Lu.rds")

Matriz_G1 <- Matriz[,G1$Sample_ID]
Matriz_G2 <- Matriz[,G2$Sample_ID]
Matriz_G3 <- Matriz[,G3$Sample_ID]

Matriz <- cbind(Matriz_G1,Matriz_G2)
Matriz <- cbind(Matriz,Matriz_G3)

# 2) Cargamos los genes correspondientes que nos interesa buscar: 

HI <- readRDS("HuellasSistemaInmune_DF.rds")
HI <- HI[,-3]
HI <- HI[!duplicated(HI$hgnc_symbol),]
rownames(HI) <- HI$hgnc_symbol

# 3) Indexamos por los genes que nos interesan del sistema inmune 

MatrizInmune <- Matriz[HI$hgnc_symbol,]
MatrizInmune <- na.omit(MatrizInmune)

# 5) Generamos vectores con el nombre de las celulas y con el nombre de los pacientes 

TipoCelularInmune <- HI$TipoCelular
TipoCelularInmune <- TipoCelularInmune[!duplicated(TipoCelularInmune)]

Pacientes <- colnames(MatrizInmune)

# 6) Generamos el objeto final de la matriz donde las lineas celulares seran la columna y los pacientes las filas 

MatrizFinalInmune <- data.frame(matrix(0,ncol=length(TipoCelularInmune),nrow=length(Pacientes)))
rownames(MatrizFinalInmune) <- Pacientes 
colnames(MatrizFinalInmune) <- TipoCelularInmune

# 7) Tenemos que generar un bucle que me busque los genes que definen a un grupo, e ir calculando la media. 

i <- 1

for(i in seq(1,length(TipoCelularInmune),1)){
  
  TipoCelular_Busqueda <- HI %>% filter(TipoCelular==TipoCelularInmune[i])
  
  MatrixTipoCelularIndexada <- MatrizInmune[TipoCelular_Busqueda$hgnc_symbol,]
  MatrixTipoCelularIndexada <- na.omit(MatrixTipoCelularIndexada)
  
  MediaPacientes <- colMeans(MatrixTipoCelularIndexada)
  names(MediaPacientes) <- Pacientes
  
  MatrizFinalInmune[,i] <- MediaPacientes
  
  
}

# 8) Convertimos el data.frame a matrix 

MatrizInmune_Matrix <- matrix(MatrizFinalInmune,ncol=ncol(MatrizFinalInmune),nrow=nrow(MatrizFinalInmune))
MatrizInmune_Matrix <- as.numeric(unlist(MatrizInmune_Matrix))
MatrizInmune_Matrix<- matrix(MatrizInmune_Matrix,ncol=ncol(MatrizFinalInmune),nrow=nrow(MatrizFinalInmune))
rownames(MatrizInmune_Matrix) <- rownames(MatrizFinalInmune)
colnames(MatrizInmune_Matrix) <- colnames(MatrizFinalInmune)

# 9) Generamos el heatmap: 

Colors <- c("green","black","red")
breaks <- seq(-1,1,by=0.01)
Pallete1 <- colorRampPalette(Colors)(length(breaks))

pheatmap(t(MatrizInmune_Matrix), 
         color = Pallete1, 
         cluster_rows =T ,
         breaks = breaks,
         cluster_cols = F,
         show_colnames = F,
         show_rownames =T,
         annotation_col = Anotaciones_HM)
