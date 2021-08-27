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
library(igraph)
library(visNetwork)
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

AnotacionesPacientesLuad_HP <- data.frame(Group=AnotacionesPacientesLuad$Cluster)
rownames(AnotacionesPacientesLuad_HP) <- AnotacionesPacientesLuad$Pacientes

AnotacionesPacientesLusc <- data.frame(Pacientes=c(G1_Lusc$Sample_ID,G2_Lusc$Sample_ID,G3_Lusc$Sample_ID),
                                       Cluster=c(G1_Lusc$nmf_subtypes,G2_Lusc$nmf_subtypes,G3_Lusc$nmf_subtypes))
AnotacionesPacientesLusc$Cluster <- paste("Group",AnotacionesPacientesLusc$Cluster)
rownames(AnotacionesPacientesLusc) <- AnotacionesPacientesLusc$Pacientes

AnotacionesPacientesLusc_HP <- data.frame(Group=AnotacionesPacientesLusc$Cluster)
rownames(AnotacionesPacientesLusc_HP) <- AnotacionesPacientesLusc$Pacientes

# A.4) Ahora vamos a generar los heatmap: 

# ··············· L U A D ·································

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

MLuad_FT_Matrix <- as.numeric(unlist(MatrizSinCeros))
MLuad_FT_Matrix<- matrix(MLuad_FT_Matrix,ncol=ncol(MatrizSinCeros),nrow=nrow(MatrizSinCeros))
rownames(MLuad_FT_Matrix) <- rownames(MatrizSinCeros)
colnames(MLuad_FT_Matrix) <- colnames(MatrizSinCeros)

colnames(MLuad_FT_Matrix) <- str_replace_all(colnames(MLuad_FT_Matrix),"\\.","-")


pheatmap(t(MLuad_FT_Matrix),
         color = Pallete1,
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F ,
         show_rownames =F,
         annotation_row =AnotacionesPacientesLuad_HP)

#············ L U S C ································
         
MLusc_FT_Matrix <- as.numeric(unlist(MLusc_FT))
MLusc_FT_Matrix<- matrix(MLusc_FT_Matrix,ncol=ncol(MLusc_FT),nrow=nrow(MLusc_FT))
rownames(MLusc_FT_Matrix) <- rownames(MLusc_FT)
colnames(MLusc_FT_Matrix) <- colnames(MLusc_FT)

Varianzas <- rowVars(MLusc_FT_Matrix)
names(Varianzas) <- rownames(MLusc_FT_Matrix)
hist(Varianzas)

VarianzasSinCeros <- Varianzas[which(Varianzas>0.2)]

MLusc_FT_Matrix_DF<- data.frame(MLusc_FT_Matrix)

MatrizSinCeros <- MLusc_FT_Matrix_DF[names(VarianzasSinCeros),]

MLusc_FT_Matrix <- as.numeric(unlist(MatrizSinCeros))
MLusc_FT_Matrix<- matrix(MLusc_FT_Matrix,ncol=ncol(MatrizSinCeros),nrow=nrow(MatrizSinCeros))
rownames(MLusc_FT_Matrix) <- rownames(MatrizSinCeros)
colnames(MLusc_FT_Matrix) <- colnames(MatrizSinCeros)

colnames(MLusc_FT_Matrix) <- str_replace_all(colnames(MLusc_FT_Matrix),"\\.","-")


pheatmap(t(MLusc_FT_Matrix),
         color = Pallete1,
         cluster_rows =F ,
         breaks = breaks,
         cluster_cols = T,
         show_colnames = F ,
         show_rownames =F,
         annotation_row =AnotacionesPacientesLusc_HP)

##############################################################################################################################

# MATRICES DE CORRELACION: Para comparar con la de las lineas celulares. 

MatrizLuad <- t(MatrizLuad)
MatrizLusc <- t(MatrizLusc)

MediaExpresion_Luad <- colMeans(MatrizLuad)
names(MediaExpresion_Luad) <- colnames(MatrizLuad)

MediaExpresion_Lusc <- colMeans(MatrizLusc)
names(MediaExpresion_Lusc) <- colnames(MatrizLusc)

MatrizAdyacencia_Luad <- matrix(0,ncol=ncol(MatrizLuad),nrow=nrow(MatrizLuad))
colnames(MatrizAdyacencia_Luad) <- colnames(MatrizLuad)
rownames(MatrizAdyacencia_Luad) <- rownames(MatrizLuad)

i <- 1 #Columna
j <- 1 #Fila

SecuenciaFilas <- seq(1,nrow(MatrizLuad),1)

for(i in seq_along(MediaExpresion_Luad)){
  MediaLinea <- MediaExpresion_Luad[i]
  
  for( j in seq_along(SecuenciaFilas)){
    X <- MatrizLuad[j,i]
    
    #Si es mayor que la media ponemos un 1 si es menor dejamos la posicion porque tenemos que la matriz final es de 0 y 1. 
    if (X >= MediaLinea){MatrizAdyacencia_Luad[j,i] <- 1}
    else(next)
  }
}

MatrizAdyacencia_Lusc <- matrix(0,ncol=ncol(MatrizLusc),nrow=nrow(MatrizLusc))
colnames(MatrizAdyacencia_Lusc) <- colnames(MatrizLusc)
rownames(MatrizAdyacencia_Lusc) <- rownames(MatrizLusc)

i <- 1 #Columna
j <- 1 #Fila

SecuenciaFilas <- seq(1,nrow(MatrizLusc),1)

for(i in seq_along(MediaExpresion_Luad)){
  MediaLinea <- MediaExpresion_Lusc[i]
  
  for( j in seq_along(SecuenciaFilas)){
    X <- MatrizLusc[j,i]
    
    #Si es mayor que la media ponemos un 1 si es menor dejamos la posicion porque tenemos que la matriz final es de 0 y 1. 
    if (X >= MediaLinea){MatrizAdyacencia_Lusc[j,i] <- 1}
    else(next)
  }
}

Correlacion_Luad <- cor(MatrizAdyacencia_Luad)

pheatmap(Correlacion_Luad)

Cor2_Luad <- Correlacion_Luad
Cor2_Luad[Cor2_Luad>0.5] <- 0.5
diag(Cor2_Luad) <- 0
pheatmap(Cor2_Luad)

Correlacion_Lusc <- cor(MatrizAdyacencia_Lusc)

pheatmap(Correlacion_Lusc)

Cor2_Lusc <- Correlacion_Lusc
Cor2_Lusc[Cor2_Lusc>0.5] <- 0.5
diag(Cor2_Lusc) <- 0
pheatmap(Cor2_Lusc)


aa <- graph_from_adjacency_matrix(Correlacion_Luad,mode="undirected")
visnet_obj <- toVisNetworkData(aa)
visNetwork(nodes = visnet_obj$nodes, 
           edges = visnet_obj$edges)
network <- graph_from_adjacency_matrix(Correlacion_Luad,mode="undirected")
V(network) #Son los vertices de nuestro 

#Eliminar los verticies que sean cero: 
net2 <- delete_vertices(network,V(network)[degree(network)==0])

visnet_obj2 <- toVisNetworkData(net2)
visNetwork(nodes = visnet_obj2$nodes, 
           edges = visnet_obj2$edges)

plot(cluster_fast_greedy(net2),net2)
