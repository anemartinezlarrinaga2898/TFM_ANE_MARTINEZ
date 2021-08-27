# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 24-06-2021

#OBJETIVO: Obtencion de la matriz de correlacion de mis lineas celulares para ver si representado de otra forma podemos obtener una mejor 
#visualizacion. 

###########################################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(pheatmap)
library(igraph)
library(visNetwork)

###########################################################################################################################################################

MatrizExpresionLineas <- readRDS("Matriz_Huellas_Lineas.rds")
MediaExpresion <- colMeans(MatrizExpresionLineas)
names(MediaExpresion) <- colnames(MatrizExpresionLineas)

MatrizAdyacencia <- matrix(0,ncol=ncol(MatrizExpresionLineas),nrow=nrow(MatrizExpresionLineas))
colnames(MatrizAdyacencia) <- colnames(MatrizExpresionLineas)
rownames(MatrizAdyacencia) <- rownames(MatrizExpresionLineas)

i <- 1 #Columna
j <- 1 #Fila

SecuenciaFilas <- seq(1,nrow(MatrizExpresionLineas),1)

for(i in seq_along(MediaExpresion)){
  MediaLinea <- MediaExpresion[i]
  
  for( j in seq_along(SecuenciaFilas)){
    X <- MatrizExpresionLineas[j,i]
    
    #Si es mayor que la media ponemos un 1 si es menor dejamos la posicion porque tenemos que la matriz final es de 0 y 1. 
    if (X >= MediaLinea){MatrizAdyacencia[j,i] <- 1}
    else(next)
  }
}

saveRDS(MatrizAdyacencia,"MatrizAdyacencia.rds")

Correlacion <- cor(MatrizAdyacencia)
pheatmap(Correlacion)

Cor2 <- Correlacion
Cor2[Cor2>0.5] <- 0.5
diag(Cor2) <- 0
pheatmap(Cor2)

MatrizAdyavencia <- Correlacion
MatrizAdyavencia[MatrizAdyavencia<0.2] <- 0
MatrizAdyavencia[MatrizAdyavencia>0.2] <- 1
diag(MatrizAdyavencia) <- 0

aa <- graph_from_adjacency_matrix(MatrizAdyavencia,mode="undirected")
visnet_obj <- toVisNetworkData(aa)
visNetwork(nodes = visnet_obj$nodes, 
           edges = visnet_obj$edges)
network <- graph_from_adjacency_matrix(MatrizAdyavencia,mode="undirected")
V(network) #Son los vertices de nuestro 

#Eliminar los verticies que sean cero: 
net2 <- delete_vertices(network,V(network)[degree(network)==0])

visnet_obj2 <- toVisNetworkData(net2)
visNetwork(nodes = visnet_obj2$nodes, 
           edges = visnet_obj2$edges)

plot(cluster_fast_greedy(net2),net2)
