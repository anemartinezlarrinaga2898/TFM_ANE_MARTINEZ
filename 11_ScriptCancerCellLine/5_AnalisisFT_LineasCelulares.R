# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-06-2021

#OBJETIVO:

    #1) Obtencion delos FT de mis lineas

#######################################################################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(pheatmap)
library(igraph)
library(visNetwork)
#######################################################################################################################################

# 1) Anotar mis lineas celulares para luego poder encontrar mis FT: 

Matriz <- readRDS("MatrizLineasCelulares_TotalGenes.rds")
FT <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
FT <- FT[!duplicated(FT$Symbol),]
rownames(FT) <- FT$ensembl_gene_id
colnames(FT)[4] <- "ensembl_gene_id"

Matriz_FT <- na.omit(Matriz[FT$ensembl_gene_id,])
FT_Lineas <- FT[rownames(Matriz_FT),]

FT_Lineas$ensembl_gene_id==rownames(Matriz_FT)

rownames(Matriz_FT) <- FT_Lineas$Symbol

# 2) Ver si con los FT clusterizan de algun modo determinado: 

MatrizCorrelacion <- cor(Matriz_FT)
pheatmap::pheatmap(MatrizCorrelacion)

Cor2 <- MatrizCorrelacion
Cor2[Cor2>0.7] <- 0.5
diag(Cor2) <- 0
pheatmap::pheatmap(Cor2)

MatrizAdyavencia <- MatrizCorrelacion
MatrizAdyavencia[MatrizAdyavencia<0.85] <- 0
MatrizAdyavencia[MatrizAdyavencia>0.85] <- 1
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

# 3) Localizamos los FT unicos de cada grupo: 

FT_G1 <- readRDS("Grupo1_LUSC_FT_Unicos.rds")
FT_G2 <- readRDS("Grupo2_LUSC_FT_Unicos.rds")
FT_G3 <- readRDS("Grupo3_LUSC_FT_Unicos.rds")

Matriz_FT_G1 <- na.omit(Matriz_FT[FT_G1$SYMBOL,])
Matriz_FT_G2 <- na.omit(Matriz_FT[FT_G2$SYMBOL,])
Matriz_FT_G3 <- na.omit(Matriz_FT[FT_G3$SYMBOL,])

# 4) Generacion del heatmap del FT_G1

Matriz_FT_G1 <- scale(as.matrix(Matriz_FT_G1))
pheatmap(Matriz_FT_G1,cluster_rows = F,cluster_cols = F)

Matriz_FT_G2 <- scale(as.matrix(Matriz_FT_G2))
pheatmap(Matriz_FT_G2,cluster_rows = F,cluster_cols = F)
