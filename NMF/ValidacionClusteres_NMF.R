# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO: 
#SCRIPT para la validacion de los clusteres generados por NMF: 

# Cluster algoritm: NMF
# Resample scheme: resample with sample we will take 70% of the sample in each run
# Number of iteratiosn: number of times that we will resample 
# Number of clusters: the number of clusters taht we are going to try with

##############################################################################################################

set.seed(1234)
library(NMF)
library(ggplot2)
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/")
muestras <- readRDS('LuadFiltrado.rds')
muestras <- muestras+0.5

NumeroCluster_inicial <- 2
NumeroCluster_final <- 10
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)
NumeroIteraccion_inicial <- 1
NumeroIteraccion_final <- 100
iteraccionesTotales <- seq(NumeroIteraccion_inicial,NumeroIteraccion_final,1)

ErrorResidual <- matrix(data=0,nrow=NumeroCluster_final,ncol=NumeroIteraccion_final)
MediaErrorResidual <- vector(mode="numeric",length=NumeroCluster_final)
MatricesConsenso <- list()
valorcophoneti <- list()

k <- 1
h <- 1
for (k in clusteresTotales) {
  MatricesDeConectividad <- list()
  cluster <- clusteresTotales[k]
  for (h in iteraccionesTotales){
    D <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ] #Aqui hacemos el boostrapping
    MenK <- nmf(D,k,method="brunet",seed="random")
    conectividad <- connectivity(Menk)
    MatricesDeConectividad[[h]] <- conectividad
    H <- coef(Menk)
    W <- basis(Menk)
    Error <- abs(D-W*H)
    ErrorResidual[k,h] <- Error
    
  }
  consenso <- consensus(MatricesDeConectividad)
  MatricesConsenso[[k]] <- consenso
  cophonetic <- cophcor(consenso)
  valorcophoneti[[k]] <- cophonetic
  MediaErrorResidual[k] <- mean(ErrorResidual[k,])
} 

