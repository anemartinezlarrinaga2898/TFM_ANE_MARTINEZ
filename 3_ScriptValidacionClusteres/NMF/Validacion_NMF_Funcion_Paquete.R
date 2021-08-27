# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO:
#SCRIPT para la validacion de los clusteres generados por NMF: 
# Validacion del RANK mediante la funcion que nos da NMF : 

set.seed(123456789)
library(NMF)
library(ggplot2)
library(tidyverse)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/")
muestras <- readRDS('unionLUAD.txt')
muestras <- muestras+0.5
muestras <- t(muestras)

NumeroCluster_inicial <- 2
NumeroCluster_final <- 10
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)

EstimationOfRank <- nmfEstimateRank(muestras,range=clusteresTotales,nrun=100,method="brunet",seed="random")


