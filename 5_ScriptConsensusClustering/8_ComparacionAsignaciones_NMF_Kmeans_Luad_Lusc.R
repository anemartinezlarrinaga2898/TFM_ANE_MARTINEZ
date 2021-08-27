# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 06-04-2021


#OBJETIVO: comparar como se divide el grupo 1,2,3 y de k=3 en NMF tanto de Luad como Lusc en Kmeans para saber que separacion de grupos coger. 

####################################################################################################

#································L U A D ·································

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(tidyverse)

GruposKmeans_Luad <- readRDS("AsignacionCluster_Kmeans_Luad.rds")
K3_KMENAS_Luad <- GruposKmeans_Luad %>% filter(k==3)
K3.3_KMENAS_Luad <- K3_KMENAS_Luad %>% filter(cluster==3)
K3.2_KMENAS_Luad <- K3_KMENAS_Luad %>% filter(cluster==2)
K3.1_KMENAS_Luad <- K3_KMENAS_Luad %>% filter(cluster==1)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/Grupos/NMF/")

K3.1_NMF_Luad <- readRDS("K3.1_NMF_Luad.rds")
K3.2_NMF_Luad <- readRDS("K3.2_NMF_Luad.rds")
K3.3_NMF_Luad <- readRDS("K3.3_NMF_Luad.rds")


IDK3.1<- str_replace_all(K3.1_NMF_Luad$Sample_ID,"-",".")
K3.1_NMF_Luad <- cbind(K3.1_NMF_Luad,IDK3.1)
K3.1_NMF_Luad <- K3.1_NMF_Luad[,-1]
colnames(K3.1_NMF_Luad)[3] <- "item"

IDK3.2<- str_replace_all(K3.2_NMF_Luad$Sample_ID,"-",".")
K3.2_NMF_Luad <- cbind(K3.2_NMF_Luad,IDK3.2)
K3.2_NMF_Luad <- K3.2_NMF_Luad[,-1]
colnames(K3.2_NMF_Luad)[3] <- "item"

IDK3.3<- str_replace_all(K3.3_NMF_Luad$Sample_ID,"-",".")
K3.3_NMF_Luad <- cbind(K3.3_NMF_Luad,IDK3.3)
K3.3_NMF_Luad <- K3.3_NMF_Luad[,-1]
colnames(K3.3_NMF_Luad)[3] <- "item"


# Grupo 3.1 NMF en Kmenas K3: 

K3.1_NMF_KMEANS_Luad <- inner_join(K3_KMENAS_Luad,K3.1_NMF_Luad,by="item")

# Grupo 3.2 NMF en Kmenas K3: 

K3.2_NMF_KMEANS_Luad <- inner_join(K3_KMENAS_Luad,K3.2_NMF_Luad,by="item")

# Grupo 3.3 NMF en Kmenas K3: 

K3.3_NMF_KMEANS_Luad <- inner_join(K3_KMENAS_Luad,K3.3_NMF_Luad,by="item")


###################################################################################################

#····················· L U S C·········································

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(tidyverse)

GruposKmeans_Lusc <- readRDS("AsignacionCluster_Kmeans_Lusc.rds")
K3_KMENAS_Lusc <- GruposKmeans_Lusc %>% filter(k==3)
K3.3_KMENAS_Lusc <- K3_KMENAS_Lusc %>% filter(cluster==3)
K3.2_KMENAS_Lusc <- K3_KMENAS_Lusc %>% filter(cluster==2)
K3.1_KMENAS_Lusc <- K3_KMENAS_Lusc %>% filter(cluster==1)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/Grupos/NMF/")

K3.1_NMF_Lusc <- readRDS("K3.1_NMF_Lusc.rds")
K3.2_NMF_Lusc <- readRDS("K3.2_NMF_Lusc.rds")
K3.3_NMF_Lusc <- readRDS("K3.3_NMF_Lusc.rds")

IDK3.1<- str_replace_all(K3.1_NMF_Lusc$Sample_ID,"-",".")
K3.1_NMF_Lusc <- cbind(K3.1_NMF_Lusc,IDK3.1)
K3.1_NMF_Lusc <- K3.1_NMF_Lusc[,-1]
colnames(K3.1_NMF_Lusc)[3] <- "item"

IDK3.2<- str_replace_all(K3.2_NMF_Lusc$Sample_ID,"-",".")
K3.2_NMF_Lusc <- cbind(K3.2_NMF_Lusc,IDK3.2)
K3.2_NMF_Lusc <- K3.2_NMF_Lusc[,-1]
colnames(K3.2_NMF_Lusc)[3] <- "item"

IDK3.3<- str_replace_all(K3.3_NMF_Lusc$Sample_ID,"-",".")
K3.3_NMF_Lusc <- cbind(K3.3_NMF_Lusc,IDK3.3)
K3.3_NMF_Lusc <- K3.3_NMF_Lusc[,-1]
colnames(K3.3_NMF_Lusc)[3] <- "item"

# Grupo 3.1 NMF en Kmenas K3: 

K3.1_NMF_KMEANS_Lusc <- inner_join(K3_KMENAS_Lusc,K3.1_NMF_Lusc,by="item")

# Grupo 3.2 NMF en Kmenas K3: 

K3.2_NMF_KMEANS_Lusc <- inner_join(K3_KMENAS_Lusc,K3.2_NMF_Lusc,by="item")

# Grupo 3.3 NMF en Kmenas K3: 

K3.3_NMF_KMEANS_Lusc <- inner_join(K3_KMENAS_Lusc,K3.3_NMF_Lusc,by="item")



