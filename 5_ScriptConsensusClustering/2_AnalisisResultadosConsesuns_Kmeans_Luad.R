# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 06-04-2021


#OBJETIVO: analizar los resultados de como estan distribuidos los pacientes en los diferentes grupos. 

####################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(tidyverse)
#################################################################################

GruposKmeans_Luad <- readRDS("AsignacionCluster_Kmeans_Luad.rds")


#Generamos los data.frame para cada valor de K: 
NumeroCluster <- c(3,4,5)
K3 <- data.frame()
K4 <- data.frame()
K5 <- data.frame()
i <- 1
for (i in seq_along(NumeroCluster)){
  Grupo <- NumeroCluster[i]
  if (Grupo==NumeroCluster[1]){K3 <- GruposKmeans_Luad %>% filter(k==Grupo)}
  else if (Grupo==NumeroCluster[2]){K4 <- GruposKmeans_Luad %>% filter(k==Grupo)} 
  else {K5 <- GruposKmeans_Luad %>% filter(k==Grupo)}
}

#Ordenamos el orden de los pacientes: 

K3 <- K3[order(K3$cluster),]
K4 <- K4[order(K4$cluster),]
K5 <- K5[order(K5$cluster),]

###################################################################################
# ··························· K = 3 ···········································

K3.1 <- K3 %>% filter(cluster==1) # 31 pacientes
K3.2 <- K3 %>% filter(cluster==2) # 19 pacientes
K3.3 <- K3 %>% filter(cluster==3) # 11 pacientes

K3.1 <- K3.1[,c(2,3)]
K3.2 <- K3.2[,c(2,3)]
K3.3 <- K3.3[,c(2,3)]

K3.1 <- K3.1[order(K3.1$item),]
K3.2 <- K3.2[order(K3.2$item),]
K3.3 <- K3.3[order(K3.3$item),]

# G R U P O  1: ································································

K3.1_4 <- inner_join(K4,K3.1,by="item")
colnames(K3.1_4) <- c("k","PertenenciaK4","consensus","PertenenciaK3.1")
Comprobacion_3.1_4 <- table(K3.1_4$PertenenciaK4)
Comprobacion_3.1_4
# Numero de grupo:        2  3 
# Cantidad en cada grupo: 19 12 

K3.1_5 <- inner_join(K5,K3.1,by="item")
colnames(K3.1_5) <- c("k","PertenenciaK5","consensus","PertenenciaK3.1")
Comprobacion_3.1_5 <- table(K3.1_5$PertenenciaK5)
Comprobacion_3.1_5
# Numero de grupo:        1  2 
# Cantidad en cada grupo: 19 12

# G R U P O  2: ································································

K3.2_4 <- inner_join(K4,K3.2,by="item")
colnames(K3.2_4) <- c("k","PertenenciaK4","consensus","PertenenciaK3.2")
Comprobacion_3.2_4 <- table(K3.2_4$PertenenciaK4)
Comprobacion_3.2_4
# Numero de grupo:        1  2
# Cantidad en cada grupo: 18 1 

K3.2_5 <- inner_join(K5,K3.2,by="item")
colnames(K3.2_5) <- c("k","PertenenciaK5","consensus","PertenenciaK3.2")
Comprobacion_3.2_5 <- table(K3.2_5$PertenenciaK5)
Comprobacion_3.2_5
# Numero de grupo:        1 3 5     
# Cantidad en cada grupo: 3 8 8 

# G R U P O  3: ································································

K3.3_4 <- inner_join(K4,K3.3,by="item")
colnames(K3.3_4) <- c("k","PertenenciaK4","consensus","PertenenciaK3.3")
Comprobacion_3.3_4 <- table(K3.3_4$PertenenciaK4)
Comprobacion_3.3_4
# Numero de grupo:        2 4
# Cantidad en cada grupo: 1 10 

K3.3_5 <- inner_join(K5,K3.3,by="item")
colnames(K3.3_5) <- c("k","PertenenciaK5","consensus","PertenenciaK3.3")
Comprobacion_3.3_5 <- table(K3.3_5$PertenenciaK5)
Comprobacion_3.3_5
# Numero de grupo:        1 4 5 
# Cantidad en cada grupo: 1 9 1 
 # K3
###################################################################################
# ··························· K = 4 ···········································

K4.1 <- K4 %>% filter(cluster==1) # 18 pacientes
K4.2 <- K4 %>% filter(cluster==2) # 21 pacientes
K4.3 <- K4 %>% filter(cluster==3) # 12 pacientes
K4.4 <- K4 %>% filter(cluster==4) # 10 pacientes

K4.1 <- K4.1[,c(2,3)]
K4.2 <- K4.2[,c(2,3)]
K4.3 <- K4.3[,c(2,3)]

K4.1 <- K4.1[order(K4.1$item),]
K4.2 <- K4.2[order(K4.2$item),]
K4.3 <- K4.3[order(K4.3$item),]

# G R U P O  1: ································································

K4.1_3 <- inner_join(K3,K4.1,by="item")
colnames(K4.1_3) <- c("k","PertenenciaK3","consensus","PertenenciaK4.1")
Comprobacion_4.1_3 <- table(K4.1_3$PertenenciaK3)
Comprobacion_4.1_3
# Numero de grupo:        2
# Cantidad en cada grupo: 18

K4.1_5 <- inner_join(K5,K4.1,by="item")
colnames(K4.1_5) <- c("k","PertenenciaK5","consensus","PertenenciaK4.1")
Comprobacion_4.1_5 <- table(K4.1_5$PertenenciaK5)
Comprobacion_4.1_5
# Numero de grupo:         1 3 5       
# Cantidad en cada grupo:  2 8 8 

# G R U P O  2: ································································

K4.2_3 <- inner_join(K3,K4.2,by="item")
colnames(K4.2_3) <- c("k","PertenenciaK3","consensus","PertenenciaK4.2")
Comprobacion_4.2_3 <- table(K4.2_3$PertenenciaK3)
Comprobacion_4.2_3
# Numero de grupo:        1  2  3        
# Cantidad en cada grupo: 19  1  1 

K4.2_5 <- inner_join(K5,K4.2,by="item")
colnames(K4.2_5) <- c("k","PertenenciaK5","consensus","PertenenciaK4.2")
Comprobacion_4.2_5 <- table(K4.2_5$PertenenciaK5)
Comprobacion_4.2_5
# Numero de grupo:        1      
# Cantidad en cada grupo: 21

# G R U P O  3: ································································

K4.3_3 <- inner_join(K3,K4.3,by="item")
colnames(K4.3_3) <- c("k","PertenenciaK3","consensus","PertenenciaK4.3")
Comprobacion_4.3_3 <- table(K4.3_3$PertenenciaK3)
Comprobacion_4.3_3
# Numero de grupo:        1        
# Cantidad en cada grupo: 12

K4.3_5 <- inner_join(K5,K4.3,by="item")
colnames(K4.3_5) <- c("k","PertenenciaK5","consensus","PertenenciaK4.3")
Comprobacion_4.3_5 <- table(K4.3_5$PertenenciaK5)
Comprobacion_4.3_5
# Numero de grupo:        2      
# Cantidad en cada grupo: 12

# G R U P O  4: ································································

K4.4_3 <- inner_join(K3,K4.4,by="item")
colnames(K4.4_3) <- c("k","PertenenciaK3","consensus","PertenenciaK4.4")
Comprobacion_4.4_3 <- table(K4.4_3$PertenenciaK3)
Comprobacion_4.4_3
# Numero de grupo:        3        
# Cantidad en cada grupo: 10

K4.4_5 <- inner_join(K5,K4.4,by="item")
colnames(K4.4_5) <- c("k","PertenenciaK5","consensus","PertenenciaK4.4")
Comprobacion_4.4_5 <- table(K4.4_5$PertenenciaK5)
Comprobacion_4.4_5
# Numero de grupo:        4 5          
# Cantidad en cada grupo: 9 1 
 # K4
###################################################################################
# ··························· K = 5 ···········································

K5.1 <- K5 %>% filter(cluster==1) # 18 pacientes
K5.2 <- K5 %>% filter(cluster==2) # 21 pacientes
K5.3 <- K5 %>% filter(cluster==3) # 12 pacientes
K5.4 <- K5 %>% filter(cluster==4) # 10 pacientes
K5.5 <- K5 %>% filter(cluster==5) # 10 pacientes

K5.1 <- K5.1[,c(2,3)]
K5.2 <- K5.2[,c(2,3)]
K5.3 <- K5.3[,c(2,3)]

K5.1 <- K5.1[order(K5.1$item),]
K5.2 <- K5.2[order(K5.2$item),]
K5.3 <- K5.3[order(K5.3$item),]

# G R U P O  1: ································································

K5.1_3 <- inner_join(K3,K5.1,by="item")
colnames(K5.1_3) <- c("k","PertenenciaK3","consensus","PertenenciaK5.1")
Comprobacion_5.1_3 <- table(K5.1_3$PertenenciaK3)
Comprobacion_5.1_3
# Numero de grupo:        1  2  3         
# Cantidad en cada grupo: 19  3  1 


K5.1_4 <- inner_join(K4,K5.1,by="item")
colnames(K5.1_4) <- c("k","PertenenciaK4","consensus","PertenenciaK5.1")
Comprobacion_5.1_4 <- table(K5.1_4$PertenenciaK4)
Comprobacion_5.1_4
# Numero de grupo:         1  2      
# Cantidad en cada grupo:  2 21 

# G R U P O  2: ································································

K5.2_3 <- inner_join(K3,K5.2,by="item")
colnames(K5.2_3) <- c("k","PertenenciaK3","consensus","PertenenciaK5.2")
Comprobacion_5.2_3 <- table(K5.2_3$PertenenciaK3)
Comprobacion_5.2_3
# Numero de grupo:        1          
# Cantidad en cada grupo: 12


K5.2_4 <- inner_join(K4,K5.2,by="item")
colnames(K5.2_4) <- c("k","PertenenciaK4","consensus","PertenenciaK5.2")
Comprobacion_5.2_4 <- table(K5.2_4$PertenenciaK4)
Comprobacion_5.2_4
# Numero de grupo:         3     
# Cantidad en cada grupo:  12

# G R U P O  3: ································································

K5.3_3 <- inner_join(K3,K5.3,by="item")
colnames(K5.3_3) <- c("k","PertenenciaK3","consensus","PertenenciaK5.3")
Comprobacion_5.3_3 <- table(K5.3_3$PertenenciaK3)
Comprobacion_5.3_3
# Numero de grupo:        2          
# Cantidad en cada grupo: 8


K5.3_4 <- inner_join(K4,K5.3,by="item")
colnames(K5.3_4) <- c("k","PertenenciaK4","consensus","PertenenciaK5.3")
Comprobacion_5.3_4 <- table(K5.3_4$PertenenciaK4)
Comprobacion_5.3_4
# Numero de grupo:         1    
# Cantidad en cada grupo:  8

# G R U P O  4: ································································

K5.4_3 <- inner_join(K3,K5.4,by="item")
colnames(K5.4_3) <- c("k","PertenenciaK3","consensus","PertenenciaK5.4")
Comprobacion_5.4_3 <- table(K5.4_3$PertenenciaK3)
Comprobacion_5.4_3
# Numero de grupo:        3          
# Cantidad en cada grupo: 9


K5.4_4 <- inner_join(K4,K5.4,by="item")
colnames(K5.4_4) <- c("k","PertenenciaK4","consensus","PertenenciaK5.4")
Comprobacion_5.4_4 <- table(K5.4_4$PertenenciaK4)
Comprobacion_5.4_4
# Numero de grupo:         4    
# Cantidad en cada grupo:  9

# G R U P O  5: ································································

K5.5_3 <- inner_join(K3,K5.5,by="item")
colnames(K5.5_3) <- c("k","PertenenciaK3","consensus","PertenenciaK5.5")
Comprobacion_5.5_3 <- table(K5.5_3$PertenenciaK3)
Comprobacion_5.5_3
# Numero de grupo:        3          
# Cantidad en cada grupo: 9


K5.5_4 <- inner_join(K4,K5.5,by="item")
colnames(K5.5_4) <- c("k","PertenenciaK4","consensus","PertenenciaK5.5")
Comprobacion_5.5_4 <- table(K5.5_4$PertenenciaK4)
Comprobacion_5.5_4
# Numero de grupo:         1  4    
# Cantidad en cada grupo:  8  1
 # K5

#ExtracionGrupo que se mantiene: 

GrupoConstante <- K5 %>% filter(cluster==2)
saveRDS(GrupoConstante,"GrupoConstante_Luad.rds")


