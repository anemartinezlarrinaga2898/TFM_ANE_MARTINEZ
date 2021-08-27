# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
# Generacion de graficos para el analisis de GSEA 

# GRAFICOS GSEA LUAD H,C2,C5,C6,C7

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)
library(ggpubr)
##############################################################################################################

Condiciones <- c("H","C2","C5","C6","C7")

# ···································· L U S C················································
# C O M P A R A C I O N   1: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_1_Filtrado_Hallamarks.rds")
RT <- RF

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}


i <- 1

Lista <- list()

for ( i in seq_along(RT)){
  R <- RT[[i]]
  RD <- R@result
  
  RCuentas <- data.frame(ID=rownames(RD),Cuentas=rep(0,nrow(RD)))
  
  j <- 1
  
  for ( j in seq(1,nrow(RD),1)){
    CoreEnrichment <- RD$core_enrichment[j]
    
    Cuentas <- ObtencionCuentas(CoreEnrichment)
    
    RCuentas[j,2] <- Cuentas
    
  }
  
  Lista[[i]] <- RCuentas
  
}

names(Lista) <- Condiciones
ID_Cuentas <- list()

i <- 1
for ( i in seq_along(Lista)){
  X <- Lista[[i]]
  Y <- RT[[i]]
  Y <- Y@result
  
  Z <- inner_join(X,Y,by="ID")
  
  ID_Cuentas[[i]] <- Z
  
}

names(ID_Cuentas) <- Condiciones

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_1.rds")

# GRAFICAMOS: 

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_1.rds")

H <- RT[[1]]
C2 <- RT[[2]]
C5 <- RT[[3]]
C6 <- RT[[4]]
C7 <- RT[[5]]

GH <- ggplot(H, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Hallmark Genes , Group 1 Lusc")
GH


GC2 <- ggplot(C2, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 1 Lusc")
GC2

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 1 Lusc")
GC5

GC6 <- ggplot(C6, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C6 Genes , Group 1 Lusc")
GC6

GC7 <- ggplot(C7, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 1 Lusc")
GC7

# C O M P A R A C I O N   2: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_2_Filtrado_Hallamarks.rds")
RT <- RF[[3]]

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}


i <- 1

Lista <- list()

for ( i in seq_along(RT)){
  R <- RT
  RD <- R@result
  
  RCuentas <- data.frame(ID=rownames(RD),Cuentas=rep(0,nrow(RD)))
  
  j <- 1
  
  for ( j in seq(1,nrow(RD),1)){
    CoreEnrichment <- RD$core_enrichment[j]
    
    Cuentas <- ObtencionCuentas(CoreEnrichment)
    
    RCuentas[j,2] <- Cuentas
    
  }
  
  Lista[[i]] <- RCuentas
  
}

names(Lista) <- Condiciones
ID_Cuentas <- list()

X <- Lista[[1]]
Y <- RT@result
Z <- inner_join(X,Y,by="ID")
ID_Cuentas <- Z
  
}

names(ID_Cuentas) <- Condiciones

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_2.rds")

# GRAFICAMOS: 

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_2.rds")

C5 <- RT

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 2 Lusc")
GC5

# C O M P A R A C I O N   3: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_3_Filtrado_Hallamarks.rds")
RT_1 <- RT[[2]]
RT_2 <- RT[[5]]

RF <- list(RT_1,RT_2)

RT <- RF

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}


i <- 1

Lista <- list()

for ( i in seq_along(RT)){
  R <- RT[[i]]
  RD <- R@result
  
  RCuentas <- data.frame(ID=rownames(RD),Cuentas=rep(0,nrow(RD)))
  
  j <- 1
  
  for ( j in seq(1,nrow(RD),1)){
    CoreEnrichment <- RD$core_enrichment[j]
    
    Cuentas <- ObtencionCuentas(CoreEnrichment)
    
    RCuentas[j,2] <- Cuentas
    
  }
  
  Lista[[i]] <- RCuentas
  
}

names(Lista) <- Condiciones
ID_Cuentas <- list()

i <- 1
for ( i in seq_along(Lista)){
  X <- Lista[[i]]
  Y <- RT[[i]]
  Y <- Y@result
  
  Z <- inner_join(X,Y,by="ID")
  
  ID_Cuentas[[i]] <- Z
  
}

names(ID_Cuentas) <- Condiciones

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_3.rds")

# GRAFICAMOS: 

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_3.rds")

C2 <- RT[[1]]
C7 <- RT[[2]]


GC2 <- ggplot(C2, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 3 Lusc")
GC2

GC7 <- ggplot(C7, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 3 Lusc")
GC7