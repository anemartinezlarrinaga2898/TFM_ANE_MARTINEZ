# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-05-2021

#OBJETIVO:Analisis de enrichment con graficos para los genes filtrados, aunque no tengo de todos. 

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)
library(ggpubr)

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

##############################################################################################################

Condiciones <- c("H","C2","C5","C6","C7")

# ···································· L U A D················································
# C O M P A R A C I O N   1: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Luad_1_Filtrado_Hallamarks.rds")
RF_1 <- RF[[2]]
RF_2 <- RF[[3]]
RF_3 <- RF[[5]]

Condiciones <- c("C2","C5","C7")

RF <- list(RF_1,RF_2,RF_3)
RT <- RF

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

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_1.rds")

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_1.rds")

C2 <- RT[[1]]
C5 <- RT[[2]]
C7 <- RT[[3]]

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 1 Luad, Filtrados")
GC2

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 1 Luad")
GC5

GC7 <- ggplot(C7[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 1 Luad")
GC7


# C O M P A R A C I O N   2: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Luad_2_Filtrado_Hallamarks.rds")
RF_1 <- RF[[2]]
RF_2 <- RF[[3]]

Condiciones <- c("C2","C5")

RF <- list(RF_1,RF_2)
RT <- RF

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

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_2.rds")

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_2.rds")

C2 <- RT[[1]]
C5 <- RT[[2]]

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 2 Luad, Filtrados")
GC2

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 2 Luad")
GC5

# C O M P A R A C I O N   3: ·································································

RF <- readRDS("AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Luad_3_Filtrado_Hallamarks.rds")
RF_1 <- RF[[2]]
RF_2 <- RF[[3]]

Condiciones <- c("C2","C5")

RF <- list(RF_1,RF_2)
RT <- RF

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

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_3.rds")

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Luad_3.rds")

C2 <- RT[[1]]
C5 <- RT[[2]]

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 3 Luad, Filtrados")
GC2

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 3 Luad")
GC5




