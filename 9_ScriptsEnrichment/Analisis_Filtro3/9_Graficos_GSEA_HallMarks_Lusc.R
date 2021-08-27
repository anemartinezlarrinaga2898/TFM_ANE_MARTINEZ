#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-06-2021

#OBJETIVO:Analisis del enrichment de LUAD y LUSC con los diferentes hallmarks con el filtro 1.5

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
#················································ L U S C ·····························································

# C O M P A R A C I O N   1: ··························································································

# Cargamos los valores: 
RF <- readRDS("Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_1_Filtrado_Hallamarks.rds")
RT <- RF


Condiciones <- c("H","C2","C5","C6","C7")

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

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_1.rds")

H <- RT[[1]]
C2 <- RT[[2]]
C5 <- RT[[3]]
C6 <- RT[[4]]
C7 <- RT[[5]]

H <- ggplot(H, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="H Genes , Group 1 Lusc, Filtrados")
H

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 1 Lusc, Filtrados")
GC2

GC6 <- ggplot(C6, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C6 Genes , Group 1 Lusc")
GC6

GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 1 Lusc")
GC5

GC7 <- ggplot(C7, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 1 Lusc")
GC7

# C O M P A R A C I O N   2: ··························································································

# Cargamos los valores: 
RF <- readRDS("Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_2_Filtrado_Hallamarks.rds")
RT <- RF

RT <- RT[c(-1,-4)]


Condiciones <- c("C2","C5","C7")

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

saveRDS(ID_Cuentas,"Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_2.rds")

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_2.rds")

C2 <- RT[[1]]
C5 <- RT[[2]]
C7 <- RT[[3]]

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 2 Lusc, Filtrados")
GC2


GC5 <- ggplot(C5[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 2 Lusc")
GC5

GC7 <- ggplot(C7, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 2 Lusc")
GC7

# C O M P A R A C I O N   3: ··························································································

# Cargamos los valores: 
RF <- readRDS("Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_3_Filtrado_Hallamarks.rds")
RT <- RF

RT <- RT[-4]

Condiciones <- c("H","C2","C5","C7")

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

RT <- readRDS("Graficos_GSEA_Halmarks/Resultados_Hallmarks_Filtrados_Lusc_3.rds")

H <- RT[[1]]
C2 <- RT[[2]]
C5 <- RT[[3]]
C7 <- RT[[4]]


GH <- ggplot(H, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="H Genes , Group 3 Lusc, Filtrados")
GH

GC2 <- ggplot(C2[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C2 Genes , Group 3 Lusc, Filtrados")
GC2


GC5 <- ggplot(C5, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C5 Genes , Group 3 Lusc")
GC5


GC7 <- ggplot(C7, aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Cuentas))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="C7 Genes , Group 3 Lusc")
GC7


