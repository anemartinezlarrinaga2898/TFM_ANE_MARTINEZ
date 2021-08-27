# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
# Generacion de graficos para el analisis de GSEA 

# GRAFICOS GSEA LUAD BP

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)
library(ggpubr)
##############################################################################################################

# ···································· L U A D················································

# B I O L O G I C A L  P R O C E S S : 

# Graficos: 
  #Filtrado: 
    # Mi dotplot representado el NES y el tamaño el numero de genes que aparecere en el core_enrichment 
    # El gsePLOT de las primeras 5 
  # Total 
    # Mi Dotplot de las primeras 15 geneSet IDs 
    # Un emmaplot de todos los valores para ver si hay funciones que se relacionan entre ellas( tendre que ver si con todos o con los principales)
    # Cnetplot para ver cuales de los genes forma parte de esos clusteres que saque 

# C O M P A R A C I O N   1: ··························································································

# Cargamos los valores: 

Filtrados <- readRDS("AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_1_Filtrados.rds")
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_1_Total.rds")

RFiltrados <- Filtrados@result
RTotales <- Totales@result

#Obtenemos el numero de genes que contribuyen a cada genset para luego establecerlo como size de nuestro grafico: 

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

R <- list(RFiltrados,RTotales)
i <- 1

CFiltrados <- data.frame(ID=rownames(RFiltrados),Cuentas=rep(0,nrow(RFiltrados)))
CTotales <- data.frame(ID=rownames(RTotales),Cuentas=rep(0,nrow(RTotales)))

for ( i in seq_along(R)){
  Res <- R[[i]]
  j <- 1
  for ( j in seq(1,nrow(Res),1)){
    CoreEnrichment <- Res$core_enrichment[j]
    
    Cuentas <- ObtencionCuentas(CoreEnrichment)
    
    if (i==1){CFiltrados[j,2] <- Cuentas}
    if(i==2){CTotales[j,2] <- Cuentas}
  }
  
  RFiltrados <- inner_join(RFiltrados,CFiltrados,by="ID")
  
  RTotales <- inner_join(RTotales,CTotales,by="ID")
 
}

RFiltrados <- RFiltrados[,-(12)]
colnames(RFiltrados)[12] <- "Count"

RTotales <- RTotales[,-(12)]
colnames(RTotales)[12] <- "Count"

saveRDS(RFiltrados,"Graficos_GSEA_BP/Resultados_Filtrados_BP_Luad_1.rds")
saveRDS(RTotales,"Graficos_GSEA_BP/Resultados_Totales_BP_Luad_1.rds")

# Ahora graficamos: 

RFiltrados <- readRDS("Graficos_GSEA_BP/Resultados_Filtrados_BP_Luad_1.rds")
RTotales <- readRDS("Graficos_GSEA_BP/Resultados_Totales_BP_Luad_1.rds")


GR_1F <- ggplot(RFiltrados[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Filtered Genes")
GR_1F

GR_1T <- ggplot(RTotales[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Total Genes")
GR_1T

figura <- ggarrange(GR_1F,GR_1T,
                    labels = c("A","B"), ncol = 2, nrow = 1)
annotate_figure(figura,
                top = text_grob("Biological Process Group 1 Luad", color = "Black", size = 14))


# C O M P A R A C I O N   2: ··························································································

# Cargamos los valores: 

Filtrados <- readRDS("AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_2_Filtrados.rds")
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_2_Total.rds")

#RFiltrados <- Filtrados@result
RTotales <- Totales@result

#Obtenemos el numero de genes que contribuyen a cada genset para luego establecerlo como size de nuestro grafico: 

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

#R <- list(RFiltrados,RTotales)
i <- 1

#CFiltrados <- data.frame(ID=rownames(RFiltrados),Cuentas=rep(0,nrow(RFiltrados)))
CTotales <- data.frame(ID=rownames(RTotales),Cuentas=rep(0,nrow(RTotales)))

Res <- RTotales
for ( j in seq(1,nrow(Res),1)){
  CoreEnrichment <- Res$core_enrichment[j]
  
  Cuentas <- ObtencionCuentas(CoreEnrichment)
  
  CTotales[j,2] <- Cuentas
}

RTotales <- inner_join(RTotales,CTotales,by="ID")
colnames(RTotales)[12] <- "Count"

saveRDS(RTotales,"Graficos_GSEA_BP/Resultados_Totales_BP_Luad_2.rds")

# Ahora graficamos: 

RTotales <- readRDS("Graficos_GSEA_BP/Resultados_Totales_BP_Luad_2.rds")

GR_1T <- ggplot(RTotales[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Biological Process, Group 2 ,Total Genes, Luad")
GR_1T


# C O M P A R A C I O N   3: ··························································································

# Cargamos los valores: 

Filtrados <- readRDS("AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_3_Filtrados.rds")
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_3_Total.rds")

#RFiltrados <- Filtrados@result
RTotales <- Totales@result

#Obtenemos el numero de genes que contribuyen a cada genset para luego establecerlo como size de nuestro grafico: 

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

#R <- list(RFiltrados,RTotales)
i <- 1

#CFiltrados <- data.frame(ID=rownames(RFiltrados),Cuentas=rep(0,nrow(RFiltrados)))
CTotales <- data.frame(ID=rownames(RTotales),Cuentas=rep(0,nrow(RTotales)))

Res <- RTotales
for ( j in seq(1,nrow(Res),1)){
    CoreEnrichment <- Res$core_enrichment[j]
    
    Cuentas <- ObtencionCuentas(CoreEnrichment)
    
    CTotales[j,2] <- Cuentas
  }

RTotales <- inner_join(RTotales,CTotales,by="ID")
colnames(RTotales)[12] <- "Count"

saveRDS(RTotales,"AnalisisGSEA_BP/Resultados_Totales_BP_Luad_3.rds")

# Ahora graficamos: 

RTotales <- readRDS("AnalisisGSEA_BP/Resultados_Totales_BP_Luad_3.rds")

GR_1T <- ggplot(RTotales[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Biological Process, Group 3 ,Total Genes, Luad")
GR_1T





