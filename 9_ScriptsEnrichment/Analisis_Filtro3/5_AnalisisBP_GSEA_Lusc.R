# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-06-2021

#OBJETIVO:Analisis GSEA con BP el analisis y los graficos. 

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

##############################################################################################################
# ······································L U A D ····················································

# ·············  C O M P A R A C I O N  1 ···················

GLF <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF1_Lusc.rds")
GSEA_BP_F <- gseGO(geneList = GLF, 
                   OrgDb = org.Hs.eg.db, 
                   ont = 'BP',  
                   minGSSize = 10, 
                   keyType = "ENSEMBL",
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   eps=0) 
saveRDS(GSEA_BP_F,"Analisis_Filtro3/Resultados/GSEA/GSEA_Grupo1_Lusc_Filtrado.rds")

Filtrados <- GSEA_BP_F
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Lusc_1_Total.rds")

RFiltrados <- Filtrados@result
RTotales <- Totales@result

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

GR_1F <- ggplot(RFiltrados[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Filtered Genes")
GR_1F

# ·············  C O M P A R A C I O N  2 ···················

GLF <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF2_Lusc.rds")
GSEA_BP_F <- gseGO(geneList = GLF, 
                   OrgDb = org.Hs.eg.db, 
                   ont = 'BP',  
                   minGSSize = 10, 
                   keyType = "ENSEMBL",
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   eps=0) 

saveRDS(GSEA_BP_F,"Analisis_Filtro3/Resultados/GSEA/GSEA_Grupo2_Lusc_Filtrado.rds")

Filtrados <- GSEA_BP_F
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Lusc_2_Total.rds")

RFiltrados <- Filtrados@result
RTotales <- Totales@result

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

GR_1F <- ggplot(RFiltrados[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Filtered Genes")
GR_1F

# ·············  C O M P A R A C I O N  3 ···················

GLF <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF3_Lusc.rds")
GSEA_BP_F <- gseGO(geneList = GLF, 
                   OrgDb = org.Hs.eg.db, 
                   ont = 'BP',  
                   minGSSize = 10, 
                   keyType = "ENSEMBL",
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   eps=0) 
saveRDS(GSEA_BP_F,"Analisis_Filtro3/Resultados/GSEA/GSEA_Grupo3_Lusc_Filtrado.rds")

Filtrados <- GSEA_BP_F
Totales <- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Lusc_3_Total.rds")

RFiltrados <- Filtrados@result
RTotales <- Totales@result

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

GR_1F <- ggplot(RFiltrados[1:15,], aes(x=NES, y=Description))+
  geom_point(aes(colour=p.adjust,size=Count))+
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05))+
  labs(title="Filtered Genes")
GR_1F