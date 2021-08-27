# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 12-05-2021

#OBJETIVO:Analisis del ENRICHMENT ANALYSIS donde vamos ha hacer varias cosas: 
# Todo se hara con el BIOLOGICAL PROCESS
# 1) Anaslisis gsea con diferentes ontologies. 
#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
#########################################################################################################

H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)

C2 <-  msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene)

C5 <-  msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene)

C6 <-  msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, ensembl_gene)

C7 <-  msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, ensembl_gene)

# ···································· L U S C················································

# COMPARACION 1: 

geneList_Filtrado <- readRDS("GeneList/GeneList_Comparacion_1_Lusc_Filtrado.rds")
geneList_Total <- readRDS("GeneList/GeneList_Comparacion_1_Lusc_Total.rds")

geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

saveRDS(ResultadosGsea_Filtrado,"AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_1_Filtrado_Hallamarks.rds")
saveRDS(ResultadosGsea_Total,"AnalisisGSEA_Halmarks/GenesTotales/AnalisisGSEA_Lusc_1_Total_Hallamarks.rds")

# COMPARACION 2: 

geneList_Filtrado <- readRDS("GeneList/GeneList_Comparacion_2_Lusc_Filtrado.rds")
geneList_Total <- readRDS("GeneList/GeneList_Comparacion_2_Lusc_Total.rds")

geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

saveRDS(ResultadosGsea_Filtrado,"AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_2_Filtrado_Hallamarks.rds")
saveRDS(ResultadosGsea_Total,"AnalisisGSEA_Halmarks/GenesTotales/AnalisisGSEA_Lusc_2_Total_Hallamarks.rds")

# COMPARACION 3: 

geneList_Filtrado <- readRDS("GeneList/GeneList_Comparacion_3_Lusc_Filtrado.rds")
geneList_Total <- readRDS("GeneList/GeneList_Comparacion_3_Lusc_Total.rds")

geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

saveRDS(ResultadosGsea_Filtrado,"AnalisisGSEA_Halmarks/GenesFitlrados/AnalisisGSEA_Lusc_3_Filtrado_Hallamarks.rds")
saveRDS(ResultadosGsea_Total,"AnalisisGSEA_Halmarks/GenesTotales/AnalisisGSEA_Lusc_3_Total_Hallamarks.rds")



