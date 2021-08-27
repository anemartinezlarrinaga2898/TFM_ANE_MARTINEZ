# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 30-07-2021

#OBJETIVO:analisis over representative de mis comunities con mis huellas

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/18_BusquedaGenes/")

################################################################################################################################

UNIVERSE <- readRDS("Huellas_Lusc_Anotadas.rds")
TARGET <- readRDS("GenesPaper_Ampliados.rds")
colnames(TARGET)[1] <- "hgnc_symbol"

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

Anotaciones <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = TARGET$hgnc_symbol,
                     mart = ensembl)

TARGET <- inner_join(Anotaciones,TARGET,by="hgnc_symbol")
TARGET <- TARGET[!duplicated(TARGET$hgnc_symbol),]

UNIVERSE <- UNIVERSE[-which(UNIVERSE$hgnc_symbol== ""), ]

target <- TARGET$ensembl_gene_id
universe <- UNIVERSE$ensembl_gene_id

set.seed(123456789)
Enrichement_Target_Community_GO <- enrichGO(gene=target,
                                  universe = universe,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENSEMBL",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = TRUE,
                                  pAdjustMethod = "BH",
                                  ont="BP")

Pvalue <- Enrichement_Target_Community_GO@result$pvalue
Terminos <-Enrichement_Target_Community_GO@result$Description
Cuentas <- Enrichement_Target_Community_GO@result$Count

Resultados <- data.frame(Term=Terminos,Count=Cuentas,Pvalue=Pvalue)





