# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 24-05-2021

#OBJETIVO:Script para obtener las huellas y anotar las mismas con diferentes anotaciones.  

#############################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/12_Huellas/")
library(biomaRt)
library(tidyverse)

############################################################################################################################################
#Variables basicas para establecer y poder realizar las anotaciones con biomaRt

ensembl <-  useMart("ensembl",dataset="hsapiens_gene_ensembl")
anotaciones <- listAttributes(ensembl)
filtros <- listFilters(ensembl)

############################################################################################################################################

# ······································· L U A D ··········································

GenesGrupo1 <- readRDS("ResultadosFiltrados_Luad_Comparacion1.rds")
GenesGrupo2 <- readRDS("ResultadosFiltrados_Luad_Comparacion2.rds")
GenesGrupo3 <- readRDS("ResultadosFiltrados_Luad_Comparacion3.rds")

Genes_1 <- rownames(GenesGrupo1)
Genes_2 <- rownames(GenesGrupo2)
Genes_3 <- rownames(GenesGrupo3)

Huella_1 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id","name_1006"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_1, 
                  mart = ensembl)

Huella_2 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_2, 
                  mart = ensembl)

Huella_3 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_3, 
                  mart = ensembl)

saveRDS(Genes_1,"Huella_Ensembel_Luad_1.rds")
saveRDS(Genes_2,"Huella_Ensembel_Luad_2.rds")
saveRDS(Genes_3,"Huella_Ensembel_Luad_3.rds")

saveRDS(Huella_1,"Huella_Ensembel_HUGO_EntrezGene_Luad_1.rds")
saveRDS(Huella_2,"Huella_Ensembel_HUGO_EntrezGene_Luad_2.rds")
saveRDS(Huella_3,"Huella_Ensembel_HUGO_EntrezGene_Luad_3.rds")

# ······································· L U S C ··········································

GenesGrupo1 <- readRDS("ResultadosFiltrados_Lusc_Comparacion1.rds")
GenesGrupo2 <- readRDS("ResultadosFiltrados_Lusc_Comparacion2.rds")
GenesGrupo3 <- readRDS("ResultadosFiltrados_Lusc_Comparacion3.rds")

Genes_1 <- rownames(GenesGrupo1)
Genes_2 <- rownames(GenesGrupo2)
Genes_3 <- rownames(GenesGrupo3)

Huella_1 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_1, 
                  mart = ensembl)

Huella_2 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_2, 
                  mart = ensembl)

Huella_3 <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), 
                  filters = 'ensembl_gene_id', 
                  values = Genes_3, 
                  mart = ensembl)

saveRDS(Genes_1,"Huella_Ensembel_Lusc_1.rds")
saveRDS(Genes_2,"Huella_Ensembel_Lusc_2.rds")
saveRDS(Genes_3,"Huella_Ensembel_Lusc_3.rds")

saveRDS(Huella_1,"Huella_Ensembel_HUGO_EntrezGene_Lusc_1.rds")
saveRDS(Huella_2,"Huella_Ensembel_HUGO_EntrezGene_Lusc_2.rds")
saveRDS(Huella_3,"Huella_Ensembel_HUGO_EntrezGene_Lusc_3.rds")

