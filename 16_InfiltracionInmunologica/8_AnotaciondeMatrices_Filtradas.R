# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 15-07-2021

#OBJETIVO:Anotacion de los genes de mis matrices con el hugo symbol de las MATRICES FILTRADAS

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/19_AnalisisMatricesOriginales/")

################################################################################################################################
Luad <- readRDS("Luad_Filtrado.rds")
rownames(Luad) <- gsub("\\..*","",rownames(Luad))
Lusc <- readRDS("LuSC_Filtrado.rds")
rownames(Lusc) <- gsub("\\..*","",rownames(Lusc))

GenesLuad <- data.frame(ensembl_gene_id=rownames(Luad))
GenesLusc <- data.frame(ensembl_gene_id=rownames(Lusc))

AnotacionesLuad <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = GenesLuad$ensembl_gene_id,
                         mart=ensembl)

AnotacionesLusc <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = GenesLusc$ensembl_gene_id,
                         mart=ensembl)

UnionAnotacionesLuad <- inner_join(GenesLuad,AnotacionesLuad,by="ensembl_gene_id")
UnionAnotacionesLusc <- inner_join(GenesLusc,AnotacionesLusc,by="ensembl_gene_id")

UnionAnotacionesLuad <- UnionAnotacionesLuad[!duplicated(UnionAnotacionesLuad$ensembl_gene_id),]
UnionAnotacionesLuad <- UnionAnotacionesLuad[!duplicated(UnionAnotacionesLuad$hgnc_symbol),]

UnionAnotacionesLusc <- UnionAnotacionesLusc[!duplicated(UnionAnotacionesLusc$ensembl_gene_id),]
UnionAnotacionesLusc <- UnionAnotacionesLusc[!duplicated(UnionAnotacionesLusc$hgnc_symbol),]


IndicesEspaciosEnBlancoLUAD <- str_which(UnionAnotacionesLuad$hgnc_symbol,pattern = "^\\s*$")
IndicesEspaciosEnBlancoLUSC <- str_which(UnionAnotacionesLusc$hgnc_symbol,pattern = "^\\s*$")

UnionAnotacionesLuad <- UnionAnotacionesLuad[-IndicesEspaciosEnBlancoLUAD,]
UnionAnotacionesLusc <- UnionAnotacionesLusc[-IndicesEspaciosEnBlancoLUSC,]

MatrizAnotadaLuad <- Luad[UnionAnotacionesLuad$ensembl_gene_id,]
rownames(MatrizAnotadaLuad) <- UnionAnotacionesLuad$hgnc_symbol

MatrizAnotadaLusc <- Lusc[UnionAnotacionesLusc$ensembl_gene_id,]
rownames(MatrizAnotadaLusc) <- UnionAnotacionesLusc$hgnc_symbol

saveRDS(MatrizAnotadaLuad,"MatrizFiltrada_LUAD_ANOTADA.rds")
saveRDS(MatrizAnotadaLusc,"MatrizFiltrada_LUSC_ANOTADA.rds")
