# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 06-07-2021

#OBJETIVO:anotar mis huellas de las celulas del sistema inmune para luego afinar la busqueda. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
############################################################################################################################################################################################################################################################################################

HS <- readRDS("HuellasSistemaInmune_DF.rds")

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

AnotacionesHuella <- getBM(attributes = c("hgnc_symbol","description"),
                            filters="hgnc_symbol",
                            values=HS$GenesInmune,
                            mart = ensembl)

colnames(HS)[1] <- "hgnc_symbol"

HuellasAnotadas <- inner_join(HS,AnotacionesHuella,by="hgnc_symbol")
HuellasAnotadas <- HuellasAnotadas[!duplicated(HuellasAnotadas$hgnc_symbol),]

saveRDS(HuellasAnotadas,"HuellasSistemaInmune_DF.rds")

NombreHuellas <- HuellasAnotadas$hgnc_symbol


X <- NombreHuellas[str_which(NombreHuellas,pattern = "^IG[HK]")]
Y <- NombreHuellas[str_which(NombreHuellas,pattern = "^CCL")]
