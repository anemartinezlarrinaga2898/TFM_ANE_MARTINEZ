# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-07-2021

#OBJETIVO:

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")

################################################################################################################################
Huella1 <- readRDS("GLF1_Luad.rds")
Huella2 <- readRDS("GLF2_Luad.rds")
Huella3 <- readRDS("GLF3_Luad.rds")

Huella1 <- data.frame(Valores=Huella1)
Huella1$ensembl_gene_id <- rownames(Huella1)
rownames(Huella1) <- seq(1,nrow(Huella1),1)

Huella2 <- data.frame(Valores=Huella2)
Huella2$ensembl_gene_id <- rownames(Huella2)
rownames(Huella2) <- seq(1,nrow(Huella2),1)

Huella3 <- data.frame(Valores=Huella3)
Huella3$ensembl_gene_id <- rownames(Huella3)
rownames(Huella3) <- seq(1,nrow(Huella3),1)

AnotacionesHuella1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella1$ensembl_gene_id,
                            mart = ensembl)
Huella1 <- inner_join(Huella1,AnotacionesHuella1,by="ensembl_gene_id")

AnotacionesHuella2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella2$ensembl_gene_id,
                            mart = ensembl)
Huella2 <- inner_join(Huella2,AnotacionesHuella2,by="ensembl_gene_id")

AnotacionesHuella3 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella3$ensembl_gene_id,
                            mart = ensembl)
Huella3 <- inner_join(Huella3,AnotacionesHuella2,by="ensembl_gene_id")

saveRDS(Huella1,"Huella_1_LUAD_Anotada.rds")
saveRDS(Huella2,"Huella_2_LUAD_Anotada.rds")
saveRDS(Huella3,"Huella_3_LUAD_Anotada.rds")

Huella1 <- readRDS("GLF1_Lusc.rds")
Huella2 <- readRDS("GLF2_Lusc.rds")
Huella3 <- readRDS("GLF3_Lusc.rds")

Huella1 <- data.frame(Valores=Huella1)
Huella1$ensembl_gene_id <- rownames(Huella1)
rownames(Huella1) <- seq(1,nrow(Huella1),1)

Huella2 <- data.frame(Valores=Huella2)
Huella2$ensembl_gene_id <- rownames(Huella2)
rownames(Huella2) <- seq(1,nrow(Huella2),1)

Huella3 <- data.frame(Valores=Huella3)
Huella3$ensembl_gene_id <- rownames(Huella3)
rownames(Huella3) <- seq(1,nrow(Huella3),1)

AnotacionesHuella1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella1$ensembl_gene_id,
                            mart = ensembl)
Huella1 <- inner_join(Huella1,AnotacionesHuella1,by="ensembl_gene_id")

AnotacionesHuella2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella2$ensembl_gene_id,
                            mart = ensembl)
Huella2 <- inner_join(Huella2,AnotacionesHuella2,by="ensembl_gene_id")

AnotacionesHuella3 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters="ensembl_gene_id",
                            values=Huella3$ensembl_gene_id,
                            mart = ensembl)
Huella3 <- inner_join(Huella3,AnotacionesHuella2,by="ensembl_gene_id")

saveRDS(Huella1,"Huella_1_LUSC_Anotada.rds")
saveRDS(Huella2,"Huella_2_LUSC_Anotada.rds")
saveRDS(Huella3,"Huella_3_LUSC_Anotada.rds")

