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

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/18_BusquedaGenes/")

################################################################################################################################

Primitivo <- c("MCM10","E2F","E2F3","TYMS","POLA1","E2F8")
Clasico <- c("GPX2","ALDH3A1","AKR1C3","GPX2","AKR1C1","TXNRD1","GSTM3","TP63")
Secretorio <- c("ARHGDIB","TNFRSF14","MUC1","SFTPC","SFTPB","SFTPD","NKX2","TTF1")
Basal <- c("LAMB3","LAMC2","COL11A1","COL17A1","ITGB4","ITGB5","CLDN1","KRT5","PSORIASIN","S100A7","GFB5","S100")

ParticipantesPrimitivo <- rep("Primitivo",length(Primitivo))
ParticipantesClasico <- rep("Clasico",length(Clasico))
ParticipantesSecretorio <- rep("Secretorio",length(Secretorio))
ParticipantesBasal <- rep("Basal",length(Basal))

GenesPaper <- data.frame(Genes=c(Primitivo,
                                 Clasico,
                                 Secretorio,
                                 Basal),
                         TipoCelular=c(ParticipantesPrimitivo, 
                                       ParticipantesClasico, 
                                       ParticipantesSecretorio,
                                       ParticipantesBasal))
GenesPaper <- GenesPaper[,c(2,1)]

saveRDS(GenesPaper,"GenesPaper_Lusc.rds")

#·······················································································································

Huella1 <- readRDS("GLF1_Lusc.rds")
Huella2 <- readRDS("GLF2_Lusc.rds")
Huella3 <- readRDS("GLF3_Lusc.rds")

Huella1 <- data.frame(ValoresExpresion=Huella1)
Huella1$ensembl_gene_id <- rownames(Huella1)
Huella1$Grupo <- rep("Grupo1",nrow(Huella1))

Huella2 <- data.frame(ValoresExpresion=Huella2)
Huella2$ensembl_gene_id <- rownames(Huella2)
Huella2$Grupo <- rep("Grupo2",nrow(Huella2))

Huella3 <- data.frame(ValoresExpresion=Huella3)
Huella3$ensembl_gene_id <- rownames(Huella3)
Huella3$Grupo <- rep("Grupo3",nrow(Huella3))

Huellas <- rbind(Huella1,Huella2)
Huellas <- rbind(Huellas,Huella3)
rownames(Huellas) <- seq(1,nrow(Huellas),1)
Huellas <- Huellas[,c(2,3,1)]

Anotaciones <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = Huellas$ensembl_gene_id,
                               mart = ensembl)

Huellas <- inner_join(Huellas,Anotaciones,by="ensembl_gene_id")
Huellas <- Huellas[,c(1,4,2,3)]

saveRDS(Huellas,"Huellas_Lusc_Anotadas.rds")

##############################################################################################################################

Huellas <- readRDS("Huellas_Lusc_Anotadas.rds")
GenesPaper <- readRDS("GenesPaper_Lusc.rds")
colnames(GenesPaper)[2] <- "hgnc_symbol"

Huellas1 <- Huellas %>% filter(Grupo=="Grupo1")
Huellas2 <- Huellas %>% filter(Grupo=="Grupo2")
Huellas3 <- Huellas %>% filter(Grupo=="Grupo3")

Union <- inner_join(GenesPaper,Huellas,by="hgnc_symbol")


Huella2_Union <- inner_join(GenesPaper,Huellas2,by="hgnc_symbol")


