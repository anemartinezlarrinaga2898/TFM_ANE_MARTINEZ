# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 26-07-2021

#OBJETIVO:Generacion de los txt para el reposicionamiento

###############################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/20_ReposicionamientodeDrogas/")
library(tidyverse)
library(biomaRt)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

###############################################################################################################################################

#············ GRUPO 1 DE LUAD ··················

R1 <- readRDS("ResultadosAnalasisExpresion_1_Luad.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                                     row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]


R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_DOWN[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R1_DF_UP_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R1_DF_DOWN_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)

#············ GRUPO 2 DE LUAD ··················

R1 <- readRDS("ResultadosAnalasisExpresion_2_Luad.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                       row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]


R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_UP[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R2_DF_UP_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R2_DF_DOWN_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)

#············ GRUPO 3 DE LUAD ··················

R1 <- readRDS("ResultadosAnalasisExpresion_3_Luad.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                       row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]


R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_UP[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R3_DF_UP_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R3_DF_DOWN_LUAD.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)

###############################################################################################################################################

#············ GRUPO 1 DE LUSC ··················

R1 <- readRDS("ResultadosAnalasisExpresion_1_Lusc.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                       row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]

R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_UP[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R1_DF_UP_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R1_DF_DOWN_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)

#············ GRUPO 2 DE LUSC ··················

R1 <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                       row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]

R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_UP[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R2_DF_UP_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R2_DF_DOWN_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)

#············ GRUPO 3 DE LUSC ··················

R1 <- readRDS("ResultadosAnalasisExpresion_3_Lusc.rds")
R1_DF <- as.data.frame(dplyr::mutate(as.data.frame(R1)),
                       row.names=rownames(R1))
R1_DF$ensembl_gene_id <- rownames(R1_DF)

Anotaciones <- Anotaciones <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                                    filters = "ensembl_gene_id",
                                    values = R1_DF$ensembl_gene_id,
                                    mart=ensembl)
Anotaciones <- na.omit(Anotaciones)

R1_DF <- inner_join(R1_DF,Anotaciones,by="ensembl_gene_id")
R1_DF <- R1_DF[,-7]

R1_DF <- R1_DF %>% filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))
R1_DF <- na.omit(R1_DF)

R1_DF <- R1_DF[!duplicated(R1_DF$entrezgene_id),]

R1_DF <- R1_DF[order(R1_DF$log2FoldChange,decreasing = TRUE),]
R1_DF_UP <- R1_DF %>% filter(log2FoldChange >1.5)
R1_DF_DOWN <- R1_DF %>% filter(log2FoldChange < -1.5)

R1_DF_UP <- R1_DF_UP[1:150,]
R1_DF_DOWN <- R1_DF_DOWN[1:150,]

R1_DF_UP <- matrix(R1_DF_UP$entrezgene_id,nrow=nrow(R1_DF_UP),ncol=1)
R1_DF_DOWN <- matrix(R1_DF_DOWN$entrezgene_id,nrow=nrow(R1_DF_DOWN),ncol=1)

write.table(R1_DF_UP,"R3_DF_UP_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)
write.table(R1_DF_DOWN,"R3_DF_DOWN_LUSC.txt",sep="\t",row.names = FALSE,quote = FALSE, col.names = FALSE)


