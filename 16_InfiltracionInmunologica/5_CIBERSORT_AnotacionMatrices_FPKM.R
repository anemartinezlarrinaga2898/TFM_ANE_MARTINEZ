# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 15-07-2021

# OBJETIVO: Anotacion de las matrices de FPKM y luego las voy a organizar por pacientes. 

##########################################################################################################################################

directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filter <- listFilters(ensembl)
attributos <- listAttributes(ensembl)

##########################################################################################################################################

Matriz <- readRDS("MatrizLuad_FPKM.rds")
dim(Matriz)
Matriz <- Matriz[apply(Matriz[,-1], 1, function(x) !all(x==0)),]
dim(Matriz)

Genes <- data.frame(ensembl_gene_id=rownames(Matriz))

AnotacionesH <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters = "ensembl_gene_id",
                            values = Genes$ensembl_gene_id,
                            mart = ensembl)

UnionAnotaciones <- inner_join(Genes,AnotacionesH,by="ensembl_gene_id")
UnionAnotaciones <- UnionAnotaciones[!duplicated(UnionAnotaciones$ensembl_gene_id),]
UnionAnotaciones <- UnionAnotaciones[!duplicated(UnionAnotaciones$hgnc_symbol),]

IndicesEspaciosEnBlanco <- str_which(UnionAnotaciones$hgnc_symbol,pattern = "^\\s*$")

UnionAnotaciones <- UnionAnotaciones[-IndicesEspaciosEnBlanco,]

MatrizAnotada <- Matriz[UnionAnotaciones$ensembl_gene_id,]
rownames(MatrizAnotada) <- UnionAnotaciones$hgnc_symbol

K1 <- readRDS("K3.1_NMF_Luad.rds")
K2 <- readRDS("K3.2_NMF_Luad.rds")
K3 <- readRDS("K3.3_NMF_Luad.rds")

G1 <- MatrizAnotada[,K1$Sample_ID]
G2 <- MatrizAnotada[,K2$Sample_ID]
G3 <- MatrizAnotada[,K3$Sample_ID]

MatrizFinal <- cbind(G1,G2)
MatrizFinal <- cbind(MatrizFinal,G3)

saveRDS(MatrizFinal,"MatrizFinal_LUAD_Anotada.rds")

NombreGenes <- data.frame(GeneSymbol=rownames(MatrizFinal))
MatrizFinal <- cbind(NombreGenes,MatrizFinal)

write.table(MatrizFinal,"MatrizFinal_LUAD_ANOTADA.txt",quote = FALSE,sep = "\t",row.names = FALSE)

G1 <- cbind(NombreGenes,G1)
G2 <- cbind(NombreGenes,G2)
G3 <- cbind(NombreGenes,G3)

write.table(G1,"G1_LUAD.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(G2,"G2_LUAD.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(G3,"G3_LUAD.txt",quote = FALSE,sep = "\t",row.names = FALSE)


##################################################################################################################################

Matriz <- readRDS("MatrizLusc_FPKM.rds")
dim(Matriz)
Matriz <- Matriz[apply(Matriz[,-1], 1, function(x) !all(x==0)),]
dim(Matriz)

Genes <- data.frame(ensembl_gene_id=rownames(Matriz))

AnotacionesH <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = Genes$ensembl_gene_id,
                      mart = ensembl)

UnionAnotaciones <- inner_join(Genes,AnotacionesH,by="ensembl_gene_id")
UnionAnotaciones <- UnionAnotaciones[!duplicated(UnionAnotaciones$ensembl_gene_id),]
UnionAnotaciones <- UnionAnotaciones[!duplicated(UnionAnotaciones$hgnc_symbol),]

IndicesEspaciosEnBlanco <- str_which(UnionAnotaciones$hgnc_symbol,pattern = "^\\s*$")

UnionAnotaciones <- UnionAnotaciones[-IndicesEspaciosEnBlanco,]

MatrizAnotada <- Matriz[UnionAnotaciones$ensembl_gene_id,]
rownames(MatrizAnotada) <- UnionAnotaciones$hgnc_symbol

K1 <- readRDS("K3.1_NMF_Lusc.rds")
K2 <- readRDS("K3.2_NMF_Lusc.rds")
K3 <- readRDS("K3.3_NMF_Lusc.rds")

G1 <- MatrizAnotada[,K1$Sample_ID]
G2 <- MatrizAnotada[,K2$Sample_ID]
G3 <- MatrizAnotada[,K3$Sample_ID]

MatrizFinal <- cbind(G1,G2)
MatrizFinal <- cbind(MatrizFinal,G3)

saveRDS(MatrizFinal,"MatrizFinal_LUSC_Anotada.rds")

NombreGenes <- data.frame(GeneSymbol=rownames(MatrizFinal))
MatrizFinal <- cbind(NombreGenes,MatrizFinal)

write.table(MatrizFinal,"MatrizFinal_LUSC_ANOTADA.txt",quote = FALSE,sep = "\t",row.names = FALSE)

G1 <- cbind(NombreGenes,G1)
G2 <- cbind(NombreGenes,G2)
G3 <- cbind(NombreGenes,G3)

write.table(G1,"G1_LUSC.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(G2,"G2_LUSC.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(G3,"G3_LUSC.txt",quote = FALSE,sep = "\t",row.names = FALSE)


