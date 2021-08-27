#Fecha: 01-07-2021

#OBJETIVO:# COMPARACION 16: Comparacion de las tres huellas para sacar los genes que son unicos de cada grupo
# HUELLA 2
###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(VennDiagram)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/15_AnalisisEnrichment/")

################################################################################################################################

GeneList_TG1 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_1_Luad_Filtrado.rds")
GeneList_TG2 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_2_Luad_Filtrado.rds")
GeneList_TG3 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_3_Luad_Filtrado.rds")

GeneList_TG1 <- data.frame(ValoresExpresion=GeneList_TG1,
                           Genes=names(GeneList_TG1))
GeneList_TG2 <- data.frame(ValoresExpresion=GeneList_TG2,
                           Genes=names(GeneList_TG2))
GeneList_TG3 <- data.frame(ValoresExpresion=GeneList_TG3,
                           Genes=names(GeneList_TG3))

GenesTG1 <- GeneList_TG1$Genes
GenesTG2 <- GeneList_TG2$Genes
GenesTG3 <- GeneList_TG3$Genes

venn.diagram(
  x = list(GenesTG1, GenesTG2,GenesTG3),
  category.names = c("Group 1" , "Group 2","Group 3"),
  filename = "Huellas_Luad.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#6988E1", '#F3679C',"#B78FDA"),
  fill = c(alpha("#6988E1",0.5), alpha('#F3679C',0.5),alpha('#B78FDA',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

OverLap_TG1_TG2 <- calculate.overlap(x=list("Grupo1"=GenesTG1,
                                            "Grupo 2"=GenesTG2,
                                            "Grupo 3"=GenesTG3))

G1 <- OverLap_TG1_TG2[[5]]
G2 <- OverLap_TG1_TG2[[6]]
G3 <- OverLap_TG1_TG2[[7]]

AnotacionesUnicos_TG1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                               filters = "ensembl_gene_id",
                               values = G1,
                               mart = ensembl)
AnotacionesUnicos_TG2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                               filters = "ensembl_gene_id",
                               values = G2,
                               mart = ensembl)

AnotacionesUnicos_TG3 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                               filters = "ensembl_gene_id",
                               values = G3,
                               mart = ensembl)

colnames(GeneList_TG1)[2] <- "ensembl_gene_id"
colnames(GeneList_TG2)[2] <- "ensembl_gene_id"
colnames(GeneList_TG3)[2] <- "ensembl_gene_id"

G1 <- inner_join(AnotacionesUnicos_TG1,GeneList_TG1,by="ensembl_gene_id")
G2 <- inner_join(AnotacionesUnicos_TG2,GeneList_TG2,by="ensembl_gene_id")
G3 <- inner_join(AnotacionesUnicos_TG3,GeneList_TG3,by="ensembl_gene_id")

saveRDS(G1,"Comparacion16/G1.rds")
saveRDS(G2,"Comparacion16/G2.rds")
saveRDS(G3,"Comparacion16/G3.rds")


