#Fecha: 01-07-2021

#OBJETIVO:# COMPARACION 17: Comparacion de las tres huellas para sacar los genes que son unicos de cada grupo
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

GeneList_TG1 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_1_Lusc_Filtrado.rds")
GeneList_TG2 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_2_Lusc_Filtrado.rds")
GeneList_TG3 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_3_Lusc_Filtrado.rds")

GeneListLusc <- c(GeneList_TG1,GeneList_TG2,GeneList_TG3)

GeneList_TG4 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_1_Luad_Filtrado.rds")
GeneList_TG5 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_2_Luad_Filtrado.rds")
GeneList_TG6 <- readRDS("DatosPartida/GSEA/GeneListFiltro2/GeneList_Comparacion_3_Luad_Filtrado.rds")

GeneListLuad <- c(GeneList_TG4,GeneList_TG5,GeneList_TG6)


GeneList_Lusc <- data.frame(ValoresExpresion=GeneListLusc,
                           Genes=names(GeneListLusc))
GeneList_Lusc <- GeneList_Lusc[!duplicated(GeneList_Lusc$Genes),]
rownames(GeneList_Lusc) <- GeneList_Lusc$Genes
GeneList_Luad <- data.frame(ValoresExpresion=GeneListLuad,
                            Genes=names(GeneListLuad))
GeneList_Luad <- GeneList_Luad[!duplicated(GeneList_Luad$Genes),]
rownames(GeneList_Luad) <- GeneList_Luad$Genes


LUAD <- GeneList_Luad$Genes
LUSC <- GeneList_Lusc$Genes

venn.diagram(
  x = list(LUAD,LUSC),
  category.names = c("LUAD" , 
                     "LUSC"),
  filename = "Huellas_Lusc_Luad_II.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#A7DA8F", 
        '#F4DC55'),
  fill = c(alpha("#A7DA8F",0.5), 
           alpha('#F4DC55',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

OverLap_TG1_TG2 <- calculate.overlap(x=list("LUAD"=LUAD,
                                            "LUSC"=LUSC))

a1 <- OverLap_TG1_TG2$a1
a2 <- OverLap_TG1_TG2$a2
a3 <- OverLap_TG1_TG2$a3

i <- 1

for(i in seq_along(a3)){
  Gen <- a3[i]
  
  IndiceGen_Grupo3<- which(a2==Gen)
  IndiceGen_Grupo1 <- which(a1==Gen)
  
  a2 <- a2[-IndiceGen_Grupo3]
  a1 <- a1[-IndiceGen_Grupo1]
}

Expresion_GenesUnicos_TG1 <- GeneList_Luad[a1,]
Expresion_GenesUnicos_TG2 <- GeneList_Lusc[a2,]

Expresion_TG1 <- GeneList_TG1 %>% mutate(Grupo=rep("Grupo 1",nrow(GeneList_Luad)))
Expresion_TG2 <- GeneList_TG2 %>% mutate(Grupo=rep("Grupo 2",nrow(GeneList_Lusc)))
ExpresionTotal <- rbind(Expresion_TG1,Expresion_TG2)

Duplicados_TG1 <- GeneList_TG1[a3,] 
Duplicados_TG2 <- GeneList_TG2[a3,]

DuplicadosTotales <- rbind(Duplicados_TG1,Duplicados_TG2)

#······································································································

AnotacionesUnicos_TG1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                               filters = "ensembl_gene_id",
                               values = Expresion_GenesUnicos_TG1$Genes,
                               mart = ensembl)
AnotacionesUnicos_TG2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                               filters = "ensembl_gene_id",
                               values = Expresion_GenesUnicos_TG2$Genes,
                               mart = ensembl)
AnotacionesDuplicados <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = DuplicadosTotales$Genes,
                               mart = ensembl)

colnames(Expresion_GenesUnicos_TG1)[2] <- "ensembl_gene_id"
colnames(Expresion_GenesUnicos_TG2)[2] <- "ensembl_gene_id"
colnames(DuplicadosTotales)[2] <- "ensembl_gene_id"

Anotaciones_NivelesExpresion_TG1 <- inner_join(AnotacionesUnicos_TG1,Expresion_GenesUnicos_TG1,by="ensembl_gene_id")
Anotaciones_NivelesExpresion_TG2 <- inner_join(AnotacionesUnicos_TG2,Expresion_GenesUnicos_TG2,by="ensembl_gene_id")
Anotaciones_NivelesExpresion_Duplicados <- inner_join(AnotacionesDuplicados,DuplicadosTotales,by="ensembl_gene_id")

saveRDS(Anotaciones_NivelesExpresion_TG1,"Comparacion18/Anotaciones_NivelesExpresion_TG1.rds")
saveRDS(Anotaciones_NivelesExpresion_TG2,"Comparacion18/Anotaciones_NivelesExpresion_TG2.rds")
saveRDS(Anotaciones_NivelesExpresion_Duplicados,"Comparacion15/Anotaciones_NivelesExpresion_Duplicados.rds")







