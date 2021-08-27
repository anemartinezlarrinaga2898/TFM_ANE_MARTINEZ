# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 28-06-2021

#OBJETIVO:# COMPARACION 9: Grupo 3 LUAD GSEA frente al grupo 3 LUSC GSEA. Ambos son grupos POCO inmunogenicos. 
# a) si hay genes que tienes perfiles de expresion opuestos , cuantos genes se repiten y los perfiles de expresion 
# b) cuales son los genes que son unicos con sus perfiles de expresion 
# c) realizar las anotaciones de los genes. 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(VennDiagram)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/15_AnalisisEnrichment/")

################################################################################################################################

TG1 <- readRDS("DatosPartida/GSEA/AnalisisGSEA_Luad_3_Total.rds")
TG2 <- readRDS("DatosPartida/GSEA/AnalisisGSEA_Lusc_3_Total.rds")

#······································································································
# go term: adaptative inmune response --> GO:0002250

AIR_GeneSet <- data.frame(Genes=TG1@geneSets$`GO:0002250`)

#······································································································
GeneList_TG1 <- readRDS("DatosPartida/GSEA/GeneList/GeneList_Comparacion_3_Luad_Total.rds")
GeneList_TG2 <- readRDS("DatosPartida/GSEA/GeneList/GeneList_Comparacion_3_Lusc_Total.rds")

GeneList_TG1 <- data.frame(ValoresExpresion=GeneList_TG1,
                           Genes=names(GeneList_TG1))
GeneList_TG2 <- data.frame(ValoresExpresion=GeneList_TG2,
                           Genes=names(GeneList_TG2))

Expresion_TG1 <- inner_join(GeneList_TG1,AIR_GeneSet, by="Genes")
Expresion_TG2 <- inner_join(GeneList_TG2,AIR_GeneSet, by="Genes")
rownames(Expresion_TG1) <- Expresion_TG1$Genes
rownames(Expresion_TG2) <- Expresion_TG2$Genes
#······································································································

#······································································································
GenesTG1 <- Expresion_TG1$Genes
GenesTG2 <- Expresion_TG2$Genes

#Aqui dibujamos el diagrama de Venn: 
venn.diagram(
  x = list(GenesTG1, GenesTG2),
  category.names = c("Grupo 3 LUAD" , "Grupo 3 LUSC "),
  filename = "Comparacion_Grupo_3_LUAD_frente_a_Grupo_3_LUSC_GSEATOTAL.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#76D7C4", '#F1948A'),
  fill = c(alpha("#76D7C4",0.5), alpha('#F1948A',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

#······································································································

OverLap_TG1_TG2 <- calculate.overlap(x=list("Grupo1"=GenesTG1,
                                            "Grupo 2"=GenesTG2))
a3 <- OverLap_TG1_TG2$a3
a2 <- OverLap_TG1_TG2$a2  # GRUPO 2 LUSC
a1 <- OverLap_TG1_TG2$a1  # GRUPO 1 LUAD

i <- 1

for(i in seq_along(a3)){
  Gen <- a3[i]
  
  IndiceGen_Grupo3<- which(a2==Gen)
  IndiceGen_Grupo1 <- which(a1==Gen)
  
  a2 <- a2[-IndiceGen_Grupo3]
  a1 <- a1[-IndiceGen_Grupo1]
}

Expresion_GenesUnicos_TG1 <- Expresion_TG1[a1,]
Expresion_GenesUnicos_TG2 <- Expresion_TG2[a2,]

#······································································································

Expresion_TG1 <- Expresion_TG1 %>% mutate(Grupo=rep("Grupo 1",nrow(Expresion_TG1)))
Expresion_TG2 <- Expresion_TG2 %>% mutate(Grupo=rep("Grupo 2",nrow(Expresion_TG2)))
ExpresionTotal <- rbind(Expresion_TG1,Expresion_TG2)

Duplicados_TG1 <- Expresion_TG1[a3,] 
Duplicados_TG2 <- Expresion_TG2[a3,]

DuplicadosTotales <- rbind(Duplicados_TG1,Duplicados_TG2)

#······································································································

AnotacionesUnicos_TG1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = Expresion_GenesUnicos_TG1$Genes,
                               mart = ensembl)
AnotacionesUnicos_TG2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
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

saveRDS(Anotaciones_NivelesExpresion_TG1,"Comparacion9/Anotaciones_NivelesExpresion_TG1.rds")
saveRDS(Anotaciones_NivelesExpresion_TG2,"Comparacion9/Anotaciones_NivelesExpresion_TG2.rds")
saveRDS(Anotaciones_NivelesExpresion_Duplicados,"Comparacion9/Anotaciones_NivelesExpresion_Duplicados.rds")



