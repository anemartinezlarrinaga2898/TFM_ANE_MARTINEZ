# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 28-06-2021

#OBJETIVO:# COMPARACION 6: GRUPO 1 LUAD ( GSEA ) frente al GRUPO 3 LUAD ( GSEA ), GO TERM: Adaptative inmune system
#El grupo 1 de LUAD es MUY INMUNOGENICO mientras que el grupo 3 de LUAD es menos inmunogenicos vamos a ver: 
# a) si hay genes que tienes perfiles de expresion opuestos , cuantos genes se repiten y los perfiles de expresion 
# b) cuales son los genes que son unicos con sus perfiles de expresion 
# c) realizar las anotaciones de los genes. 
#Todo esto lo hacemos con el paquete VennDiagram 

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(VennDiagram)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/15_AnalisisEnrichment/")

################################################################################################################################

TG1 <- readRDS("DatosPartida/GSEA/AnalisisGSEA_Luad_1_Total.rds")
TG3 <- readRDS("DatosPartida/GSEA/AnalisisGSEA_Luad_3_Total.rds")

# go term: adaptative inmune response --> GO:0002250

AIR_GeneSet <- data.frame(Genes=TG1@geneSets$`GO:0002250`) #De aqui saco todos los genes que estan dentro de este GeneSet

GeneList_TG1 <- readRDS("DatosPartida/GSEA/GeneList/GeneList_Comparacion_1_Luad_Total.rds")
GeneList_TG3 <- readRDS("DatosPartida/GSEA/GeneList/GeneList_Comparacion_3_Luad_Total.rds")

GeneList_TG1 <- data.frame(ValoresExpresion=GeneList_TG1,
                           Genes=names(GeneList_TG1))
GeneList_TG3 <- data.frame(ValoresExpresion=GeneList_TG3,
                           Genes=names(GeneList_TG3))

#En este paso obtenemos los genes que estan dentro del GenSet de interes junto con sus valores de expresion: 
  # Gene list es lo que hemos usado para hacer el enrichment 
  # AIR es el gene set de interes en este caso el adaptative inmune response
Expresion_TG1 <- inner_join(GeneList_TG1,AIR_GeneSet, by="Genes")
Expresion_TG3 <- inner_join(GeneList_TG3,AIR_GeneSet, by="Genes")
rownames(Expresion_TG1) <- Expresion_TG1$Genes
rownames(Expresion_TG3) <- Expresion_TG3$Genes

GenesTG1 <- Expresion_TG1$Genes
GenesTG3 <- Expresion_TG3$Genes

#Aqui dibujamos el diagrama de Venn: 
venn.diagram(
  x = list(GenesTG1, GenesTG3),
  category.names = c("Grupo 1" , "Grupo3 "),
  filename = "Comparacion_Grupo_1_LUAD_frente_a_Grupo_3_LUAD_GSEATOTAL.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.5), alpha('#21908dff',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

# RESULTADOS:  Tenemos que 158 genes son comunes entre ambos grupos mientras que el GRUPO 1 tiene 87 genes UNICOS 
# y el GRUPO 3 tiene 95 genes UNICOS. 

#Con esta funcion obtenemos el conjunto de genes que estan en cada circunferencia. 
OverLap_TG1_TG3 <- calculate.overlap(x=list("Grupo1"=GenesTG1,
                                            "Grupo 3"=GenesTG3))

#El resultado de overlap tiene: 
  # a1) Es el conjunto de genes del grupo 1
  # a2) Es el conjunto de genes del grupo 3
  # a3) Es el conjunto de genes en los que hay interseccion entre ambos 

a3 <- OverLap_TG1_TG3$a3
a2 <- OverLap_TG1_TG3$a2 
a1 <- OverLap_TG1_TG3$a1 

#En este bucle lo que hacemos es de la circunferencias a1 y a2 eliminamos aquellos genes que 
#se encuentran en la interseccion para quedarnos con los genes unicos de cada uno de los conjuntos. 
#Con la funcion which identificamos el indice del gen en mis vectores con los nombres de los genes 
# y luego con eso indexo. 

i <- 1

for(i in seq_along(a3)){
  Gen <- a3[i]
  
  IndiceGen_Grupo3<- which(a2==Gen)
  IndiceGen_Grupo1 <- which(a1==Gen)
  
  a2 <- a2[-IndiceGen_Grupo3]
  a1 <- a1[-IndiceGen_Grupo1]
}

#Ahora como tenemos que en los data.frame de expresion el nombre de las filas es el de los genes podemos idnexar por el nombre del gen. 
#De esta forma tenemos tanto los duplicados, como los genes unicos de cada uno de los conjuntos. 
Expresion_GenesUnicos_TG1 <- Expresion_TG1[a1,]
Expresion_GenesUnicos_TG3 <- Expresion_TG3[a2,]

#Para obtener el conjunto: 

  # 1) Añadimos una columna al data.frame de expresion de la pertenencia a cada grupo 

Expresion_TG1 <- Expresion_TG1 %>% mutate(Grupo=rep("Grupo 1",nrow(Expresion_TG1)))
Expresion_TG3 <- Expresion_TG3 %>% mutate(Grupo=rep("Grupo 3",nrow(Expresion_TG3)))

  # 2) Los juntamos todos en el mismo data.frame 

ExpresionTotal <- rbind(Expresion_TG1,Expresion_TG3)

  # 3) Indexamos con el nombre de los genes

Duplicados_TG1 <- Expresion_TG1[a3,] #Solo son los del grupo 1
Duplicados_TG3 <- Expresion_TG3[a3,]

DuplicadosTotales <- rbind(Duplicados_TG1,Duplicados_TG3)

#Ahora vamos a anotar los genes unicos: 

AnotacionesUnicos_TG1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = Expresion_GenesUnicos_TG1$Genes,
                               mart = ensembl)
AnotacionesUnicos_TG3 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = Expresion_GenesUnicos_TG3$Genes,
                               mart = ensembl)
AnotacionesDuplicados <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = "ensembl_gene_id",
                               values = DuplicadosTotales$Genes,
                               mart = ensembl)

colnames(Expresion_GenesUnicos_TG1)[2] <- "ensembl_gene_id"
colnames(Expresion_GenesUnicos_TG3)[2] <- "ensembl_gene_id"
colnames(DuplicadosTotales)[2] <- "ensembl_gene_id"

Anotaciones_NivelesExpresion_TG1 <- inner_join(AnotacionesUnicos_TG1,Expresion_GenesUnicos_TG1,by="ensembl_gene_id")
Anotaciones_NivelesExpresion_TG3 <- inner_join(AnotacionesUnicos_TG3,Expresion_GenesUnicos_TG3,by="ensembl_gene_id")
Anotaciones_NivelesExpresion_Duplicados <- inner_join(AnotacionesDuplicados,DuplicadosTotales,by="ensembl_gene_id")

# Guardamos Elementos de interes

saveRDS(Anotaciones_NivelesExpresion_TG1,"Comparacion6/Anotaciones_NivelesExpresion_TG1.rds")
saveRDS(Anotaciones_NivelesExpresion_TG3,"Comparacion6/Anotaciones_NivelesExpresion_TG3.rds")
saveRDS(Anotaciones_NivelesExpresion_Duplicados,"Comparacion6/Anotaciones_NivelesExpresion_Duplicados.rds")

