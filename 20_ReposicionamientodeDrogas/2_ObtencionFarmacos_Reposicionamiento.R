# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 28-07-2021

#OBJETIVO:Analisis de los resultados del reposicionamiento. 

###############################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/20_ReposicionamientodeDrogas/")
library(tidyverse)
###############################################################################################################################################

#······························································· L U A D ······································································

GRUPO1_LUAD <- read_delim("Resultados/GRUPO1 LUAD.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO1_LUAD <- GRUPO1_LUAD[-1,]
GRUPO1_LUAD <- GRUPO1_LUAD[order(GRUPO1_LUAD$raw_cs,decreasing = TRUE),]

GRUPO1_LUAD_ValoresNegativos <- GRUPO1_LUAD %>% filter(raw_cs<0)
GRUPO1_LUAD_ValoresNegativos <- GRUPO1_LUAD_ValoresNegativos[order(GRUPO1_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G1 <- GRUPO1_LUAD_ValoresNegativos[1:10,]
TOP_10_G1$Grupo <- "GRUPO 1"

GRUPO2_LUAD <- read_delim("Resultados/GRUPO2 LUAD.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO2_LUAD <- GRUPO2_LUAD[-1,]
GRUPO2_LUAD <- GRUPO2_LUAD[order(GRUPO2_LUAD$raw_cs,decreasing = TRUE),]

GRUPO2_LUAD_ValoresNegativos <- GRUPO2_LUAD %>% filter(raw_cs<0)
GRUPO2_LUAD_ValoresNegativos <- GRUPO2_LUAD_ValoresNegativos[order(GRUPO2_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G2 <- GRUPO2_LUAD_ValoresNegativos[1:10,]
TOP_10_G2$Grupo <- "GRUPO 2"

GRUPO3_LUAD <- read_delim("Resultados/GRUPO3 LUAD.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO3_LUAD <- GRUPO3_LUAD[-1,]
GRUPO3_LUAD <- GRUPO3_LUAD[order(GRUPO3_LUAD$raw_cs,decreasing = TRUE),]

GRUPO3_LUAD_ValoresNegativos <- GRUPO3_LUAD %>% filter(raw_cs<0)
GRUPO3_LUAD_ValoresNegativos <- GRUPO3_LUAD_ValoresNegativos[order(GRUPO3_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G3 <- GRUPO3_LUAD_ValoresNegativos[1:10,]
TOP_10_G3$Grupo <- "GRUPO 3"

TOP_30_LUAD <- rbind(TOP_10_G1,TOP_10_G2)
TOP_30_LUAD <- rbind(TOP_30_LUAD,TOP_10_G3)

Farmacos <- TOP_30_LUAD[,c(13,14,17)]
Farmacos$desc <- gsub(":.*","",Farmacos$desc)
Farmacos <- Farmacos[order(Farmacos$desc),]

FarmacosDuplicados <- Farmacos[duplicated(Farmacos$desc),]

Farmacos$raw_cs <- as.numeric(Farmacos$raw_cs)
Farmacos$Grupo <- as.factor(Farmacos$Grupo)

write.table(Farmacos,"FARMACOS_LUAD_TOP30.csv",sep="\tab",quote = FALSE,row.names = FALSE)
write.table(TOP_30_LUAD,"FARMACOS_LUAD_TOP30_Completo.csv",sep="\tab",quote = FALSE,row.names = FALSE)

#······························································· L U S C ······································································

GRUPO1_LUAD <- read_delim("Resultados/GRUPO1 LUSC.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO1_LUAD <- GRUPO1_LUAD[-1,]
GRUPO1_LUAD <- GRUPO1_LUAD[order(GRUPO1_LUAD$raw_cs,decreasing = TRUE),]

GRUPO1_LUAD_ValoresNegativos <- GRUPO1_LUAD %>% filter(raw_cs<0)
GRUPO1_LUAD_ValoresNegativos <- GRUPO1_LUAD_ValoresNegativos[order(GRUPO1_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G1 <- GRUPO1_LUAD_ValoresNegativos[1:10,]
TOP_10_G1$Grupo <- "GRUPO 1"

GRUPO2_LUAD <- read_delim("Resultados/GRUPO2 LUSC.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO2_LUAD <- GRUPO2_LUAD[-1,]
GRUPO2_LUAD <- GRUPO2_LUAD[order(GRUPO2_LUAD$raw_cs,decreasing = TRUE),]

GRUPO2_LUAD_ValoresNegativos <- GRUPO2_LUAD %>% filter(raw_cs<0)
GRUPO2_LUAD_ValoresNegativos <- GRUPO2_LUAD_ValoresNegativos[order(GRUPO2_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G2 <- GRUPO2_LUAD_ValoresNegativos[1:10,]
TOP_10_G2$Grupo <- "GRUPO 2"

GRUPO3_LUAD <- read_delim("Resultados/GRUPO3 LUSC.gct", 
                          "\t", escape_double = FALSE, 
                          locale = locale(decimal_mark = ",",grouping_mark = "."), 
                          trim_ws = TRUE, 
                          skip = 2)
GRUPO3_LUAD <- GRUPO3_LUAD[-1,]
GRUPO3_LUAD <- GRUPO3_LUAD[order(GRUPO3_LUAD$raw_cs,decreasing = TRUE),]

GRUPO3_LUAD_ValoresNegativos <- GRUPO3_LUAD %>% filter(raw_cs<0)
GRUPO3_LUAD_ValoresNegativos <- GRUPO3_LUAD_ValoresNegativos[order(GRUPO3_LUAD_ValoresNegativos$raw_cs,decreasing=TRUE),]

TOP_10_G3 <- GRUPO3_LUAD_ValoresNegativos[1:10,]
TOP_10_G3$Grupo <- "GRUPO 3"

TOP_30_LUAD <- rbind(TOP_10_G1,TOP_10_G2)
TOP_30_LUAD <- rbind(TOP_30_LUAD,TOP_10_G3)

Farmacos <- TOP_30_LUAD[,c(13,14,17)]
Farmacos$desc <- gsub(":.*","",Farmacos$desc)
Farmacos <- Farmacos[order(Farmacos$desc),]

FarmacosDuplicados <- Farmacos[duplicated(Farmacos$desc),]

Farmacos$raw_cs <- as.numeric(Farmacos$raw_cs)
Farmacos$Grupo <- as.factor(Farmacos$Grupo)

write.table(Farmacos,"FARMACOS_LUSC_TOP30.csv",sep="\tab",quote = FALSE,row.names = FALSE)
write.table(TOP_30_LUAD,"FARMACOS_LUSC_TOP30_Completo.csv",sep="\tab",quote = FALSE,row.names = FALSE)


