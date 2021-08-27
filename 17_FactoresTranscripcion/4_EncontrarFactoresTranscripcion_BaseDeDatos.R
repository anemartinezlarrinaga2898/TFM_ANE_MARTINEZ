# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 07-07-2021

#OBJETIVO: Generar los data.frames del analisis DESEQ y luego comprobar que factores de transcripcion tenemos. 
###########################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/19_AnalisisMatricesOriginales/")
library(tidyverse)
###############################################################################################

G1 <- readRDS("ResultadosAnalasisExpresion_1_Luad.rds")
G1 <- as.data.frame(dplyr::mutate(as.data.frame(G1)),
                                     row.names=rownames(G1))
G1$ensembl <- rownames(G1)

G2 <- readRDS("ResultadosAnalasisExpresion_2_Luad.rds")
G2 <- as.data.frame(dplyr::mutate(as.data.frame(G2)),
                    row.names=rownames(G2))
G2$ensembl <- rownames(G2)

G3 <- readRDS("ResultadosAnalasisExpresion_3_Luad.rds")
G3 <- as.data.frame(dplyr::mutate(as.data.frame(G3)),
                    row.names=rownames(G3))
G3$ensembl <- rownames(G3)

FT <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
FT <- FT[,-1]
colnames(FT)[3] <- "ensembl"

FT_G1 <- inner_join(G1,FT,by="ensembl")
FT_G1$Grupo <- "Grupo1"
FT_G2 <- inner_join(G2,FT,by="ensembl")
FT_G2$Grupo <- "Grupo2"
FT_G3 <- inner_join(G3,FT,by="ensembl")
FT_G3$Grupo <- "Grupo3"

FT_Totales <- rbind(FT_G1,FT_G2)
FT_Totales <- rbind(FT_Totales,FT_G3)

FT_Totales <- FT_Totales[order(FT_Totales$SYMBOL),]

FT_UP_G1 <- FT_G1[order(FT_G2$log2FoldChange,decreasing = TRUE),]
FT_UP_G1_5 <- FT_UP_G1[1:5,]

FT_down_G1 <- FT_G1[order(FT_G2$log2FoldChange,decreasing = FALSE),]
FT_down_G1 <- FT_down_G1[1:5,]

FT_G1_UP_DOWN <- rbind(FT_UP_G1_5,FT_down_G1)
write.table(FT_G1_UP_DOWN,"FT_G2_UP_DOWN.csv",sep ="\t",row.names = FALSE)


saveRDS(FT_Totales,"FT_TOTALES_LUAD.rds")

#·····························································································································

G1 <- readRDS("ResultadosAnalasisExpresion_1_Lusc.rds")
G1 <- as.data.frame(dplyr::mutate(as.data.frame(G1)),
                    row.names=rownames(G1))
G1$ensembl <- rownames(G1)

G2 <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
G2 <- as.data.frame(dplyr::mutate(as.data.frame(G2)),
                    row.names=rownames(G2))
G2$ensembl <- rownames(G2)

G3 <- readRDS("ResultadosAnalasisExpresion_3_Lusc.rds")
G3 <- as.data.frame(dplyr::mutate(as.data.frame(G3)),
                    row.names=rownames(G3))
G3$ensembl <- rownames(G3)

FT <- readRDS("FactoresTranscripcion_Symbol_Entrez_Ensembel.rds")
FT <- FT[,-1]
colnames(FT)[3] <- "ensembl"

FT_G1 <- inner_join(G1,FT,by="ensembl")
FT_G1$Grupo <- "Grupo1"
FT_G2 <- inner_join(G2,FT,by="ensembl")
FT_G2$Grupo <- "Grupo2"
FT_G3 <- inner_join(G3,FT,by="ensembl")
FT_G3$Grupo <- "Grupo3"

FT_Totales <- rbind(FT_G1,FT_G2)
FT_Totales <- rbind(FT_Totales,FT_G3)

FT_Totales <- FT_Totales[order(FT_Totales$SYMBOL),]

saveRDS(FT_Totales,"FT_TOTALES_LUSC.rds")

##############################################################################################################################

#Ahora buscamos en LUSC los FT que son de cada grupo en toda la huella. 

FT <- readRDS("FT_TOTALES_LUSC.rds")
FT_Paper <- readRDS("GenesPaper_Lusc.rds")
FT_Paper[21,2] <-"NKX2-1"
colnames(FT_Paper)[2] <- "Symbol"

Paper <- inner_join(FT,FT_Paper, by="Symbol")

##############################################################################################################################
#Buscar los genes en todos los genes que se encuentran diferencialmente expresados. 

GP <- readRDS("GenesPaper_Lusc.rds")
GP[21,2] <-"NKX2-1"
colnames(GP)[2] <- "SYMBOL"

H1 <- readRDS("ResultadosComparacion1_Lusc_DF.rds")
H1$ensembl <- rownames(H1)
H1$Grupo <- "Grupo1"

H2 <- readRDS("ResultadosComparacion2_Lusc_DF.rds")
H2$ensembl <- rownames(H2)
H2$Grupo <- "Grupo2"

H3 <- readRDS("ResultadosComparacion3_Lusc_DF.rds")
H3$ensembl <- rownames(H3)
H3$Grupo <- "Grupo3"

H <- rbind(H1,H2)
H <- rbind(H,H3)

PAPER <- inner_join(H,GP,by="SYMBOL")
PAPER <- PAPER[order(PAPER$SYMBOL),]

HP1 <- PAPER %>% filter(Grupo=="Grupo1")
HP2 <- PAPER %>% filter(Grupo=="Grupo2")
HP3 <- PAPER %>% filter(Grupo=="Grupo3")

###############################################################################################

H1 <- readRDS("ResultadosComparacion1_Lusc_DF.rds")
H1$ensembl <- rownames(H1)
H1$Grupo <- "Grupo1"

H2 <- readRDS("ResultadosComparacion2_Lusc_DF.rds")
H2$ensembl <- rownames(H2)
H2$Grupo <- "Grupo2"

H3 <- readRDS("ResultadosComparacion3_Lusc_DF.rds")
H3$ensembl <- rownames(H3)
H3$Grupo <- "Grupo3"

H <- rbind(H1,H2)
H <- rbind(H,H3)

GP <- readRDS("GenesPaper_Ampliados.rds")

PAPER <- inner_join(H,GP,by="SYMBOL")
PAPER <- PAPER[order(PAPER$SYMBOL),]

HP1 <- PAPER %>% filter(Grupo=="Grupo1")
HP2 <- PAPER %>% filter(Grupo=="Grupo2")
HP3 <- PAPER %>% filter(Grupo=="Grupo3")


