# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 29-06-2021

#OBJETIVO: Obtencion de un nuevo manifiesto del FPKM de mis pacientes. 

#####################################################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)
#####################################################################################################################################################################

manifiesto <- read.delim("Manifiesto_Luad_FPKM.txt")
sampleSheet <- read.delim("SampleSheet_Luad.tsv")

names(sampleSheet) <- c("file_ID","file_name","data_category","data_type","proyect_ID","case_ID","sample_ID","sample_type")
names(manifiesto) <- c("file_ID","file_name","md5","size","state")

UnionDatos <- inner_join(sampleSheet,manifiesto,by=c("file_ID","file_name"))
ManifiestoLusc_Completo_Unico<- UnionDatos %>% distinct(case_ID, .keep_all = TRUE)

Huella1 <- readRDS("K3.1_NMF_Luad.rds")
names(Huella1)[1] <- "case_ID"
Huella2 <- readRDS("K3.2_NMF_Luad.rds")
names(Huella2)[1] <- "case_ID"
Huella3 <- readRDS("K3.3_NMF_Luad.rds")
names(Huella3)[1] <- "case_ID"

ManifiestoHuella1 <- inner_join(ManifiestoLusc_Completo_Unico,Huella1,by="case_ID")
ManifiestoHuella1 <- ManifiestoHuella1[,-c(13,12)]

ManifiestoHuella3 <- inner_join(ManifiestoLusc_Completo_Unico,Huella3,by="case_ID")
ManifiestoHuella3 <- ManifiestoHuella3[,-c(13,12)]

ManifiestoHuella2 <- inner_join(ManifiestoLusc_Completo_Unico,Huella2,by="case_ID")
ManifiestoHuella2 <- ManifiestoHuella2[,-c(13,12)]

ManifiestoTotal <- rbind(ManifiestoHuella1,ManifiestoHuella2)
ManifiestoTotal <- rbind(ManifiestoTotal,ManifiestoHuella3)

write_delim(ManifiestoTotal,"Manifiesto_Huella_1_2_3_Luad.txt",delim="\t") #Estan que los primeros 18 son del grupo 1 y los segundos 25 son del grupo 3

#####################################################################################################################################################################

manifiesto <- read.delim("Manifiesto_Lusc_FPKM.txt")
sampleSheet <- read.delim("SampleSheet_Lusc.tsv")

names(sampleSheet) <- c("file_ID","file_name","data_category","data_type","proyect_ID","case_ID","sample_ID","sample_type")
names(manifiesto) <- c("file_ID","file_name","md5","size","state")

UnionDatos <- inner_join(sampleSheet,manifiesto,by=c("file_ID","file_name"))
ManifiestoLusc_Completo_Unico<- UnionDatos %>% distinct(case_ID, .keep_all = TRUE)

Huella1 <- readRDS("K3.1_NMF_Lusc.rds")
names(Huella1)[1] <- "case_ID"
Huella2 <- readRDS("K3.2_NMF_Lusc.rds")
names(Huella2)[1] <- "case_ID"
Huella3 <- readRDS("K3.3_NMF_Lusc.rds")
names(Huella3)[1] <- "case_ID"

ManifiestoHuella1 <- inner_join(ManifiestoLusc_Completo_Unico,Huella1,by="case_ID")
ManifiestoHuella1 <- ManifiestoHuella1[,-c(13,12)]

ManifiestoHuella2 <- inner_join(ManifiestoLusc_Completo_Unico,Huella2,by="case_ID")
ManifiestoHuella2 <- ManifiestoHuella2[,-c(13,12)]

ManifiestoHuella3 <- inner_join(ManifiestoLusc_Completo_Unico,Huella3,by="case_ID")
ManifiestoHuella3 <- ManifiestoHuella3[,-c(13,12)]

ManifiestoTotal <- rbind(ManifiestoHuella1,ManifiestoHuella2)
ManifiestoTotal <- rbind(ManifiestoTotal,ManifiestoHuella3)

write_delim(ManifiestoTotal,"Manifiesto_Huella_1_2_3_Lusc.txt",delim="\t") #Estan que los primeros 18 son del grupo 1 y los segundos 25 son del grupo 3

##############################################################################################################################

ManifiestoReferencia <- read.delim("Manifiesto_Luad_FPKM.txt")

ManifiestoTotal_Luad <- read.delim("Manifiesto_Huella_1_2_3_Luad.txt")
ManifiestoTotal_Lusc <- read.delim("Manifiesto_Huella_1_2_3_Lusc.txt")

NombresManifiesto <- colnames(ManifiestoReferencia)

ManifiestoLuad <- ManifiestoTotal_Luad[,c(1,2,9,10,11)]
ManifiestoLusc <- ManifiestoTotal_Lusc[,c(1,2,9,10,11)]

colnames(ManifiestoLuad) <- NombresManifiesto
colnames(ManifiestoLusc) <- NombresManifiesto

write_delim(ManifiestoLuad,"ManifiestoLuad.txt",delim="\t")
write_delim(ManifiestoLusc,"ManifiestoLusc.txt",delim="\t")
