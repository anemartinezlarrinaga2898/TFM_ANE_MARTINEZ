# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-03-2021

#OBJETIVO: Obtencion del manifiesto junto con los case ID: 

library(tidyverse)

directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos")
manifiestoLusc <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/manifiesto_LUSC_21_03_21.txt")
sampleSheet <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/sample_sheet_archivos_selecionados_lusc_21_03_21.tsv")

names(sampleSheet) <- c("file_ID","file_name","data_category","data_type","proyect_ID","case_ID","sample_ID","sample_type")
names(manifiestoLusc) <- c("file_ID","file_name","md5","size","state")

UnionDatos <- inner_join(sampleSheet,manifiestoLusc,by=c("file_ID","file_name"))
ManifiestoLusc_Completo_Unico<- UnionDatos %>% distinct(case_ID, .keep_all = TRUE)

ManifiestoDescargarLusc <- ManifiestoLusc_Completo_Unico[,c(1,2,9,10,11)]
names(ManifiestoDescargarLusc) <- c("id","filename","md5","size","state")


write_delim(UnionDatos, "ManifiestoLusc_Completo.txt", delim = "\t")
write_delim(ManifiestoLusc_Completo_Unico, "ManifiestoLusc_Completo_Unico.txt", delim = "\t")
write_delim(ManifiestoDescargarLusc,"ManifiestoDescargar_Lusc.txt",delim="\t")

