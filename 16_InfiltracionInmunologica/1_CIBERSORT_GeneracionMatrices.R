# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 30-06 -2021

# OBJETIVO: Generacion de una matriz de cuentas, de LUAD donde el nombre de las columnas sea el case ID.

##########################################################################################################################################

directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
library(tidyverse)

##########################################################################################################################################

directoryLUAD <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/LUAD/")

files <- list.files(pattern="*.FPKM.txt",full.names = TRUE) 
files <- str_replace_all(files,pattern = "^./",replacement = "")

DatosLuad_Manifiesto <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/Manifiesto_Huella_1_2_3_Luad.txt")

secuenciasMuestrasTotales <- seq(from=1,to=61,by=1)
DatosFiles <- data.frame(NumeroPaciente=secuenciasMuestrasTotales,
                         file_name=files)

UnionDatos <- inner_join(DatosLuad_Manifiesto,DatosFiles,by="file_name")

filesDescargar <- UnionDatos$file_name
dataList <- lapply(filesDescargar,read.table, header = FALSE)

tamañoLista <- length(dataList)
secuenciaMuestras <- seq(from=3,to=tamañoLista,by=1)

sample1 <- dataList[[1]]
sample2 <- dataList[[2]]

sample <- full_join(sample1,sample2,by="V1")

i <- 1
for (i in seq(3,61,1)){
  sampleNuevo <- dataList[[i]]
  unionLuad <- full_join(sample,sampleNuevo,by="V1")
  sample <- unionLuad
}

rownames(unionLuad) <- unionLuad$V1
unionLuad <- unionLuad[,-1]

i <- 1
for (i in seq(filesDescargar)){
  nombre <- UnionDatos$case_ID[i]
  names(unionLuad)[i] <- nombre
}

rownames(unionLuad) <- gsub("\\..*","",rownames(unionLuad))
directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
saveRDS(unionLuad,"MatrizLuad_FPKM.rds")

##########################################################################################################################################

directoryLUSC <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/LUSC/")

files <- list.files(pattern="*.FPKM.txt",full.names = TRUE) 
files <- str_replace_all(files,pattern = "^./",replacement = "")

DatosLusc_Manifiesto <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/Manifiesto_Huella_1_2_3_Lusc.txt")

secuenciasMuestrasTotales <- seq(from=1,to=495,by=1)
DatosFiles <- data.frame(NumeroPaciente=secuenciasMuestrasTotales,
                         file_name=files)

UnionDatos <- inner_join(DatosLusc_Manifiesto,DatosFiles,by="file_name")

filesDescargar <- UnionDatos$file_name
dataList <- lapply(filesDescargar,read.table, header = FALSE)

tamañoLista <- length(dataList)
secuenciaMuestras <- seq(from=3,to=tamañoLista,by=1)

sample1 <- dataList[[1]]
sample2 <- dataList[[2]]

sample <- full_join(sample1,sample2,by="V1")

i <- 1
for (i in seq(3,495,1)){
  sampleNuevo <- dataList[[i]]
  unionLuad <- full_join(sample,sampleNuevo,by="V1")
  sample <- unionLuad
}

rownames(unionLuad) <- unionLuad$V1
unionLuad <- unionLuad[,-1]

i <- 1
for (i in seq(filesDescargar)){
  nombre <- UnionDatos$case_ID[i]
  names(unionLuad)[i] <- nombre
}

rownames(unionLuad) <- gsub("\\..*","",rownames(unionLuad))
directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
saveRDS(unionLuad,"MatrizLusc_FPKM.rds")






