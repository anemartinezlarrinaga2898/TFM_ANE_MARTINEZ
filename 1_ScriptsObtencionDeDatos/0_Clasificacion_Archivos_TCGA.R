# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 21-03-2021

# TITULO: Clasificacion de los diferentes datos que ofrece la TCGA: 
  #Objetivo_Script_ El objetivo de este script es, clasificar los diferentes documentos que puedes descargarte
  #de la TCGA. Vamos a trabajar con el caso de LUAD pero en el caso de LUSC se haria igual.

#Establecezco el directorio de trabajo: 
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/")

#Primero trabajo con los datos descargados a fecha del 21-03-2021
ManifiestoLuad_Nuevo <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/Manifiesto_Luad_21_03_21.txt", 
                             row.names=NULL)
#Manifiesto tiene 547 files que corresponden a 475 casos

ManifiestoLuad_Viejo <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/Genes/ManifiestoORIGINAL.txt")

#Este sample sheet lo he descargado del carrito donde pone SampleSheet. 
  #A la hora de descargar todos los datos. 

RNA_SEQ_ArchivosSelecionados <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/RNA_SEQ_ArchivosSelecionados.tsv")

#Dentro de la carpeta de Bioespecimen la hoja de sample:(del carrito, el bioespecimen que se descarga desde la hoja del carrito)

sample_carrito <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/CarpetaBioespecimen_Carrito/sample.tsv")

#El bioespecimen que me descargo desde el repository: 

sample_repository <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/CarpetaBioespecimen_Repository/sample_21_03_2021_todos_los_archivos.tsv")

#Todos los archivos relacionados con transcriptome profiling: 

TranscrioptomeProfiling <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/NuevasCuentas/Archivos_TCGA/LUAD/TranscrioptomeProfiling.tsv")

#Una vez que tenemos los datos voy a intentar ver si en el transcriptome profiling tengo los datos del manifiesto original
  #Si los datos estan no tengo que empezar desde cero
  #Si los datos no estan, tengo que empezar desde cero

#Los archivos que necesito son: 

TranscrioptomeProfiling #Voy a selecionar las columnas que me interesan ( luego si sale bien tendre que cambiarle los nombres)
ManifiestoLuad_Viejo

TranscriptomeProfilingPrueba <- TranscrioptomeProfiling[,c(1,2,6,7)]
names(TranscriptomeProfilingPrueba) <- c("file_ID","file_name","case_ID","sample_ID")
names(ManifiestoLuad_Viejo) <- c("file_ID","file_name","md5","size","state")

TranscriptomeProfilingPruebaParteDeLosDatos <- TranscriptomeProfilingPrueba[1:10,]

library(tidyverse)
#Parte de los datos
UnionFull_joinParte <- full_join(TranscriptomeProfilingPruebaParteDeLosDatos,ManifiestoLuad_Viejo,by=c("file_ID","file_name"))
#UnionConTodosLosDatos: 
UnionInner_join <- inner_join(TranscriptomeProfilingPrueba,ManifiestoLuad_Viejo,by=c("file_ID","file_name"))
UnionFull_joinTodo <- full_join(RNA_SEQ_ArchivosSelecionados,ManifiestoLuad_Viejo,by=c("file_ID","file_name"))

caseIDNombres <- UnionInner_join[,3]
nombresUnicos <- unique(caseIDNombres)


#Quiero probar con otros manifiesto: 
  #Estan Filtrado
    
Manifiesto_551_Files <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/Genes/Manifiesto_551_Files.txt")
Manifiesto_547_Files <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/Genes/Manifiesto_547_Files.txt")

Luad_Filtrado_70 <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/Genes/Luad_Filtrado.txt")
Luad_Filtrado_Sin_Fusiones <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/TXT_Cuentas/Genes/Luad_Filtrado_Sin_Fusiones.txt")


