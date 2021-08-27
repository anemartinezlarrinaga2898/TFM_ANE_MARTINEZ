# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-03-2021


# OBJETIVO: Generacion de una matriz de cuentas, de LUAD donde el nombre de las columnas sea el case ID.

##########################################################################################################################################

directory <- setwd( "/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/2_ScriptsPreProcesamiento/LUAD/Datos/")
library(tidyverse)

##########################################################################################################################################

# Cargamos los archivos, todos aquellos que tengan un *count.  
files <- list.files(pattern="*.counts",full.names = TRUE) 

#Limpiamos el nombre de los archivos para poder trabajar luego ya que al principio tienen un ./ que me va a molestar a la hora de generar la matriz 
files <- str_replace_all(files,pattern = "^./",replacement = "")

#Cargamos el manifiesto completo: 

DatosLuad_Manifiesto <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/ManifiestoFinales/ManifiestoLuad_Completo_Unico.txt")

#Relacionar los files con el manifiesto para poder obtener los case ID
  #Primero generamo un data.frame con la lista de los files. 
secuenciasMuestrasTotales <- seq(from=1,to=61,by=1)
DatosFiles <- data.frame(NumeroPaciente=secuenciasMuestrasTotales,
                         file_name=files)
#unimos los files a nuestro manifiesto para saber donde se localiza cada uno 
  #Como estan desorganizados, cogemos el orden establecido en el manifiesto. 
UnionDatos <- inner_join(DatosLuad_Manifiesto,DatosFiles,by="file_name")
  #Vector con el nombre de todos los files
filesDescargar <- UnionDatos$file_name
  #Leemos los files
dataList <- lapply(filesDescargar,read.table, header = FALSE)

#Realizamos un bucle para generar primero la matriz de cuentas donde vamos a poner cada paciente en una columna
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

#Colocamos el nombre de los genes en las filas
rownames(unionLuad) <- unionLuad$V1
unionLuad <- unionLuad[,-1]

#Colocamos el nombre de las columnas el CASE ID de la tcga
i <- 1
for (i in seq(filesDescargar)){
  nombre <- UnionDatos$case_ID[i]
  names(unionLuad)[i] <- nombre
}
  

  

  