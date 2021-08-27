# Generacion de una matriz de cuentas

setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/DATA/DatosLuad")

# Cargar los archivos: 
files <- list.files(pattern="*.counts",full.names = TRUE) 
dataList <- lapply(files, read.table, header = FALSE)

tamañoLista <- length(dataList)
secuenciaMuestras <- seq(from=3,to=tamañoLista,by=1)
secuenciasMuestrasCompleta <- seq(from=0,to=tamañoLista,by=1)

sample1 <- dataList[[1]]
sample2 <- dataList[[2]]

sample <- full_join(sample1,sample2,by="V1")

i=1
for (i in seq_along(secuenciaMuestras)){
  sampleNuevo <- dataList[[i]]
  unionLuad <- full_join(sampleNuevo,sample,by="V1")
  sample <- unionLuad
}


i=1
for (i in seq(secuenciasMuestrasCompleta)){
  nombre <- secuenciasMuestrasCompleta[i]
  names(unionLuad)[i] <- nombre
}
  
names(unionLuad)[1] <- "Gen" 
library(tidyverse)
write_delim(unionLuad,"unionLUAD")
  
  