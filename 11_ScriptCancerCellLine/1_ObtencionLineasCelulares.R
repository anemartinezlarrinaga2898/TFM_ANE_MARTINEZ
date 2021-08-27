# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 25-05-2021

#OBJETIVO:Script para obtener las lineas celulares de mis huellas. 

################################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")

################################################################################################################################################
#Cargamos los datos
HUELLAS <- readRDS("Huellas_LineasCelulares_Lusc.rds")
LINES <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/Cell_lines_annotations_20181226.txt")

#Buscamos cuantas lineas UNICAS hay para concer sus nombres
AnotacionesLineas <- unique(LINES$PATHOLOGIST_ANNOTATION)

#Obtenemos las lineas que corresponden a LUSC
LineasLusc <- LINES[which(LINES$PATHOLOGIST_ANNOTATION=="Lung:NSCLC_Squamous"),]

#Los ID de las lineas celulares de LUSC que nos interesan
NombresLineas <- LineasLusc$CCLE_ID

#Obtener las columnas en mis huellas: 
i <- 1
Huellas_Lineas_Lusc <- list()
Nombres <- c("Grupo 1","Grupo 2","Grupo 3")

for ( i in seq_along(HUELLAS)){
  Huella <- HUELLAS[[i]]
  
  j <- 1
  
  X <- data.frame(X=rep(0,nrow(Huella))) #Es un objeto intermedio para poder ir agrupando las lineas
  rownames(X) <- rownames(Huella) #Para poder luego hacer el cbind y mantener el ID de los genes con los que estamos trabajando
  X <- X[,-1]
  for ( j in seq_along(NombresLineas)){
    Linea <- NombresLineas[j]
    
    linea_huella <- Huella[,which(colnames(Huella)==Linea)] #Aquie encontramos la columna de nustra linea de interes
    
    X <- cbind(X,linea_huella) #Para unrilo con lo anterior
    colnames(X)[j] <- Linea #Para mantener el ID de la cell line con la que estamos trabajando
  }
  
  X <- t(X) #Tenemos que usar la traspuesta porque queremos estandarizar los valores de un gen en las diferentes lineas celulares. 
  X <- scale(X,center=TRUE,scale=TRUE)
  X <- t(X)
  
  Huellas_Lineas_Lusc[[i]] <- X
  names(Huellas_Lineas_Lusc)[[i]] <- Nombres[i]
}

mean(Huellas_Lineas_Lusc[[1]][1,])
sd(Huellas_Lineas_Lusc[[1]][1,])

saveRDS(Huellas_Lineas_Lusc,"Huellas_LineasCelulares_Lusc.rds")
