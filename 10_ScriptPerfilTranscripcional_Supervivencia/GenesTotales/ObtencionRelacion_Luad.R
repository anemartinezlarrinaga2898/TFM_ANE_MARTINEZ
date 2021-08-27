# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:GUna vez tenemos las matrices de cuentas generaas vamos a sumar la epxresion de los genes en los diferentes genes. 
#Tambien vamos a generar las matrices porque ahora las tengo en una lista 
#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
#########################################################################################################

MatrizGenesPerfilTranscripcional_Luad <- readRDS("MatrizGenesPerfilTranscripcional_Luad.rds")

#·······························  L U A D ·······························································

Matriz_Luad_1 <- MatrizGenesPerfilTranscripcional_Luad[[1]]
Matriz_Luad_2 <- MatrizGenesPerfilTranscripcional_Luad[[2]]
Matriz_Luad_3 <- MatrizGenesPerfilTranscripcional_Luad[[3]]

Matriz_Luad_1 <- t(Matriz_Luad_1)
Matriz_Luad_2 <- t(Matriz_Luad_2)
Matriz_Luad_3 <- t(Matriz_Luad_3)

saveRDS(Matriz_Luad_1,"MatrizTranscripcion_Luad_1.rds")
saveRDS(Matriz_Luad_2,"MatrizTranscripcion_Luad_2.rds")
saveRDS(Matriz_Luad_3,"MatrizTranscripcion_Luad_3.rds")

Sumatorio_Matriz_1 <- t(rowSums(Matriz_Luad_1))
Sumatorio_Matriz_2 <- t(rowSums(Matriz_Luad_2))
Sumatorio_Matriz_3 <- t(rowSums(Matriz_Luad_3))

Sumatorio_Matriz_1 <- round(Sumatorio_Matriz_1,2)
Sumatorio_Matriz_2 <- round(Sumatorio_Matriz_2,2)
Sumatorio_Matriz_3 <- round(Sumatorio_Matriz_3,2)

saveRDS(Sumatorio_Matriz_1,"Sumatorio_Matriz_1_Luad.rds")
saveRDS(Sumatorio_Matriz_2,"Sumatorio_Matriz_2_Luad.rds")
saveRDS(Sumatorio_Matriz_3,"Sumatorio_Matriz_3_Luad.rds")

Media_1 <- round(mean(Sumatorio_Matriz_1),2)
Mediana_1 <- median(Sumatorio_Matriz_1)

Media_2 <- round(mean(Sumatorio_Matriz_2),2)
Mediana_2 <- median(Sumatorio_Matriz_2)

Media_3 <- round(mean(Sumatorio_Matriz_3),2)
Mediana_3 <- median(Sumatorio_Matriz_3)

MediasTotales <- c(Media_1,Media_2,Media_3)
MedianasTotales <- c(Mediana_1,Mediana_2,Mediana_3)

# BUCLE PARA VER SI ESTA POR ENCIMA O POR DEBAJO TANTO DE LA MEDIA COMO DE LA MEDIANA: 

Informacion_1 <- data.frame(Media=rep(0,length(Sumatorio_Matriz_1)),
                          Mediana=rep(0,length(Sumatorio_Matriz_1)))
Informacion_2 <- data.frame(Media=rep(0,length(Sumatorio_Matriz_1)),
                            Mediana=rep(0,length(Sumatorio_Matriz_1)))
Informacion_3 <- data.frame(Media=rep(0,length(Sumatorio_Matriz_1)),
                            Mediana=rep(0,length(Sumatorio_Matriz_1)))

rownames(Informacion_1) <- rownames(Matriz_Luad_1)
rownames(Informacion_2) <- rownames(Matriz_Luad_2)
rownames(Informacion_3) <- rownames(Matriz_Luad_3)

TotalSumatorio <- list(Sumatorio_Matriz_1,Sumatorio_Matriz_2,Sumatorio_Matriz_3)
TotalInformacion <- list(Informacion_1,Informacion_2,Informacion_3)

i <- 1

for ( i in seq_along(TotalSumatorio)){
  
  Sumatorio <- TotalSumatorio[[i]]
  Informacion <- TotalInformacion[[i]]
  Media <- MediasTotales[i]
  Mediana <- MedianasTotales[i]
  
  j <- 1
  
  for (j in seq_along(Sumatorio)){
    Valor <- Sumatorio[j]
    
    if (Valor > Media){Informacion[j,1] <- 1}
    else{Informacion[j,1] <- 0}
    
    if(Valor > Mediana){Informacion[j,2] <- 1}
    else(Informacion[j,2] <- 0)
  }
  
  if (i == 1){TotalInformacion[[i]] <- Informacion}
  if (i == 2){TotalInformacion[[i]] <- Informacion}
  if (i == 3){TotalInformacion[[i]] <- Informacion}
  
}

saveRDS(TotalInformacion,"InformacionTotal_SinClinica_Luad.rds")

# OBNTENCION DE LA INFORMACION CLINICA: 

Informacion <- readRDS("InformacionTotal_SinClinica_Luad.rds")
Clinica <- readRDS("ClinicalData_PacientesLuad_cBio.rds")

Clinica <- as.data.frame(Clinica[,c(2,67,68)])
colnames(Clinica) <- gsub(" ","",colnames(Clinica))
Clinica$OverallSurvivalStatus <- gsub("\\:.*","",Clinica$OverallSurvivalStatus)

Informacion_1 <- Informacion[[1]]
Informacion_2 <- Informacion[[2]]
Informacion_3 <- Informacion[[3]]

PacientesInformacion <- rownames(Informacion_1)

Informacion_1 <- cbind(Informacion_1,PacientesInformacion)
Informacion_2 <- cbind(Informacion_2,PacientesInformacion)
Informacion_3 <- cbind(Informacion_3,PacientesInformacion)

Informacion_1 <- Informacion_1[,c(3,1,2)]
Informacion_2 <- Informacion_2[,c(3,1,2)]
Informacion_3 <- Informacion_3[,c(3,1,2)]

colnames(Informacion_1)[1] <- colnames(Clinica)[1]
colnames(Informacion_2)[1] <- colnames(Clinica)[1]
colnames(Informacion_3)[1] <- colnames(Clinica)[1]

Informacion_Clinica_1 <- full_join(Informacion_1,Clinica,by=colnames(Clinica)[1])
Informacion_Clinica_2 <- full_join(Informacion_2,Clinica,by=colnames(Clinica)[1])
Informacion_Clinica_3 <- full_join(Informacion_3,Clinica,by=colnames(Clinica)[1])

saveRDS(Informacion_Clinica_1,"Informacion_Clinica_Luad_1.rds")
saveRDS(Informacion_Clinica_2,"Informacion_Clinica_Luad_2.rds")
saveRDS(Informacion_Clinica_3,"Informacion_Clinica_Luad_3.rds")


