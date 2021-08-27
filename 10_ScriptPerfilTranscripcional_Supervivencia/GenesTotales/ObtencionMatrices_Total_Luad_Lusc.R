# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 13-05-2021

#OBJETIVO:Generar las matricecs de cuentas donde tenga los genes filtrados : 
# Con los genes que tengo de el filtro dos los cojo y los obtengo de mi matriz de cuentas. 

#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")

Luad <- readRDS("Luad_Filtrado.rds")
Lusc <- readRDS("Lusc_Filtrado.rds")

#Hay que eliminar el .Version del ENSEMBEL porque si no, no van a encajar. 
rownames(Luad) <- gsub("\\..*","",rownames(Luad))
rownames(Lusc) <- gsub("\\..*","",rownames(Lusc))


#························· L U A D ············································

#Filtros 2: p.adjust de 0.05 y un LFC de [-infinito,-2]U[2,+infinito]

Filtrados_1 <- readRDS("ResultadosFiltrados_Luad_Comparacion1.rds")
Filtrados_2 <- readRDS("ResultadosFiltrados_Luad_Comparacion2.rds")
Filtrados_3 <- readRDS("ResultadosFiltrados_Luad_Comparacion3.rds")

#Del filtrao queremos el nombre de los genes: 
Genes_1 <- rownames(Filtrados_1) 
Genes_2 <- rownames(Filtrados_2)
Genes_3 <- rownames(Filtrados_3)

GenesTotales <- list(Genes_1,Genes_2,Genes_3)
MatrizTotal <- list()
i <- 1
j <- 1

for ( i in seq_along(GenesTotales)){
  Genes <- GenesTotales[[i]]
  MatrizGenesPerfilTranscripcional <- data.frame()
  for ( j in seq_along(Genes)){
    Gen <- Genes[j]
    GenPaciente <- Luad[which(rownames(Luad)==Gen),]
    MatrizGenesPerfilTranscripcional <- rbind(MatrizGenesPerfilTranscripcional,GenPaciente)
  }
  
  if (i==1){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
  if (i==2){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
  if (i==3){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
}

saveRDS(MatrizTotal,"MatrizGenesPerfilTranscripcional_Luad.rds")

#························· L U S C ············································

Filtrados_1 <- readRDS("ResultadosFiltrados_Lusc_Comparacion1.rds")
Filtrados_2 <- readRDS("ResultadosFiltrados_Lusc_Comparacion2.rds")
Filtrados_3 <- readRDS("ResultadosFiltrados_Lusc_Comparacion3.rds")

Genes_1 <- rownames(Filtrados_1)
Genes_2 <- rownames(Filtrados_2)
Genes_3 <- rownames(Filtrados_3)

GenesTotales <- list(Genes_1,Genes_2,Genes_3)

MatrizGenesPerfilTranscripcional <- data.frame()
MatrizTotal <- list()
i <- 1
j <- 1

for ( i in seq_along(GenesTotales)){
  Genes <- GenesTotales[[i]]
  MatrizGenesPerfilTranscripcional <- data.frame()
  for ( j in seq_along(Genes)){
    Gen <- Genes[j]
    GenPaciente <- Lusc[which(rownames(Lusc)==Gen),]
    MatrizGenesPerfilTranscripcional <- rbind(MatrizGenesPerfilTranscripcional,GenPaciente)
  }
  
  if (i==1){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
  if (i==2){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
  if (i==3){MatrizTotal[[i]] <- MatrizGenesPerfilTranscripcional}
}

saveRDS(MatrizTotal,"MatrizGenesPerfilTranscripcional_Lusc.rds")


