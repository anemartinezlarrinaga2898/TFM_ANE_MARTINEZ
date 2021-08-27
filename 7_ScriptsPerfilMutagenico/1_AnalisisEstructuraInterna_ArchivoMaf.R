#Script para analizar el archvio MAF
#PRUEBA
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")

MutacionesLuad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/Datos/Luad/gdc_download_20210421_083318.951054/0458c57f-316c-4a7c-9294-ccd11c97c2f9/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf", header=FALSE, comment.char="#")

Columna16 <- unique(MutacionesLuad$V16)

NombresIDLuad <- readRDS("K3_NMF_Luad.rds")

NombresIDLuadPrueba <- "TCGA-99-7458"

Nombre <- NombresIDLuad[which(NombresIDLuadPrueba==NombresIDLuad$Sample_ID),]


#########################################################################################################################
#PRUEBA
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")
MutacionesLuad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/Datos/Luad/gdc_download_20210421_083318.951054/0458c57f-316c-4a7c-9294-ccd11c97c2f9/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf", header=FALSE, comment.char="#")

ColumnaID <- unique(MutacionesLuad$V16)
ColumnaID <- ColumnaID[-1]

Prueba <- ColumnaID[1]
Prueba2 <- ColumnaID[2]
Prueba3 <- ColumnaID[3]

Prueba1 <- gsub("(-01A).*","",Prueba)

ColumnasID_Limpias <- gsub("(-01A).*|(-01B).*","",ColumnaID)
ColumnasID_Limpias

#Nueva matrix con la columna buena 
MutacionesLuad$V16 <- gsub("(-01A).*|(-01B).*","",MutacionesLuad$V16)
NombresIDLuad <- readRDS("K3_NMF_Luad.rds")

Mutaciones_k3_Luad <- MutacionesLuad[which(MutacionesLuad$V16==NombresIDLuad$Sample_ID[1]),]


#########################################################################################################################
#FECHA: 22-04-2021:

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")
MutacionesLuad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/Datos/Luad/gdc_download_20210421_083318.951054/0458c57f-316c-4a7c-9294-ccd11c97c2f9/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf", header=FALSE, comment.char="#")

MutacionesLuad$V16 <- gsub("(-01A).*|(-01B).*","",MutacionesLuad$V16)
colnames(MutacionesLuad) <- MutacionesLuad[1,]#Aqui estan todos los pacientes, y necesito tener solo a mis 61. 
saveRDS(MutacionesLuad,"TotalMutacionesLuad.rds")
GrupoK3_Luad<- readRDS("K3_NMF_Luad.rds")

PacientesLuad <- GrupoK3_Luad$Sample_ID

i <- 1
MutacionesPacienteLuad <- data.frame()
for (i in seq_along(PacientesLuad)){
  Paciente <- PacientesLuad[i]
  MutacionesPaciente_i <- MutacionesLuad[which(MutacionesLuad[,16]==Paciente),]
  MutacionesPacienteLuad_Interes <- rbind(MutacionesPaciente_i,MutacionesPacienteLuad)
  MutacionesPacienteLuad <- MutacionesPacienteLuad_Interes
}

saveRDS(MutacionesPacienteLuad_Interes,"MutacionesPacientesLuad.rds")

 #Busqueda de mutaciones: 
HugoNames <- c("TP53","STK11","ATM","KEAP1","CDKN2A")
EntrezID <- c(7157,6794,472,9817,1029)

Prueba <- MutacionesPacienteLuad_Interes[which(MutacionesPacienteLuad_Interes[,1]==HugoNames[4]),]
Prueba <- MutacionesLuad[which(MutacionesLuad$Entrez_Gene_Id==EntrezID[1]),]

Prueba <- MutacionesPacienteLuad_Interes[which(MutacionesPacienteLuad_Interes$Entrez_Gene_Id==EntrezID[1]),]


GenesMutados <- data.frame()
i <- 1

for (i in seq_along(EntrezID)){
  Mutacion <- EntrezID[i]
  Mutacion_i <- MutacionesPacienteLuad_Interes[which(MutacionesPacienteLuad_Interes$Entrez_Gene_Id==Mutacion),]
  GenesMutadosInteres <- rbind(Mutacion_i,GenesMutados)
  GenesMutados <- GenesMutadosInteres
}


