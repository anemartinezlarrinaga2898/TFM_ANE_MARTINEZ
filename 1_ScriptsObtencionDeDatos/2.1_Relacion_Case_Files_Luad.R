# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 21-03-2021

#OBJETIVO: Tras haber analizado los datos que disponia guardados y los descargados a dia de hoy, voy
#a analizar y contrarestar el workflow que segui en su dia. Primero, filtraremos ambos manifiestos el que obtuve
#en diciembre respecto al que he obtenido en marzo. Luego filtraremos por los oncodrivers, luego juntaremos
# con el nombre de los cases para ver que pacientes obtenemos.Haremos dos checks: 1) antes de haber filtrados, porque juntaremos antes y
#2) despues de haber filtrado. 

################################################################################################################################################
library(tidyverse)

##################################################################################################################################################
# MANIFIESTO DE DICIEMBRE: 551 files. 

Manifiesto551 <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/Manifiesto_551_Files.txt")
TranscriptomeProfiling_Luad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/TranscrioptomeProfiling.tsv")

#Tenemos que cambiar e nombre para que luego al poder unirlo sea igual: 
names(TranscriptomeProfiling_Luad) <- c("file_ID","file_name","data_category","data_type","proyect_ID","case_ID","sample_ID","sample_type")
names(Manifiesto551) <- c("file_ID","file_name","md5","size","state")

#Unimos los datos: 

UnionDatos_Files_Case_551 <- inner_join(TranscriptomeProfiling_Luad,Manifiesto551,by=c("file_ID","file_name"))

#Numeros de cases: 

NumeroCases_551 <- unique(UnionDatos_Files_Case_551$case_ID) 
print(length(NumeroCases_551))

#MANIFIESTO DE DICIEMBRE 547 files: 

Manifiesto547 <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/Manifiesto_547_Files.txt")
names(Manifiesto547) <- c("file_ID","file_name","md5","size","state")
UnionDatos_Files_Case_547 <- inner_join(TranscriptomeProfiling_Luad,Manifiesto547,by=c("file_ID","file_name"))

NumeroCases_547 <- unique(UnionDatos_Files_Case_547$case_ID) 
print(length(NumeroCases_547))

#En este caso es lo mismo. 

#MANIFIESTO DE MARZO 547 files

Manifiesto210321 <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/Manifiesto_Luad_21_03_21.txt")
names(Manifiesto210321) <- c("file_ID","file_name","md5","size","state")
UnionDatos_Files_Case_210321 <- inner_join(TranscriptomeProfiling_Luad,Manifiesto210321,by=c("file_ID","file_name"))

NumeroCases_210321 <- unique(UnionDatos_Files_Case_210321$case_ID) 
print(length(NumeroCases_210321))

#En los tres casos obtenemos 475 pacientes. 

#Vamos a comprobar si en los tres el numero de paciente , el ID del paciente coincide en los tres: 

length(intersect(NumeroCases_210321,NumeroCases_547)) #Son iguales

##################################################################################################################################################################################

#En esta parte vamos ha hacer la filtracion de los files que esten presentes en los manifiestos de los oncodrives.
#El objetivo es comprobar cuantos pacientes nos quedan. 

#Primero uno el manifiesto con el que he estado trabajando al transcriptomeprofiling, para tener los datos de los pacientes y ver cuantos saco. 

Luad_Filtrado <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/Luad_Filtrado.txt")
names(Luad_Filtrado) <- c("file_ID","file_name","md5","size","state")

UnionDatos_LuadFiltrado <- inner_join(TranscriptomeProfiling_Luad,Luad_Filtrado,by=c("file_ID","file_name"))
NumeroCases_LuadFiltrado <- unique(UnionDatos_LuadFiltrado$case_ID) 
print(length(NumeroCases_LuadFiltrado))

#En el manifiesto me quedo con 61 pacientes en vez de con 70. 

#Cargar los datos de los oncodrivers: 

directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/GenesOncoDrivers/")

kras <- read.table("kras.txt",header=TRUE)
egfr <- read.table("EGFR.txt",header=TRUE)
met <- read.table("MET.txt",header=TRUE)
her2 <- read.table("HER2.txt",header=TRUE)
braf <- read.table("BRAF.txt",header=TRUE)
ntrk <- read.table("ntk1.txt",header=TRUE)
pik3ca <- read.table("pik3ca.txt",header=TRUE)
mek1 <- read.table("MEK1.txt",header=TRUE)

alk <-read.table("ALK.txt",header=TRUE)
ret <-read.table("RET.txt",header=TRUE)
ros <- read.table("ROS1.txt",header=TRUE)

#Eliminacion de los genes, vamos quitando al anterior que hemos eliminado

  #Eliminacion al manifiesto de 551 files: 

manifiesto <- Manifiesto551
names(manifiesto)[1] <- "id"

manifiestoSINKRAS <- anti_join(manifiesto,kras,by="id")
manifiestoSINEGFR <- anti_join(manifiestoSINKRAS,egfr,by="id")
manifiestoSINMET <- anti_join(manifiestoSINEGFR,met,by="id")

manifiestoSINHER2 <- anti_join(manifiestoSINMET,her2,by="id")
manifiestoSINBRAF <- anti_join(manifiestoSINHER2,braf,by="id")
manifiestoSINNTRK <- anti_join(manifiestoSINBRAF,ntrk,by="id")

manifiestoSINPIKC <- anti_join(manifiestoSINNTRK,pik3ca,by="id")
manifiestoSINMEK <- anti_join(manifiestoSINPIKC,mek1,by="id")

#Eliminar fusiones: 

manifiestoSINROS1 <- anti_join(manifiestoSINMEK,ros,by="id")
manifiestoSINRET <- anti_join(manifiestoSINROS1,ret,by="id")
manifiestoSINALK_Manifiesto551 <- anti_join(manifiestoSINRET,alk,by="id")

names(manifiestoSINALK_Manifiesto551) <- c("file_ID","file_name","md5","size","state")

UnionDatos_Manifiesto551 <- inner_join(TranscriptomeProfiling_Luad,manifiestoSINALK_Manifiesto551,by=c("file_ID","file_name"))
NumeroCases_Manifiesto551 <- unique(UnionDatos_Manifiesto551$case_ID) 
print(length(NumeroCases_Manifiesto551))

  #Eliminacion del manifiesto de 547 files: 

manifiesto <- Manifiesto547
names(manifiesto)[1] <- "id"

manifiestoSINKRAS <- anti_join(manifiesto,kras,by="id")
manifiestoSINEGFR <- anti_join(manifiestoSINKRAS,egfr,by="id")
manifiestoSINMET <- anti_join(manifiestoSINEGFR,met,by="id")

manifiestoSINHER2 <- anti_join(manifiestoSINMET,her2,by="id")
manifiestoSINBRAF <- anti_join(manifiestoSINHER2,braf,by="id")
manifiestoSINNTRK <- anti_join(manifiestoSINBRAF,ntrk,by="id")

manifiestoSINPIKC <- anti_join(manifiestoSINNTRK,pik3ca,by="id")
manifiestoSINMEK <- anti_join(manifiestoSINPIKC,mek1,by="id")

#Eliminar fusiones: 

manifiestoSINROS1 <- anti_join(manifiestoSINMEK,ros,by="id")
manifiestoSINRET <- anti_join(manifiestoSINROS1,ret,by="id")
manifiestoSINALK_Manifiesto547 <- anti_join(manifiestoSINRET,alk,by="id")

names(manifiestoSINALK_Manifiesto547) <- c("file_ID","file_name","md5","size","state")

UnionDatos_Manifiesto547 <- inner_join(TranscriptomeProfiling_Luad,manifiestoSINALK_Manifiesto547,by=c("file_ID","file_name"))
NumeroCases_Manifiesto547 <- unique(UnionDatos_Manifiesto547$case_ID) 
print(length(NumeroCases_Manifiesto547))

  #Eliminacion del ultimo manifiesto descargado: 

manifiesto <- Manifiesto210321
names(manifiesto)[1] <- "id"

manifiestoSINKRAS <- anti_join(manifiesto,kras,by="id")
manifiestoSINEGFR <- anti_join(manifiestoSINKRAS,egfr,by="id")
manifiestoSINMET <- anti_join(manifiestoSINEGFR,met,by="id")

manifiestoSINHER2 <- anti_join(manifiestoSINMET,her2,by="id")
manifiestoSINBRAF <- anti_join(manifiestoSINHER2,braf,by="id")
manifiestoSINNTRK <- anti_join(manifiestoSINBRAF,ntrk,by="id")

manifiestoSINPIKC <- anti_join(manifiestoSINNTRK,pik3ca,by="id")
manifiestoSINMEK <- anti_join(manifiestoSINPIKC,mek1,by="id")

#Eliminar fusiones: 

manifiestoSINROS1 <- anti_join(manifiestoSINMEK,ros,by="id")
manifiestoSINRET <- anti_join(manifiestoSINROS1,ret,by="id")
manifiestoSINALK_Manifiesto210321 <- anti_join(manifiestoSINRET,alk,by="id")

names(manifiestoSINALK_Manifiesto210321) <- c("file_ID","file_name","md5","size","state")

UnionDatos_Manifiesto210321 <- inner_join(TranscriptomeProfiling_Luad,manifiestoSINALK_Manifiesto210321,by=c("file_ID","file_name"))
NumeroCases_Manifiesto210321 <- unique(UnionDatos_Manifiesto210321$case_ID) 
print(length(NumeroCases_Manifiesto210321))


########################################################################################################################################




