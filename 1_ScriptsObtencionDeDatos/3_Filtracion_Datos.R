# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-03-2021

#OBJETIVO: Con los datos obtenidos el 21-03-2021 generar un manifiesto para descargarme los datos correspondientes. 

########################################################################################################################
library(tidyverse)

########################################################################################################################

directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/GenesOncoDrivers")

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

directoryII <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos")

manifiestoLuad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/Manifiesto_Luad_21_03_21.txt")

#Filtracion de los datos: 

manifiestoSINKRAS <- anti_join(manifiestoLuad,kras,by="id")
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
manifiestoSINALK<- anti_join(manifiestoSINRET,alk,by="id")
manifiestoLuad_Filtrado <- manifiestoSINALK #Manifiesto al que le hemos quitado los oncodrivers. 

#Tenemos los datos del paciente el case ID 
Sample_Luad <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/1_ScriptsObtencionDeDatos/RNA_SEQ_ArchivosSelecionados.tsv")

#Cambiamos los nombres para poder juntarlos: 
names(Sample_Luad) <- c("file_ID","file_name","data_category","data_type","proyect_ID","case_ID","sample_ID","sample_type")
names(manifiestoLuad_Filtrado) <- c("file_ID","file_name","md5","size","state")

ManifiestoLuad_Completo<- inner_join(Sample_Luad,manifiestoLuad_Filtrado,by=c("file_ID","file_name"))

CaseID_Duplicados <- duplicated(ManifiestoLuad_Completo$case_ID)

ManifiestoLuad_Completo_Unico <- ManifiestoLuad_Completo[!CaseID_Duplicados,]
ManifiestoLuad_Completo_Duplicados <- ManifiestoLuad_Completo[CaseID_Duplicados,]

ManifiestoDescargarLuad <- ManifiestoLuad_Completo_Unico[,c(1,2,9,10,11)]
names(ManifiestoDescargarLuad) <- c("id","filename","md5","size","state")



write_delim(manifiestoLuad_Filtrado, "ManifiestoLuad_Filtrado_No_OncoDrivers.txt", delim = "\t")
write_delim(ManifiestoLuad_Completo, "ManifiestoLuad_Completo.txt", delim = "\t")
write_delim(ManifiestoLuad_Completo_Unico, "ManifiestoLuad_Completo_Unico.txt", delim = "\t")
write_delim(ManifiestoDescargarLuad, "ManifiestoDescargarLuad.txt", delim = "\t")





