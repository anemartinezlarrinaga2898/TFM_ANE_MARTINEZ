# Descargados los archivos de los diferentes genes y he hecho la extraccion de los ID que nome interesan.
# Por otro lado cargamos el manifiesto LUAD para cargalos a la base de datos

library(tidyverse)

a <- data.frame(x=seq(1,10,1))
b <- data.frame(x=c(1,5,7,8,1,11,12,15,4,10))

newDataFrame <- anti_join(a,b)

setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/DATA")

manifiesto <- read.table("manifiestoORIGINAL.txt",header=TRUE)

kras <- read.table("kras.txt",header=TRUE)
egfr <- read.table("EGFR.txt",header=TRUE)
met <- read.table("MET.txt",header=TRUE)
her2 <- read.table("HER2.txt",header=TRUE)
braf <- read.table("BRAF.txt",header=TRUE)
ntrk <- read.table("ntk1.txt",header=TRUE)
pik3ca <- read.table("pik3ca_prueba.txt",header=TRUE)
mek1 <- read.table("MEK1_2.txt",header=TRUE)

manifiestoSINKRAS <- anti_join(manifiesto,kras,by="id")
manifiestoSINEGFR <- anti_join(manifiestoSINKRAS,egfr,by="id")
manifiestoSINMET <- anti_join(manifiestoSINEGFR,met,by="id")

manifiestoSINHER2 <- anti_join(manifiestoSINMET,her2,by="id")
manifiestoSINBRAF <- anti_join(manifiestoSINHER2,braf,by="id")
manifiestoSINNTRK <- anti_join(manifiestoSINBRAF,ntrk,by="id")

manifiestoSINPIKC <- anti_join(manifiestoSINNTRK,pik3ca,by="id")
manifiestoSINMEK <- anti_join(manifiestoSINPIKC,mek1,by="id")




# COMPROBACION: El nombre de los archivos el nombre del gen seguido de un _2:
kras2 <- read.table("KRAS_2.txt",header=TRUE)
egfr2 <- read.table("EGFR_2.txt",header=TRUE)
met2 <- read.table("MET_2.txt",header=TRUE)
her2_2 <- read.table("HER2_2.txt",header=TRUE)
braf2 <- read.table("BRAF_2.txt",header=TRUE)
ntrk2 <- read.table("NTRK_2.txt",header=TRUE)
pik3ca2 <- read.table("PIK3CA_2.txt",header=TRUE)
mek12 <- read.table("MEK1_2.txt",header=TRUE)

# Fusiones 
alk <-read.table("ALK_2.txt",header=TRUE)
ret <-read.table("RET_2.txt",header=TRUE)
ros <- read.table("ROS_1.txt",header=TRUE)

#Eliminacion de los genes, vamos quitando al anterior que hemos eliminado

manifiestoSINKRAS <- anti_join(manifiesto,kras2,by="id")
manifiestoSINEGFR <- anti_join(manifiestoSINKRAS,egfr2,by="id")
manifiestoSINMET <- anti_join(manifiestoSINEGFR,met2,by="id")

manifiestoSINHER2 <- anti_join(manifiestoSINMET,her2_2,by="id")
manifiestoSINBRAF <- anti_join(manifiestoSINHER2,braf2,by="id")
manifiestoSINNTRK <- anti_join(manifiestoSINBRAF,ntrk2,by="id")

manifiestoSINPIKC <- anti_join(manifiestoSINNTRK,pik3ca2,by="id")
manifiestoSINMEK <- anti_join(manifiestoSINPIKC,mek12,by="id")

#Eliminar fusiones: 

manifiestoSINROS1 <- anti_join(manifiestoSINMEK,ros,by="id")
manifiestoSINRET <- anti_join(manifiestoSINROS1,ret,by="id")
manifiestoSINALK <- anti_join(manifiestoSINRET,alk,by="id")

LUAD <- manifiestoSINALK
write_delim(iris, "./datos/iris_3.txt", delim = "\t")
write_delim(LUAD, "LUAD.txt", delim = "\t")




