# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-04-2021

#OBJETIVO: El objetivo de este script es generar data.frames de cada uno de los grupos con sus cuentas. Es decir las cuentas o los niveles de expresion
#del grupo 1 , las cuentas del grupo 2 y las cuentas del grupo 3.

#################################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)
#################################################################################################
# ··················· LUAD ······································

Matriz <- readRDS("Luad_Filtrado.rds")

#···················· GRUPOS ·····················

K3.1 <- readRDS("K3.1_NMF_Luad.rds") #18 pacientes
K3.2 <- readRDS("K3.2_NMF_Luad.rds") #18 pacientes
K3.3 <- readRDS("K3.3_NMF_Luad.rds") #25 pacientes

Nombres_K3.1 <- rownames(K3.1)
Nombres_K3.2 <- rownames(K3.2)
Nombres_K3.3 <- rownames(K3.3)

Grupo1 <- Matriz[,Nombres_K3.1]
Grupo2 <- Matriz[,Nombres_K3.2]
Grupo3 <- Matriz[,Nombres_K3.3]

colnames(Grupo1)==Nombres_K3.1 #Para comprobar que estan en el mismo orden. 

saveRDS(Grupo1,"Grupo3.1_NMF_Luad.rds")
saveRDS(Grupo2,"Grupo3.2_NMF_Luad.rds")
saveRDS(Grupo3,"Grupo3.3_NMF_Luad.rds")

# ··········· C O M P A R A C I O N   1  ··················

Comparacion1 <- cbind(Grupo1,Grupo2)
Comparacion1 <- cbind(Comparacion1,Grupo3)

Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="Tratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)
Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="NoTratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="NoTratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)

Comparacion1_ColData <- rbind(Grupo1_ColData,Grupo2_ColData)
Comparacion1_ColData <- rbind(Comparacion1_ColData,Grupo3_ColData)

rownames(Comparacion1_ColData)==colnames(Comparacion1)

Comparacion1_ColData$Grupo <- as.factor(Comparacion1_ColData$Grupo)
Comparacion1_ColData$Grupo

saveRDS(Comparacion1,"Comparacion1_Luad.rds")
saveRDS(Comparacion1_ColData,"Comparacion1ColData_Luad.rds")

# NO HACE FALTA HACER EL RELEVEL PORQUE AL HABERLO PUESTO EN CASTELLANO NOTRATADO VA ANTES

# ··········· C O M P A R A C I O N   2  ··················

Comparacion2 <- cbind(Grupo2,Grupo1)
Comparacion2 <- cbind(Comparacion2,Grupo3)

Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="Tratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="NoTratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)
Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="NoTratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)

Comparacion2_ColData <- rbind(Grupo2_ColData,Grupo1_ColData)
Comparacion2_ColData <- rbind(Comparacion2_ColData,Grupo3_ColData)

rownames(Comparacion2_ColData)==colnames(Comparacion2)

Comparacion2_ColData$Grupo <- as.factor(Comparacion2_ColData$Grupo)
Comparacion2_ColData$Grupo

saveRDS(Comparacion2,"Comparacion2_Luad.rds")
saveRDS(Comparacion2_ColData,"Comparacion2ColData_Luad.rds")

# ··········· C O M P A R A C I O N   3  ··················

Comparacion3 <- cbind(Grupo3,Grupo2)
Comparacion3 <- cbind(Comparacion3,Grupo1)

Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="Tratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)
Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="NoTratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="NoTratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)

Comparacion3_ColData <- rbind(Grupo3_ColData,Grupo2_ColData)
Comparacion3_ColData <- rbind(Comparacion3_ColData,Grupo1_ColData)

rownames(Comparacion3_ColData)==colnames(Comparacion3)

Comparacion3_ColData$Grupo <- as.factor(Comparacion3_ColData$Grupo)
Comparacion3_ColData$Grupo

saveRDS(Comparacion3,"Comparacion3_Luad.rds")
saveRDS(Comparacion3_ColData,"Comparacion3ColData_Luad.rds")

#################################################################################################
# ··················· LUSC ······································

Matriz <- readRDS("Lusc_Filtrado.rds")

#···················· GRUPOS ·····················

K3.1 <- readRDS("K3.1_NMF_Lusc.rds") #18 pacientes
K3.2 <- readRDS("K3.2_NMF_Lusc.rds") #18 pacientes
K3.3 <- readRDS("K3.3_NMF_Lusc.rds") #18 pacientes

Nombres_K3.1 <- rownames(K3.1)
Nombres_K3.2 <- rownames(K3.2)
Nombres_K3.3 <- rownames(K3.3)

Grupo1 <- Matriz[,Nombres_K3.1]
Grupo2 <- Matriz[,Nombres_K3.2]
Grupo3 <- Matriz[,Nombres_K3.3]

colnames(Grupo1)==Nombres_K3.1 #Para comprobar que estan en el mismo orden. 

saveRDS(Grupo1,"Grupo3.1_NMF_Lusc.rds")
saveRDS(Grupo2,"Grupo3.2_NMF_Lusc.rds")
saveRDS(Grupo3,"Grupo3.3_NMF_Lusc.rds")

# ··········· C O M P A R A C I O N   1  ··················

Comparacion1 <- cbind(Grupo1,Grupo2)
Comparacion1 <- cbind(Comparacion1,Grupo3)

Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="Tratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)
Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="NoTratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="NoTratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)

Comparacion1_ColData <- rbind(Grupo1_ColData,Grupo2_ColData)
Comparacion1_ColData <- rbind(Comparacion1_ColData,Grupo3_ColData)

rownames(Comparacion1_ColData)==colnames(Comparacion1)

Comparacion1_ColData$Grupo <- as.factor(Comparacion1_ColData$Grupo)
Comparacion1_ColData$Grupo

saveRDS(Comparacion1,"Comparacion1_Lusc.rds")
saveRDS(Comparacion1_ColData,"Comparacion1ColData_Lusc.rds")

# NO HACE FALTA HACER EL RELEVEL PORQUE AL HABERLO PUESTO EN CASTELLANO NOTRATADO VA ANTES

# ··········· C O M P A R A C I O N   2  ··················

Comparacion2 <- cbind(Grupo2,Grupo1)
Comparacion2 <- cbind(Comparacion2,Grupo3)

Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="Tratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="NoTratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)
Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="NoTratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)

Comparacion2_ColData <- rbind(Grupo2_ColData,Grupo1_ColData)
Comparacion2_ColData <- rbind(Comparacion2_ColData,Grupo3_ColData)

rownames(Comparacion2_ColData)==colnames(Comparacion2)

Comparacion2_ColData$Grupo <- as.factor(Comparacion2_ColData$Grupo)
Comparacion2_ColData$Grupo

saveRDS(Comparacion2,"Comparacion2_Lusc.rds")
saveRDS(Comparacion2_ColData,"Comparacion2ColData_Lusc.rds")

# ··········· C O M P A R A C I O N   3  ··················

Comparacion3 <- cbind(Grupo3,Grupo2)
Comparacion3 <- cbind(Comparacion3,Grupo1)

Grupo3_ColData <- data.frame(ID_Sample=colnames(Grupo3),Grupo="Tratado")
rownames(Grupo3_ColData) <- colnames(Grupo3)
Grupo2_ColData <- data.frame(ID_Sample=colnames(Grupo2),Grupo="NoTratado")
rownames(Grupo2_ColData) <- colnames(Grupo2)
Grupo1_ColData <- data.frame(ID_Sample=colnames(Grupo1),Grupo="NoTratado")
rownames(Grupo1_ColData) <- colnames(Grupo1)

Comparacion3_ColData <- rbind(Grupo3_ColData,Grupo2_ColData)
Comparacion3_ColData <- rbind(Comparacion3_ColData,Grupo1_ColData)

rownames(Comparacion3_ColData)==colnames(Comparacion3)

Comparacion3_ColData$Grupo <- as.factor(Comparacion3_ColData$Grupo)
Comparacion3_ColData$Grupo

saveRDS(Comparacion3,"Comparacion3_Lusc.rds")
saveRDS(Comparacion3_ColData,"Comparacion3ColData_Lusc.rds")

