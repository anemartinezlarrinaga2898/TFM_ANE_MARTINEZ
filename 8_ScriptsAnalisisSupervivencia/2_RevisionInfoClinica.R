# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: Encontrar las diferencias que hay entre las dos bases de datos. 
################################################################################################################################
#La idea es que cargue las dos bases de datos e ir uniendo por paciente la informacion que me interesa. 
#No hace falta que lo divida entre grupos porque tengo que ver que pacientes tienen informacion diferente 

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/8_ScriptsAnalisisSupervivencia/")
library(tidyverse)
library(survival)
library(survminer)

DatosClinicos_Cbio_Luad <- readRDS("DatosClinicos/DatosClinicosLuadTotales_ConPertenencia_cBio.rds")
DatosClinicos_TGCA_Luad <- readRDS("DatosClinicos/DatosClinicosLuadTotales_ConPertenencia.rds")

DatosClinicos_Cbio_Lusc <- readRDS("DatosClinicos/DatosClinicosLuscTotales_ConPertenencia_cBio.rds")
DatosClinicos_TGCA_Lusc <- readRDS("DatosClinicos/DatosClinicosLuscTotales_ConPertenencia.rds")

################################################################################################################################
#Primero vamos a localizar aquellos pacientes que si que estan en el cBio de los que no estan. 

# ··································L U A D ······································································

CBIO <- DatosClinicos_Cbio_Luad
TGCA <- DatosClinicos_TGCA_Luad[,c(2,1,3,4)]

ElementosUnicos_CBIO <- CBIO[!CBIO$SampleID %in% TGCA$SampleID,]

# ··································L U S C ·····································································

CBIO <- DatosClinicos_Cbio_Lusc
TGCA <- DatosClinicos_TGCA_Lusc[,c(2,1,3,4)]

ElementosUnicos_CBIO <- CBIO[!CBIO$SampleID %in% TGCA$SampleID,]

################################################################################################################################

CBIO_Luad <- na.omit(DatosClinicos_Cbio_Luad)
CBIO_Lusc <- na.omit(DatosClinicos_Cbio_Lusc)

TGCA_Luad <- DatosClinicos_TGCA_Luad
TGCA_Lusc <- DatosClinicos_TGCA_Lusc

EU_CB_Luad <- CBIO_Luad[!CBIO_Luad$SampleID %in% TGCA_Luad$SampleID,]
EU_CB_Lusc<- CBIO_Lusc[!CBIO_Lusc$SampleID %in% TGCA_Lusc$SampleID,]

CBIO_Luad$times <- CBIO_Luad$times*30
CBIO_Lusc$times <- CBIO_Lusc$times*30

# ················ L U A D ··································
SV_CB_Luad <- survfit(Surv(times,status)~Grupo,data=CBIO_Luad)

ggsurvplot(SV_CB_Luad,
           title="Analisis Supervivencia Luad cBIO",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#95D1F5", "#C6A2F5","#F5E97F"))

SV_TC_Luad <- survfit(Surv(times,status)~Grupo,data=TGCA_Luad)

ggsurvplot(SV_TC_Luad,
           title="Survival Analysis",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#95D1F5", "#C6A2F5","#F5E97F"))

# ················ L U S C ··································
SV_CB_Lusc <- survfit(Surv(times,status)~Grupo,data=CBIO_Lusc)

ggsurvplot(SV_CB_Lusc,
           title="Analisis Supervivencia Lusc cBIO",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#F59D8C", "#C5F56F","#5864A8"))

SV_TC_Lusc <- survfit(Surv(times,status)~Grupo,data=TGCA_Lusc)

ggsurvplot(SV_TC_Lusc,
           title="Survival Analysis",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#F59D8C", "#C5F56F","#5864A8"))
################################################################################################################################

#En esta parte vamos a unir los datos clinicos de los pacientes con las dos bases de datos para ver en cuales no coinciden. 

CBIO_Luad <- na.omit(DatosClinicos_Cbio_Luad)
CBIO_Lusc <- na.omit(DatosClinicos_Cbio_Lusc)

CBIO_Luad$times <- CBIO_Luad$times*30
CBIO_Lusc$times <- CBIO_Lusc$times*30

TGCA_Luad <- DatosClinicos_TGCA_Luad[,c(2,1,3,4)]
TGCA_Lusc <- DatosClinicos_TGCA_Lusc[,c(2,1,3,4)]

Columnas_CBIO <- c("SampleID","Times_CBIO","Status_CBIO","Grupo_CBIO")
Columnas_TGCA <- c("SampleID","Times_TGCA","Status_TGCA","Grupo_TGCA")

# ······················ L U A D ·······························

colnames(CBIO_Luad) <- Columnas_CBIO
colnames(TGCA_Luad) <- Columnas_TGCA

UnionBD_CB_TG_Luad <- full_join(CBIO_Luad,TGCA_Luad,by="SampleID")

# ······················ L U S C ·······························

colnames(CBIO_Lusc) <- Columnas_CBIO
colnames(TGCA_Lusc) <- Columnas_TGCA

UnionBD_CB_TG_Lusc <- full_join(CBIO_Lusc,TGCA_Lusc,by="SampleID")

saveRDS(UnionBD_CB_TG_Luad,"EncontrarDiferencias/UnionBD_CB_TG_Luad.rds")
saveRDS(UnionBD_CB_TG_Lusc,"EncontrarDiferencias/UnionBD_CB_TG_Lusc.rds")

saveRDS(CBIO_Luad,"DatosClinicos/CBIO_LUAD.rds")
saveRDS(CBIO_Lusc,"DatosClinicos/CBIO_LUSC.rds")

saveRDS(TGCA_Luad,"DatosClinicos/TGCA_LUAD.rds")
saveRDS(TGCA_Lusc,"DatosClinicos/TGCA_LUSC.rds")

library("xlsx")

write.xlsx(UnionBD_CB_TG_Luad, "Excel_UNION_Luad.xlsx", sheetName = "Sheet1", 
           col.names = T, row.names = T, append = F)
write.xlsx(UnionBD_CB_TG_Lusc, "Excel_UNION_Lusc.xlsx", sheetName = "Sheet1", 
           col.names = T, row.names = T, append = F)
