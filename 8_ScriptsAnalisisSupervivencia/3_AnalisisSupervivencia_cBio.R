# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: Analisis de supervivencia con los datos del cBio
###########################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/8_ScriptsAnalisisSupervivencia/")
library(tidyverse)
library(survival)
library(survminer)

###########################################################################################

# ···························· L U A D ········································

Clinica <- readRDS("ClinicalData_PacientesLuad_cBio.rds")
Clinica <- Clinica[,c(2,67,68)]

#Cambiar el nombre de las columnas para que no tengan espacios
colnames(Clinica) -> Columnas
Columnas <- str_replace_all(Columnas," ","")
colnames(Clinica) <- Columnas
colnames(Clinica)[1] <- "SampleID"

#Cambiar la columna de estado para que sea 0 o 1. 
Clinica$OverallSurvivalStatus -> Estados 
Estados <- gsub("\\:.*","",Estados)
Clinica$OverallSurvivalStatus <- Estados

#Los pacientes con sus correspondientes pertenencias: 
PacientesTotales <- readRDS("K3_NMF_Luad.rds")
PacientesTotales$nmf_subtypes <- as.numeric(PacientesTotales$nmf_subtypes)
Pacientes <- PacientesTotales$Sample_ID

i <- 1
NumeroGrupos <- as.numeric(c(1,2,3))
DatosClinicos_Pertenencia<- data.frame()
for (i in seq_along(NumeroGrupos)){
  NumeroGrupo <- NumeroGrupos[i]
  
  PacientesGrupo <- PacientesTotales %>% filter(nmf_subtypes==NumeroGrupo)
  
  ID_Paciente <- PacientesGrupo$Sample_ID
  
  j <- 1
  DatosClinicosGrupos <- data.frame()
  for (j in seq_along(ID_Paciente)){
    Paciente <- ID_Paciente[j]
    Clinico <- Clinica[which(Clinica$SampleID==Paciente),]
    DatosClinicosGrupos <- rbind(DatosClinicosGrupos,Clinico)
  }
  DatosClinicosGrupos <- DatosClinicosGrupos %>% mutate(Grupo=NumeroGrupo)
  DatosClinicos_Pertenencia <- rbind(DatosClinicos_Pertenencia,DatosClinicosGrupos)
}
colnames(DatosClinicos_Pertenencia)[c(2,3)] <- c("times","status")
DatosClinicos_Pertenencia$status <- as.numeric(DatosClinicos_Pertenencia$status)

saveRDS(DatosClinicos_Pertenencia,"DatosClinicosLuadTotales_ConPertenencia_cBio.rds")
saveRDS(Clinica,"DatosLuad_cBio_Resumen.rds")

# A n a l i s i s  s u p e r v i v e n c i a  e n  Meses: 

sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
summary(sfit,times=seq(0,365*5*365)) # for the #rst 5 years survival rates.
ggsurvplot(sfit,
           title="Analisis Supervivencia Luad",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#95D1F5", "#C6A2F5","#F5E97F"))

# A n a l i s i s  s u p e r v i v e n c i a  e n  Dias:

DatosClinicos_Pertenencia$times <- DatosClinicos_Pertenencia$times*30
sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
ggsurvplot(sfit)

#Test estadistico de si realmente hay diferencia entre los tres grupos 

surv_diff <- survdiff(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
surv_diff
# p.value =  0.6 no hay diferencia entre los tres grupos

###########################################################################################

# ···························· L U S C ········································

Clinica <- readRDS("ClinicalData_PacientesLusc_cBio.rds")
Clinica <- Clinica[,c(2,67,68)]

colnames(Clinica) -> Columnas
Columnas <- str_replace_all(Columnas," ","")
colnames(Clinica) <- Columnas
colnames(Clinica)[1] <- "SampleID"

#Cambiar la columna de estado para que sea 0 o 1. 
Clinica$OverallSurvivalStatus -> Estados 
Estados <- gsub("\\:.*","",Estados)
Clinica$OverallSurvivalStatus <- Estados

#Los pacientes con sus correspondientes pertenencias: 
PacientesTotales <- readRDS("K3_NMF_Lusc.rds")
PacientesTotales$nmf_subtypes <- as.numeric(PacientesTotales$nmf_subtypes)
Pacientes <- PacientesTotales$Sample_ID

i <- 1
NumeroGrupos <- as.numeric(c(1,2,3))
DatosClinicos_Pertenencia<- data.frame()
for (i in seq_along(NumeroGrupos)){
  NumeroGrupo <- NumeroGrupos[i]
  
  PacientesGrupo <- PacientesTotales %>% filter(nmf_subtypes==NumeroGrupo)
  
  ID_Paciente <- PacientesGrupo$Sample_ID
  
  j <- 1
  DatosClinicosGrupos <- data.frame()
  for (j in seq_along(ID_Paciente)){
    Paciente <- ID_Paciente[j]
    Clinico <- Clinica[which(Clinica$SampleID==Paciente),]
    DatosClinicosGrupos <- rbind(DatosClinicosGrupos,Clinico)
  }
  DatosClinicosGrupos <- DatosClinicosGrupos %>% mutate(Grupo=NumeroGrupo)
  DatosClinicos_Pertenencia <- rbind(DatosClinicos_Pertenencia,DatosClinicosGrupos)
}
colnames(DatosClinicos_Pertenencia)[c(2,3)] <- c("times","status")
DatosClinicos_Pertenencia$status <- as.numeric(DatosClinicos_Pertenencia$status)

saveRDS(DatosClinicos_Pertenencia,"DatosClinicosLuscTotales_ConPertenencia_cBio.rds")
saveRDS(Clinica,"DatosLusc_cBio_Resumen.rds")

# A n a l i s i s  s u p e r v i v e n c i a  e n  Meses: 

sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
summary(sfit,times=seq(0,365*5*365)) # for the #rst 5 years survival rates.
ggsurvplot(sfit)

ggsurvplot(sfit,
           title="Analisis Supervivencia Lusc",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#F59D8C", "#C5F56F","#5864A8"))

# A n a l i s i s  s u p e r v i v e n c i a  e n  Dias:

DatosClinicos_Pertenencia$times <- DatosClinicos_Pertenencia$times*30
sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
ggsurvplot(sfit)

surv_diff <- survdiff(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia)
surv_diff

# P.value es 0.006 es significativo!