# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: Analisis de supervivencia con los datos del RTCGA
###########################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/8_ScriptsAnalisisSupervivencia/")
library(tidyverse)
library(survival)
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)
library("survminer")

###########################################################################################

# ···························· L U A D ········································

DatosClinicos <- survivalTCGA(LUAD.clinical)
colnames(DatosClinicos)[c(2,3)] <- c("SampleID","status")

PacientesLuadTotales <- readRDS("K3_NMF_Luad.rds")
PacientesLuadTotales$nmf_subtypes <- as.numeric(PacientesLuadTotales$nmf_subtypes)
PacientesLuad <- PacientesLuadTotales$Sample_ID

DatosClinicosLuad <- data.frame() #Hay tres pacientes que no tienen informacion clinica
i <- 1

for (i in seq_along(PacientesLuad)){
  Paciente <- PacientesLuad[i]
  Clinico <- DatosClinicos[which(DatosClinicos$SampleID==Paciente),]
  DatosClinicosLuad <- rbind(DatosClinicosLuad,Clinico)
}

i <- 1
NumeroGrupos <- as.numeric(c(1,2,3))
DatosClinicos_Pertenencia_Luad <- data.frame()
for (i in seq_along(NumeroGrupos)){
  NumeroGrupo <- NumeroGrupos[i]
  
  PacientesGrupo <- PacientesLuadTotales %>% filter(nmf_subtypes==NumeroGrupo)
  
  ID_Paciente <- PacientesGrupo$Sample_ID
  
  j <- 1
  DatosClinicosGrupos <- data.frame()
  for (j in seq_along(ID_Paciente)){
    Paciente <- ID_Paciente[j]
    Clinico <- DatosClinicosLuad[which(DatosClinicosLuad$SampleID==Paciente),]
    DatosClinicosGrupos <- rbind(DatosClinicosGrupos,Clinico)
  }
  DatosClinicosGrupos <- DatosClinicosGrupos %>% mutate(Grupo=NumeroGrupo)
  DatosClinicos_Pertenencia_Luad <- rbind(DatosClinicos_Pertenencia_Luad,DatosClinicosGrupos)
}

# El data frame de interes es el de DatosClinicos_Pertenencia_Luad
xtabs(~Grupo+status, data=DatosClinicos_Pertenencia_Luad) %>% addmargins()
sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia_Luad)
summary(sfit,times=seq(0,365*5*365)) # for the #rst 5 years survival rates.
ggsurvplot(sfit,
           title="Analisis Supervivencia Luad",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#95D1F5", "#C6A2F5","#F5E97F"))

surv_diff <- survdiff(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia_Luad)
surv_diff

#p.value =1

saveRDS(DatosClinicosLuad,"DatosClinicosLuadTotales_SinPertenencia.rds")
saveRDS(DatosClinicos_Pertenencia_Luad,"DatosClinicosLuadTotales_ConPertenencia.rds")

###########################################################################################

# ···························· L U S C ········································

DatosClinicos <- survivalTCGA(LUSC.clinical)
colnames(DatosClinicos)[c(2,3)] <- c("SampleID","status")

PacientesLuscTotales <- readRDS("K3_NMF_Lusc.rds")
PacientesLuscTotales$nmf_subtypes <- as.numeric(PacientesLuscTotales$nmf_subtypes)
PacientesLusc <- PacientesLuscTotales$Sample_ID

DatosClinicosLusc <- data.frame() #Hay tres pacientes que no tienen informacion clinica
i <- 1

for (i in seq_along(PacientesLusc)){
  Paciente <- PacientesLusc[i]
  Clinico <- DatosClinicos[which(DatosClinicos$SampleID==Paciente),]
  DatosClinicosLusc <- rbind(DatosClinicosLusc,Clinico)
}

i <- 1
NumeroGrupos <- as.numeric(c(1,2,3))
DatosClinicos_Pertenencia_Lusc <- data.frame()
for (i in seq_along(NumeroGrupos)){
  NumeroGrupo <- NumeroGrupos[i]
  
  PacientesGrupo <- PacientesLuscTotales %>% filter(nmf_subtypes==NumeroGrupo)
  
  ID_Paciente <- PacientesGrupo$Sample_ID
  
  j <- 1
  DatosClinicosGrupos <- data.frame()
  for (j in seq_along(ID_Paciente)){
    Paciente <- ID_Paciente[j]
    Clinico <- DatosClinicosLusc[which(DatosClinicosLusc$SampleID==Paciente),]
    DatosClinicosGrupos <- rbind(DatosClinicosGrupos,Clinico)
  }
  DatosClinicosGrupos <- DatosClinicosGrupos %>% mutate(Grupo=NumeroGrupo)
  DatosClinicos_Pertenencia_Lusc <- rbind(DatosClinicos_Pertenencia_Lusc,DatosClinicosGrupos)
}

saveRDS(DatosClinicosLusc,"DatosClinicosLuscTotales_SinPertenencia.rds")
saveRDS(DatosClinicos_Pertenencia_Lusc,"DatosClinicosLuscTotales_ConPertenencia.rds")

# El data frame de interes es el de DatosClinicos_Pertenencia_Lusc
xtabs(~Grupo+status, data=DatosClinicos_Pertenencia_Lusc) %>% addmargins()
sfit <- survfit(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia_Lusc)
summary(sfit,times=seq(0,365*5*365)) # for the #rst 5 years survival rates.
ggsurvplot(sfit,
           title="Analisis Supervivencia Lusc",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#F59D8C", "#C5F56F","#5864A8"))


surv_diff <- survdiff(Surv(times,status)~Grupo,data=DatosClinicos_Pertenencia_Lusc)
surv_diff

#p.value es de 0.07
