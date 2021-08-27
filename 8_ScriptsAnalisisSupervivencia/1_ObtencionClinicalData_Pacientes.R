# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 04-05-2021

#OBJETIVO: Con los datos que tengo de la TGCA enlazar mis pacientes con sus datos. Luego hacer lo mismo con los datos que me descargo del cBIO. 

##################################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/8_ScriptsAnalisisSupervivencia/")

# ································· L U A D -- Datos TGCA -- ················································

library(readr)
ClinicalTotal <- read_delim("Clinical_Luad_TGCA.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

Grupos <- readRDS("K3_NMF_Luad.rds")
Pacientes <- Grupos[,1]

Clinical_Pacientes <- data.frame()
i <- 1

for (i in seq_along(Pacientes)){
  Paciente <- Pacientes[i]
  Paciente_Clinical <- ClinicalTotal[which(ClinicalTotal$case_submitter_id==Paciente),]
  Clinical_Pacientes <- rbind(Clinical_Pacientes,Paciente_Clinical)
}

saveRDS(Clinical_Pacientes,"ClinicalData_PacientesLuad_TGCA.rds")

# ································· L U S C -- Datos TGCA -- ················································

library(readr)
ClinicalTotal <- read_delim("Clinical_Lusc_TGCA.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

Grupos <- readRDS("K3_NMF_Lusc.rds")
Pacientes <- Grupos[,1]

Clinical_Pacientes <- data.frame()
i <- 1

for (i in seq_along(Pacientes)){
  Paciente <- Pacientes[i]
  Paciente_Clinical <- ClinicalTotal[which(ClinicalTotal$case_submitter_id==Paciente),]
  Clinical_Pacientes <- rbind(Clinical_Pacientes,Paciente_Clinical)
}

saveRDS(Clinical_Pacientes,"ClinicalData_PacientesLusc_TGCA.rds")

#########################################################################################################

#En esta parte vamos a trabajar con los datos disponibles en el cBio portal. 

# ································· L U A D -- Datos cBIO -- ················································

ClinicalTotal <- read_delim("Clinical_Luad_cBio.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

Grupos <- readRDS("K3_NMF_Luad.rds")
Pacientes <- Grupos[,1]

Clinical_Pacientes <- data.frame()
i <- 1

for (i in seq_along(Pacientes)){
  Paciente <- Pacientes[i]
  Paciente_Clinical <- ClinicalTotal[which(ClinicalTotal$`Patient ID`==Paciente),]
  Clinical_Pacientes <- rbind(Clinical_Pacientes,Paciente_Clinical)
}

saveRDS(Clinical_Pacientes,"ClinicalData_PacientesLuad_cBio.rds")

# ································· L U S C -- Datos cBIO -- ················································

ClinicalTotal <- read_delim("Clinical_Lusc_cBio.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

Grupos <- readRDS("K3_NMF_Lusc.rds")
Pacientes <- Grupos[,1]

Clinical_Pacientes <- data.frame()
i <- 1

for (i in seq_along(Pacientes)){
  Paciente <- Pacientes[i]
  Paciente_Clinical <- ClinicalTotal[which(ClinicalTotal$`Patient ID`==Paciente),]
  Clinical_Pacientes <- rbind(Clinical_Pacientes,Paciente_Clinical)
}

saveRDS(Clinical_Pacientes,"ClinicalData_PacientesLusc_cBio.rds")

