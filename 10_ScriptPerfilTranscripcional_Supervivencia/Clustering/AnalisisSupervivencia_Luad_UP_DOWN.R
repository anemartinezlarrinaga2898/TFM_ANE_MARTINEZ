# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-05-2021

#OBJETIVO:Script para realizar el analisis de supervivencia con los nuevos grupos que hemos generado. LUAD

##############################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/10_ScriptPerfilTranscripcional_Supervivencia/")
library(tidyverse)
library(survival)
library(survminer)

#########################################################################################################

#······································ L U A D ···············································

Clinica <- readRDS("Clustering/Clinica_Total_Luad.rds")

# C O M P A R A C I O N  1 : ·································································

UP_1 <- Clinica[[1]]
DOWN_1 <- Clinica[[2]]

rownames(UP_1) <- UP_1[,1]
rownames(DOWN_1) <- UP_1[,1]

UP_1 <- UP_1[,-c(1,2)]
DOWN_1 <- DOWN_1[,-c(1,2)]

UP_1$Grupo <- ifelse(UP_1$Grupo=="UP",1,0)
DOWN_1$Grupo <- ifelse(DOWN_1$Grupo=="UP",1,0)

UP_1$Grupo <- as.factor(UP_1$Grupo)
DOWN_1$Grupo <- as.factor(DOWN_1$Grupo)

colnames(UP_1)[c(2,3)] <- c("times","status")
colnames(DOWN_1)[c(2,3)] <- c("times","status")

UP_1$status <- as.numeric(UP_1$status)
DOWN_1$status <- as.numeric(DOWN_1$status)

Survival_1_UP <- survfit(Surv(times,status)~Grupo,data=UP_1)
Survival_1_DOWN <- survfit(Surv(times,status)~Grupo,data=DOWN_1)

Plot_Survival_1_UP <- ggsurvplot(Survival_1_UP,
           title="Group 1 UP Regulated Genes",
           legend.title = "Clustering",
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           palette = c("#95D1F5", "#C6A2F5"))
Plot_Survival_1_UP

Plot_Survival_1_DOWN <- (ggsurvplot(Survival_1_DOWN,
                                 title="Group 1 DOWN Regulated Genes",
                                 legend.title = "Clustering",
                                 pval = TRUE,
                                 risk.table = TRUE, 
                                 risk.table.y.text.col = TRUE,
                                 risk.table.col = "strata",
                                 palette = c("#95D1F5", "#C6A2F5")))
Plot_Survival_1_DOWN

# C O M P A R A C I O N  2 : ·································································

UP_2 <- Clinica[[3]]
DOWN_2 <- Clinica[[4]]

rownames(UP_2) <- UP_2[,1]
rownames(DOWN_2) <- UP_2[,1]

UP_2 <- UP_2[,-c(1,2)]
DOWN_2 <- DOWN_2[,-c(1,2)]

UP_2$Grupo <- ifelse(UP_2$Grupo=="UP",1,0)
DOWN_2$Grupo <- ifelse(DOWN_2$Grupo=="UP",1,0)

UP_2$Grupo <- as.factor(UP_2$Grupo)
DOWN_2$Grupo <- as.factor(DOWN_2$Grupo)

colnames(UP_2)[c(2,3)] <- c("times","status")
colnames(DOWN_2)[c(2,3)] <- c("times","status")

UP_2$status <- as.numeric(UP_2$status)
DOWN_2$status <- as.numeric(DOWN_2$status)

Survival_2_UP <- survfit(Surv(times,status)~Grupo,data=UP_2)
Survival_2_DOWN <- survfit(Surv(times,status)~Grupo,data=DOWN_2)

Plot_Survival_2_UP <- ggsurvplot(Survival_2_UP,
                                 title="Group 2 UP Regulated Genes",
                                 legend.title = "Clustering",
                                 pval = TRUE,
                                 risk.table = TRUE, 
                                 risk.table.y.text.col = TRUE,
                                 risk.table.col = "strata",
                                 palette = c("#9AC4F8", "#E36588"))
Plot_Survival_2_UP

Plot_Survival_2_DOWN <- (ggsurvplot(Survival_2_DOWN,
                                    title="Group 2 DOWN Regulated Genes",
                                    legend.title = "Clustering",
                                    pval = TRUE,
                                    risk.table = TRUE, 
                                    risk.table.y.text.col = TRUE,
                                    risk.table.col = "strata",
                                    palette = c("#9AC4F8", "#E36588")))
Plot_Survival_2_DOWN

# C O M P A R A C I O N  3 : ·································································

UP_3 <- Clinica[[5]]
DOWN_3 <- Clinica[[6]]

rownames(UP_3) <- UP_3[,1]
rownames(DOWN_3) <- UP_3[,1]

UP_3 <- UP_3[,-c(1,2)]
DOWN_3 <- DOWN_3[,-c(1,2)]

UP_3$Grupo <- ifelse(UP_3$Grupo=="UP",1,0)
DOWN_3$Grupo <- ifelse(DOWN_3$Grupo=="UP",1,0)

UP_3$Grupo <- as.factor(UP_3$Grupo)
DOWN_3$Grupo <- as.factor(DOWN_3$Grupo)

colnames(UP_3)[c(2,3)] <- c("times","status")
colnames(DOWN_3)[c(2,3)] <- c("times","status")

UP_3$status <- as.numeric(UP_3$status)
DOWN_3$status <- as.numeric(DOWN_3$status)

Survival_3_UP <- survfit(Surv(times,status)~Grupo,data=UP_3)
Survival_3_DOWN <- survfit(Surv(times,status)~Grupo,data=DOWN_3)

Plot_Survival_3_UP <- ggsurvplot(Survival_3_UP,
                                 title="Group 3 UP Regulated Genes",
                                 legend.title = "Clustering",
                                 pval = TRUE,
                                 risk.table = TRUE, 
                                 risk.table.y.text.col = TRUE,
                                 risk.table.col = "strata",
                                 palette = c("#8CB369", "#5B8E7D"))
Plot_Survival_3_UP

Plot_Survival_3_DOWN <- (ggsurvplot(Survival_3_DOWN,
                                    title="Group 3 DOWN Regulated Genes",
                                    legend.title = "Clustering",
                                    pval = TRUE,
                                    risk.table = TRUE, 
                                    risk.table.y.text.col = TRUE,
                                    risk.table.col = "strata",
                                    palette = c("#8CB369", "#5B8E7D")))
Plot_Survival_3_DOWN


