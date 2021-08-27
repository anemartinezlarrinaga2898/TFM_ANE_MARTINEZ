# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-04-2021

#OBJETIVO: Localizacion de los genes que se encuentran mutados dentro de nuestros pacientes en el caso de luad, especificamente en el grupo que me interesa de los 12 pacientes. Los genes 
#con los que vamos a trabajar son: p53,lkb1,atm,keap1,cdkn2a

########################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")
library(tidyverse)

MutacionesLuad <- readRDS("TotalMutacionesLuad.rds")
GrupoConstante <- readRDS("GrupoConstante_Luad.rds")

MutacionesLuad$Tumor_Sample_Barcode <- str_replace_all(MutacionesLuad$Tumor_Sample_Barcode,"-",".")

PacientesLuad <- GrupoConstante$item

i <- 1
MutacionesPacienteLuad <- data.frame()
for (i in seq_along(PacientesLuad)){
  Paciente <- PacientesLuad[i]
  MutacionesPaciente_i <- MutacionesLuad[which(MutacionesLuad[,16]==Paciente),]
  MutacionesPacienteLuad_Interes <- rbind(MutacionesPaciente_i,MutacionesPacienteLuad)
  MutacionesPacienteLuad <- MutacionesPacienteLuad_Interes
}

saveRDS(MutacionesPacienteLuad_Interes,"MutacionesPacientesLuad_GrupoConstante.rds")

HugoNames <- c("TP53","STK11","ATM","KEAP1","CDKN2A")
EntrezID <- c(7157,6794,472,9817,1029)

GenesMutados <- data.frame()
i <- 1

for (i in seq_along(EntrezID)){
  Mutacion <- EntrezID[i]
  Mutacion_i <- MutacionesPacienteLuad_Interes[which(MutacionesPacienteLuad_Interes$Entrez_Gene_Id==Mutacion),]
  GenesMutadosInteres <- rbind(Mutacion_i,GenesMutados)
  GenesMutados <- GenesMutadosInteres
}

saveRDS(GenesMutadosInteres,"GenesMutadosPacientes_GrupoConstante.rds")


TP53 <- GenesMutadosInteres %>% filter(Hugo_Symbol=="TP53")
TP53 <- TP53[!duplicated(TP53$Tumor_Sample_Barcode),]

STK11 <- GenesMutadosInteres %>% filter(Hugo_Symbol=="STK11")
STK11 <- STK11[!duplicated(STK11$Tumor_Sample_Barcode),]

KEAP1 <- GenesMutadosInteres%>% filter(Hugo_Symbol=="KEAP1")
KEAP1 <- KEAP1[!duplicated(KEAP1$Tumor_Sample_Barcode),]

MatrizPerfilMutagenico <- matrix(0,nrow=3,ncol=length(PacientesLuad))
MatrizPerfilMutagenico <- as.data.frame(MatrizPerfilMutagenico)
colnames(MatrizPerfilMutagenico) <- PacientesLuad
rownames(MatrizPerfilMutagenico) <- c("TP53","STK11","KEAP1")

Mutaciones <- rbind(TP53,STK11)
Mutaciones <- rbind(Mutaciones,KEAP1)

i <- 1
j <- 1

HugoNames <- rownames(MatrizPerfilMutagenico)

for (i in seq_along(HugoNames)){
  Mutacion <- HugoNames[i]
  DataFrameMutaciones <- Mutaciones[which(Mutaciones$Hugo_Symbol==Mutacion),]
  
  for (j in seq_along(PacientesLuad)){
    Paciente <- PacientesLuad[j]
    PacientePresenteEnMutacion <- DataFrameMutaciones[which(DataFrameMutaciones$Tumor_Sample_Barcode==Paciente),]
    
    Comprobacion <- isTRUE(PacientePresenteEnMutacion$Tumor_Sample_Barcode==Paciente)
    
    if (Comprobacion==TRUE){
      Comprobacion <- 1
    } else{
      Comprobacion <- 0
    }
    
    if (Comprobacion==1){MatrizPerfilMutagenico[i,j] <- 1} 
    else {next}
  }
}


MatrizPerfilMutagenico <- MatrizPerfilMutagenico[-c(4,5),]
ComprobacionMatriz_TP53 <- sum(as.numeric(MatrizPerfilMutagenico[1,]))

saveRDS(MatrizPerfilMutagenico,"MatrizPerfilMutagenico_GrupoConstante.rds")

# ······················ M A T R I Z  P E R F I L   M U T A G E N I C O (Generacion) ·················

# 1) Generacion de la matriz con el paquete de ggplot: ···············································

MatrizPerfilMutagenico <- readRDS("MatrizPerfilMutagenico_GrupoConstante.rds")

#Hacemos la traspuesta de la matriz para tener los genes en las columnas y a los pacientes en las filas. 
Matriz <- data.frame(t(MatrizPerfilMutagenico))
#creamos una nueva columna con el nombre de los pacientes en las filas


#Eliminamos la pertenencia para poder hacer el pivot longer
Matriz_SinPertenencia <- Matriz

#Cambiamos la organizacion de la matriz
Matriz_Longer <- Matriz_SinPertenencia %>% mutate(Row = rownames(Matriz_SinPertenencia)) %>%
  pivot_longer(-Row) 
#Cambiamos el nombre de la columna del nombre del paciente para luego poder hacer el inner Join y recuperar la pertenencia
colnames(Matriz_Longer)[1] <- "SampleID"

PerfilMutagenico <- ggplot(data = Matriz_Longer, 
                           mapping = aes(x = SampleID,
                                         y = reorder(name, desc(name)),
                                         fill = as.factor(value))) +
  geom_tile() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(name = "Presencia Mutacion", values = c("grey","black"))+
  xlab(label = "Sample") +
  ylab(label="Supresor Genes")+
  ggtitle(label = "Perfil Mutagenico Grupo Constante")

PerfilMutagenico 

#2) Generacion de grafico con pheatmap: ········································


library(pheatmap)

MatrizPerfilMutagenico <- readRDS("MatrizPerfilMutagenico_GrupoConstante.rds")

#Hacemos la traspuesta de la matriz para tener los genes en las columnas y a los pacientes en las filas. 
Matriz <- data.frame(t(MatrizPerfilMutagenico))
#creamos una nueva columna con el nombre de los pacientes en las filas


#Eliminamos la pertenencia para poder hacer el pivot longer
Matriz_SinPertenencia <- Matriz

#Cambiamos la organizacion de la matriz
Matriz_Longer <- Matriz_SinPertenencia %>% mutate(Row = rownames(Matriz_SinPertenencia)) %>%
  pivot_longer(-Row) 
#Cambiamos el nombre de la columna del nombre del paciente para luego poder hacer el inner Join y recuperar la pertenencia
colnames(Matriz_Longer)[1] <- "SampleID"

MatrizDatos <- Matriz
MatrizDatos[,1] <- as.numeric(MatrizDatos[,1])
MatrizDatos[,2] <- as.numeric(MatrizDatos[,2])
MatrizDatos[,3] <- as.numeric(MatrizDatos[,3])

MatrizDatos <- as.matrix(MatrizDatos)
MatrizDatos <- t(MatrizDatos)

Pertenencia <- as.data.frame(Matriz_Longer_Pertenencia[,c(1,4)])
Pertenencia <- Pertenencia[!duplicated(Pertenencia$SampleID),]
rownames(Pertenencia) <- Pertenencia$SampleID
PertenenciaDF <- data.frame(Pertenencia[,-1])
rownames(PertenenciaDF) <- Pertenencia$SampleID
colnames(PertenenciaDF) <- "Cluster "

out <- pheatmap(MatrizDatos, 
                show_rownames=T,
                show_colnames=F,
                cluster_cols=F, 
                cluster_rows=F,
                border_color=FALSE,
                annotation_legend = T,
                main="Genes grouped by mutations")





