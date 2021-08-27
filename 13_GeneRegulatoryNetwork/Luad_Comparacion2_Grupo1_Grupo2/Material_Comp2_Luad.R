# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-06-2021

#OBJETIVO: Generacion del material para correr el Trare

###########################################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/13_GeneRegulatoryNetwork/")
library(tidyverse)

###########################################################################################################################################################

# Generar la matriz correspondiente con mis pacientes, del grupo 1 y del grupo 2. 

Lusc <- readRDS("MatrizLuad_Genes_HugoSymbol.rds")

Grupo1 <- readRDS("K3.1_NMF_Luad.rds")
Grupo2 <- readRDS("K3.2_NMF_Luad.rds")

Grupo1 <- Grupo1$Sample_ID
Grupo2 <- Grupo2$Sample_ID

Pacientes <- list(Grupo1,Grupo2)
names(Pacientes) <- c("Pacientes Grupo 1","Pacientes Grupo 2")

#Obtencion de las cuentas de todos los genes de los pacientes de interes: 

MatrizCuentasFinal1 <- matrix(0,nrow=nrow(Lusc),ncol=length(Grupo1))
rownames(MatrizCuentasFinal1) <- rownames(Lusc)
colnames(MatrizCuentasFinal1) <- Grupo1
MatrizCuentasFinal2 <- matrix(0,nrow=nrow(Lusc),ncol=length(Grupo2))
rownames(MatrizCuentasFinal2) <- rownames(Lusc)
colnames(MatrizCuentasFinal2) <- Grupo2
Matrices <- list(MatrizCuentasFinal1,MatrizCuentasFinal2)
names(Matrices) <- c("Pacientes Grupo 1","Pacientes Grupo 2")

i <- 1

for (i in seq_along(Pacientes)){
  Grupo <- Pacientes[[i]]
  
  Matriz <- Matrices[[i]]
  
  j <- 1
  
  for (j in seq_along(Grupo)){
    Paciente <- Grupo[j]
    
    CuentasPaciente <- Lusc[,which(colnames(Lusc)==Paciente)]
    
    Matriz[,j] <- CuentasPaciente
  }
  
  Matrices[[i]] <- Matriz
  
}

Matriz1 <- Matrices[[1]]
Matriz2 <- Matrices[[2]]

MatrizCuentasLusc <- cbind(Matriz1,Matriz2) #Esta es la matriz final, donde tengo a mis 120 pacientes. 

saveRDS(MatrizCuentasLusc,"Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_SinNormalizar.rds")

###########################################################################################################################################################

# Ahora tenemos que noramlizar, es decir aquella normalizacion que estabilize mejor la varianza. 
# Son tres metodos de estabilizacion de la varianza 
# Normalization of low-count-filtered raw count matrixes

Matriz <- readRDS("Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_SinNormalizar.rds") #Los samples ( pacientes ) tienen que estar en las columnas. 
dim(Matriz)

library(DESeq2)
library(vsn) #Para hacer los graficos coger la funcion vsn::meansSDPlot
library(hexbin) #Tambien es necesario para 
library(egg)
#Hacemos las normalizacion con la VST que es una funcion que se encuentra en el paquete de Deseq2
Normalizacion_Parametrica<-varianceStabilizingTransformation((Matriz),fitType = "parametric")
Normalizacion_Media<-varianceStabilizingTransformation((Matriz),fitType = "mean")
Normalizacion_Local<-varianceStabilizingTransformation((Matriz),fitType = "local")

pdf("normalization types.pdf")
a<-vsn::meanSdPlot(Normalizacion_Parametrica)$gg + ggtitle("Parametric Dispersion Fitting")
b<-vsn::meanSdPlot(Normalizacion_Local)$gg + ggtitle("Local Dispersion Fitting")
c<-vsn::meanSdPlot(Normalizacion_Media)$gg + ggtitle("Mean Dispersion Fitting")
egg::ggarrange(a,b,c,nrow=2,ncol=2)

dev.off()

#Cogemos la segunda normalizacion, la LOCAL. 
saveRDS(Normalizacion_Local,"Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_Normalizada.rds")


############################################################################################################

##################################################################################################################################

Matriz <- readRDS("Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_Normalizar.rds")

# Grupo 1: NR
# Grupo 2: R

DF_FenotipoPacientes <- data.frame(Samples=colnames(Matriz),Class=rep(0,ncol(Matriz)))
DF_FenotipoPacientes[1:18,2] <- "NR"
DF_FenotipoPacientes[19:36,2] <- "R"

DF_FenotipoPacientes$Class <- as.factor(DF_FenotipoPacientes$Class)
DF_FenotipoPacientes$Class

table(DF_FenotipoPacientes$Class)


saveRDS(DF_FenotipoPacientes,"Luad_Comparacion2_Grupo1_Grupo2/DF_Fenotipo_Comp2_Luad.rds")

# OBJETOS REWIRING: 

Matriz <- readRDS("Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_Normalizada.rds")
write.table(Matriz,"Luad_Comparacion2_Grupo1_Grupo2/Matriz_Comp2_Luad_Normalizada.tsv",sep="\t")

Fenotipo <- readRDS("Luad_Comparacion2_Grupo1_Grupo2/DF_Fenotipo_Comp2_Luad.rds")
write.table(Fenotipo,"Luad_Comparacion2_Grupo1_Grupo2/Fenotipo_Comp2_Luad.tsv",sep="\t",row.names = F)



