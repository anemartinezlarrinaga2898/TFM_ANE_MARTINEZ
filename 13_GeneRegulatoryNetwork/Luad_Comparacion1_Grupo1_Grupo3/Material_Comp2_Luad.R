# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-06-2021

#OBJETIVO: Generacion del material para correr el Trare

###########################################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/13_GeneRegulatoryNetwork/")
library(tidyverse)

###########################################################################################################################################################

# Generar la matriz correspondiente con mis pacientes, del grupo 1 y del grupo 3. 

Lusc <- readRDS("MatrizLuad_Genes_HugoSymbol.rds")

Grupo1 <- readRDS("K3.1_NMF_Luad.rds")
Grupo3 <- readRDS("K3.3_NMF_Luad.rds")

Grupo1 <- Grupo1$Sample_ID
Grupo3 <- Grupo3$Sample_ID

Pacientes <- list(Grupo1,Grupo3)
names(Pacientes) <- c("Pacientes Grupo 1","Pacientes Grupo 3")

#Obtencion de las cuentas de todos los genes de los pacientes de interes: 

MatrizCuentasFinal1 <- matrix(0,nrow=nrow(Lusc),ncol=length(Grupo1))
rownames(MatrizCuentasFinal1) <- rownames(Lusc)
colnames(MatrizCuentasFinal1) <- Grupo1
MatrizCuentasFinal3 <- matrix(0,nrow=nrow(Lusc),ncol=length(Grupo3))
rownames(MatrizCuentasFinal3) <- rownames(Lusc)
colnames(MatrizCuentasFinal3) <- Grupo3
Matrices <- list(MatrizCuentasFinal1,MatrizCuentasFinal3)
names(Matrices) <- c("Pacientes Grupo 1","Pacientes Grupo 3")

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
Matriz3 <- Matrices[[2]]

MatrizCuentasLusc <- cbind(Matriz1,Matriz3) #Esta es la matriz final, donde tengo a mis 120 pacientes. 

saveRDS(MatrizCuentasLusc,"Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_SinNormalizar.rds")

###########################################################################################################################################################

# Ahora tenemos que noramlizar, es decir aquella normalizacion que estabilize mejor la varianza. 
# Son tres metodos de estabilizacion de la varianza 
# Normalization of low-count-filtered raw count matrixes

Matriz <- readRDS("Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_SinNormalizar.rds") #Los samples ( pacientes ) tienen que estar en las columnas. 
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
saveRDS(Normalizacion_Local,"Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_Normalizada.rds")


############################################################################################################

##################################################################################################################################

Matriz <- readRDS("Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_SinNormalizar.rds")

# Grupo 1: NR
# Grupo 3: R

DF_FenotipoPacientes <- data.frame(Samples=colnames(Matriz),Class=rep(0,ncol(Matriz)))
DF_FenotipoPacientes[1:18,2] <- "NR"
DF_FenotipoPacientes[(19:43),2] <- "R"

DF_FenotipoPacientes$Class <- as.factor(DF_FenotipoPacientes$Class)
DF_FenotipoPacientes$Class

table(DF_FenotipoPacientes$Class)

saveRDS(DF_FenotipoPacientes,"Luad_Comparacion1_Grupo1_Grupo3/DF_Fenotipo_Comp1_Luad.rds")

# OBJETOS REWIRING: 

Matriz <- readRDS("Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_Normalizada.rds")
write.table(Matriz,"Luad_Comparacion1_Grupo1_Grupo3/Matriz_Comp1_Luad_Normalizada.tsv",sep="\t")

Fenotipo <- readRDS("Luad_Comparacion1_Grupo1_Grupo3/DF_Fenotipo_Comp1_Luad.rds")
write.table(Fenotipo,"Luad_Comparacion1_Grupo1_Grupo3/Fenotipo_Comp1_Luad.tsv",sep="\t",row.names = F)

FT <- readRDS("Luad_Comparacion1_Grupo1_Grupo3/DF_FT_Luad.rds")
write.table(FT,"Luad_Comparacion1_Grupo1_Grupo3/DF_FT_Comp1_Luad.tsv",sep="\t",row.names = F)


