# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-04-2021

#OBJETIVO: Generacion del analisis diferencial de mis muestras, vamos ha hacerlo con el GRUPO 1 de LUAD frente
#a los grupos 2 y 3 para ver el flujo entero y asegurarme de que esta bien hecho. 


#####################################################################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(DESeq2)

#####################################################################################################################################

# 1) Generacion del coldata y de las matrices en el orden correspondiente: 

Matriz <- readRDS("Luad_Filtrado.rds")

# 1.1) En este paso generamos las matrices independientes de cada grupo.

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

# 1.2) Generamos el Col data de la comparacion 1 es que el grupo 1 es el tratado y el 2 y 3 es el de referencia. 

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

# 1.3) Generamos el objeto deseq y hacemos el analisis diferencial. 

Matriz <- Comparacion1
ColData <- Comparacion1_ColData

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL


dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

#En este paso hacemos el analisis diferencial 
DDS_1_Luad <- DESeq(dds)
DDS_1_Luad_Resultados <- results(DDS_1_Luad)

DDS_1_Luad_Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                             keys=rownames(DDS_1_Luad_Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")

ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(DDS_1_Luad_Resultados)),
                                     row.names=rownames(DDS_1_Luad_Resultados))
Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

saveRDS(Filtracion2,"ComprobacionFiltracion2_Luad.rds")

##################################################################################################################################

Filtracion2_Comprobacion <- readRDS("ComprobacionFiltracion2_Luad.rds")
Filtracion2_HuellaAntigua <- readRDS("Huella_Ensembel_Luad_1.rds")


ObjetosUnicosComprobacion <- setdiff(rownames(Filtracion2_Comprobacion),Filtracion2_HuellaAntigua) #Me dan los mismos. 



