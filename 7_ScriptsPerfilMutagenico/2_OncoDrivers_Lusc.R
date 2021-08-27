# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: matriz de las mutaciones que puede tener LUSC en los diferentes genes ONCO-DRIVERS. 

##########################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")
library(tidyverse)

Mutaciones <- readRDS("TotalMutacionesLusc.rds")
MutacionesLusc <- readRDS("MutacionesPacientesLusc.rds")
GrupoK3_Lusc<- readRDS("K3_NMF_Lusc.rds")
GrupoK3_Lusc <- GrupoK3_Lusc[order(GrupoK3_Lusc$nmf_subtypes),]
GrupoK3_Lusc <- GrupoK3_Lusc[,-3]



HugoNames <- c("PIK3CA","PIK3CB","NOTCH1","NF1","MYCN","EGFR","FGF6","FGF3","FGF10","SOX2","MYC")
EntrezID <- c(5290,5291,4851,4763,4613,1956,2251,2248,2255,6657,4609)

GenesMutados <- data.frame()
i <- 1

for (i in seq_along(EntrezID)){
  Mutacion <- EntrezID[i]
  Mutacion_i <- MutacionesLusc[which(MutacionesLusc$Entrez_Gene_Id==Mutacion),]
  GenesMutadosInteres <- rbind(Mutacion_i,GenesMutados)
  GenesMutados <- GenesMutadosInteres
}

colnames(GenesMutadosInteres)[16] <- "Sample_ID"
GenesMutadosInteres_PertenenciaCluster <- inner_join(GrupoK3_Lusc,GenesMutadosInteres,by="Sample_ID")

PIK3CA <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="PIK3CA")
PIK3CA <- PIK3CA[order(PIK3CA$nmf_subtypes),]
PIK3CA <- PIK3CA[!duplicated(PIK3CA$Sample_ID),]

PIK3CB <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="PIK3CB")
PIK3CB <- PIK3CB[order(PIK3CB$nmf_subtypes),]
PIK3CB <- PIK3CB[!duplicated(PIK3CB$Sample_ID),]

NOTCH1 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="NOTCH1")
NOTCH1 <- NOTCH1[order(NOTCH1$nmf_subtypes),]
NOTCH1 <- NOTCH1[!duplicated(NOTCH1$Sample_ID),]

NF1 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="NF1")
NF1 <- NF1[order(NF1$nmf_subtypes),]
NF1 <- NF1[!duplicated(NF1$Sample_ID),]

MYCN <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="MYCN")
MYCN <- MYCN[order(MYCN$nmf_subtypes),]
MYCN <- MYCN[!duplicated(MYCN$Sample_ID),]

EGFR <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="EGFR")
EGFR <- EGFR[order(EGFR$nmf_subtypes),]
EGFR <- EGFR[!duplicated(EGFR$Sample_ID),]

FGF6 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="FGF6")
FGF6 <- FGF6[order(FGF6$nmf_subtypes),]
FGF6 <- FGF6[!duplicated(FGF6$Sample_ID),]

FGF3 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="FGF3")
FGF3 <- FGF3[order(FGF3$nmf_subtypes),]
FGF3 <- FGF3[!duplicated(FGF3$Sample_ID),]

FGF6 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="FGF6")
FGF6 <- FGF6[order(FGF6$nmf_subtypes),]
FGF6 <- FGF6[!duplicated(FGF6$Sample_ID),]

FGF10 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="FGF10")
FGF10 <- FGF10[order(FGF10$nmf_subtypes),]
FGF10 <- FGF10[!duplicated(FGF10$Sample_ID),]

SOX2 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="SOX2")
SOX2 <- SOX2[order(SOX2$nmf_subtypes),]
SOX2 <- SOX2[!duplicated(SOX2$Sample_ID),]

MYC <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="MYC")
MYC <- MYC[order(MYC$nmf_subtypes),]
MYC <- MYC[!duplicated(MYC$Sample_ID),]

PacientesLusc <- GrupoK3_Lusc$Sample_ID
MatrizPerfilMutagenico <- matrix(0,nrow=length(HugoNames)+1,ncol=length(PacientesLusc))
MatrizPerfilMutagenico <- as.data.frame(MatrizPerfilMutagenico)
colnames(MatrizPerfilMutagenico) <- PacientesLusc
rownames(MatrizPerfilMutagenico) <- c("PertenenciaCluster",HugoNames)
MatrizPerfilMutagenico[1,] <- GrupoK3_Lusc$nmf_subtypes

Mutaciones <- rbind(PIK3CA,PIK3CB)
Mutaciones <- rbind(Mutaciones,NOTCH1)
Mutaciones <- rbind(Mutaciones,NF1)
Mutaciones <- rbind(Mutaciones,MYCN)
Mutaciones <- rbind(Mutaciones,EGFR)
Mutaciones <- rbind(Mutaciones,FGF6)
Mutaciones <- rbind(Mutaciones,FGF3)
Mutaciones <- rbind(Mutaciones,FGF10)
Mutaciones <- rbind(Mutaciones,SOX2)
Mutaciones <- rbind(Mutaciones,MYC)

i <- 1
j <- 1
for (i in seq_along(HugoNames)){
  Mutacion <- HugoNames[i]
  DataFrameMutaciones <- Mutaciones[which(Mutaciones$Hugo_Symbol==Mutacion),]
  
  for (j in seq_along(PacientesLusc)){
    Paciente <- PacientesLusc[j]
    PacientePresenteEnMutacion <- DataFrameMutaciones[which(DataFrameMutaciones$Sample_ID==Paciente),]
    
    Comprobacion <- isTRUE(PacientePresenteEnMutacion$Sample_ID==Paciente)
    
    if (Comprobacion==TRUE){
      Comprobacion <- 1
    } else{
      Comprobacion <- 0
    }
    
    if (Comprobacion==1){MatrizPerfilMutagenico[i+1,j] <- 1} 
    else {next}
  }
}


ComprobacionMatriz_TP53 <- sum(as.numeric(MatrizPerfilMutagenico[2,]))

library(pheatmap)


#Hacemos la traspuesta de la matriz para tener los genes en las columnas y a los pacientes en las filas. 
Matriz <- data.frame(t(MatrizPerfilMutagenico))
#creamos una nueva columna con el nombre de los pacientes en las filas


#Eliminamos la pertenencia para poder hacer el pivot longer
Matriz_SinPertenencia <- Matriz[,-1]

#Cambiamos la organizacion de la matriz
Matriz_Longer <- Matriz_SinPertenencia %>% mutate(Row = rownames(Matriz_SinPertenencia)) %>%
  pivot_longer(-Row) 
#Cambiamos el nombre de la columna del nombre del paciente para luego poder hacer el inner Join y recuperar la pertenencia
colnames(Matriz_Longer)[1] <- "SampleID"

Matriz <- Matriz %>% mutate(SampleID=rownames(Matriz))
Matriz_Longer_Pertenencia <- inner_join(Matriz_Longer,Matriz,by="SampleID")
Matriz_Longer_Pertenencia <- Matriz_Longer_Pertenencia[,-c(5:15)]



MatrizDatos <- Matriz[,c(2:12)]
MatrizDatos[,1] <- as.numeric(MatrizDatos[,1])
MatrizDatos[,2] <- as.numeric(MatrizDatos[,2])
MatrizDatos[,3] <- as.numeric(MatrizDatos[,3])
MatrizDatos[,4] <- as.numeric(MatrizDatos[,4])
MatrizDatos[,5] <- as.numeric(MatrizDatos[,5])
MatrizDatos[,6] <- as.numeric(MatrizDatos[,6])
MatrizDatos[,7] <- as.numeric(MatrizDatos[,7])
MatrizDatos[,8] <- as.numeric(MatrizDatos[,8])
MatrizDatos[,9] <- as.numeric(MatrizDatos[,9])
MatrizDatos[,10] <- as.numeric(MatrizDatos[,10])
MatrizDatos[,11] <- as.numeric(MatrizDatos[,11])

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
                annotation_col =PertenenciaDF,
                annotation_legend = T,
                main="Genes grouped by mutations")

NombresClusteres <- c("Cluster 1","Cluster 2","Cluster 3")
names(NombresClusteres) <- c("1","2","3")
PerfilMutagenico <- ggplot(data = Matriz_Longer_Pertenencia, 
                           mapping = aes(x = SampleID,
                                         y = reorder(name, desc(name)),
                                         fill = as.factor(value))) +
  geom_tile() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(name = "Presencia Mutacion", values = c("grey","black"))+
  xlab(label = "Sample") +
  ylab(label="Supresor Genes")+
  facet_grid(~ PertenenciaCluster, 
             scales = "free_x", 
             space = "free_x",
             labeller = labeller(PertenenciaCluster = NombresClusteres))+
  ggtitle(label = "Perfil Mutagenico Lusc")

PerfilMutagenico 


