# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 22-04-2021

#OBJETIVO: Localizacion de los genes que se encuentran mutados dentro de nuestros pacientes en el caso de luad. Los genes 
#con los que vamos a trabajar son: p53,lkb1,atm,keap1,cdkn2a

########################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/")
library(tidyverse)
MutacionesLusc <- read.delim("~/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/7_ScriptsPerfilMutagenico/Datos/Lusc/gdc_download_20210421_082747.056397/95258183-63ea-4c97-ae29-1bae9ed06334/TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf", comment.char="#")

#Primero limpiamos el nombre del sample ID para poder hacer indexacion despues
MutacionesLusc$Tumor_Sample_Barcode <- gsub("(-01A).*|(-01B).*","",MutacionesLusc[,16]) #En este data.frame tengo a todos mis pacientes
saveRDS(MutacionesLusc,"TotalMutacionesLusc.rds")

GrupoK3_Lusc<- readRDS("K3_NMF_Lusc.rds")
GrupoK3_Lusc <- GrupoK3_Lusc[order(GrupoK3_Lusc$nmf_subtypes),]
GrupoK3_Lusc <- GrupoK3_Lusc[,-3]

PacientesLusc <- GrupoK3_Lusc$Sample_ID #Nombres de los pacientes de LUAD

# ······················· B U S Q U E D A  P A C I E N T E S ·····················
# B U C L E: ··································································
# Objetivo: en este bucle lo que hacemos es ir indexando en la tabla total de mutaciones con nuestros 61 pacientes. A la vez que vamos indexando vamos
# juntadolo en un unico data.frame para obtener el resultado final. 
i <- 1
MutacionesPacienteLusc <- data.frame()
for (i in seq_along(PacientesLusc)){
  Paciente <- PacientesLusc[i]
  MutacionesPaciente_i <- MutacionesLusc[which(MutacionesLusc[,16]==Paciente),]
  MutacionesPacienteLusc_Interes <- rbind(MutacionesPaciente_i,MutacionesPacienteLusc)
  MutacionesPacienteLusc <- MutacionesPacienteLusc_Interes
}

saveRDS(MutacionesPacienteLusc_Interes,"MutacionesPacientesLusc.rds")

# ······················· B U S Q U E D A  M U T A C I O N E S ·····················
# B U C L E: ··································································
#Objetivo: al igual que con los pacientes lo que hacemos es ir indexando por los diferentes identificadores de nuestros genes supresores de tumores
#de interes. En este caso hemos filtrado por el entrezID, pero podria haber sido por el HUGO SYMBOL o por el ENSEMBEL

HugoNames <- c("TP53","STK11","PTEN","KEAP1","CDKN2A")
EntrezID <- c(7157,6794,5728,9817,1029)

GenesMutados <- data.frame()
i <- 1

for (i in seq_along(EntrezID)){
  Mutacion <- EntrezID[i]
  Mutacion_i <- MutacionesPacienteLusc_Interes[which(MutacionesPacienteLusc_Interes$Entrez_Gene_Id==Mutacion),]
  GenesMutadosInteres <- rbind(Mutacion_i,GenesMutados)
  GenesMutados <- GenesMutadosInteres
}

saveRDS(GenesMutadosInteres,"GenesMutadosPacientes_Lusc.rds")

# Por lo tanto tenemos el DATA.FRAME final donde tengo a mis pacientes con las mutaciones que presentan.

#Voy a unir la matriz de los genes con la de la pertenencia a los clusteres, para tener a que grupo pertenece cada uno. 

colnames(GenesMutadosInteres)[16] <- "Sample_ID"
GenesMutadosInteres_PertenenciaCluster <- inner_join(GrupoK3_Lusc,GenesMutadosInteres,by="Sample_ID")

# Obtencion de un data.frame por cada uno de los genes alterados los cuales organizamos para tener la pertenencia a cada cluster. 
#1) Filtramos por gen
#2) Ordenamos por pertenencia a cluster 
#3) Eliminamos duplicados. 


TP53 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="TP53")
TP53 <- TP53[order(TP53$nmf_subtypes),]
TP53 <- TP53[!duplicated(TP53$Sample_ID),]

STK11 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="STK11")
STK11 <- STK11[order(STK11$nmf_subtypes),]
STK11 <- STK11[!duplicated(STK11$Sample_ID),]

PTEN <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="PTEN")
PTEN <- PTEN[order(PTEN$nmf_subtypes),]
PTEN <- PTEN[!duplicated(PTEN$Sample_ID),]

KEAP1 <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="KEAP1")
KEAP1 <- KEAP1[order(KEAP1$nmf_subtypes),]
KEAP1 <- KEAP1[!duplicated(KEAP1$Sample_ID),]

CDKN2A <- GenesMutadosInteres_PertenenciaCluster %>% filter(Hugo_Symbol=="CDKN2A")
CDKN2A <- CDKN2A[order(CDKN2A$nmf_subtypes),]
CDKN2A <- CDKN2A[!duplicated(CDKN2A$Sample_ID),]


# ······················ M A T R I Z  P E R F I L   M U T A G E N I C O (obtencion) ·················
MatrizPerfilMutagenico <- matrix(0,nrow=length(HugoNames)+1,ncol=length(PacientesLusc))
MatrizPerfilMutagenico <- as.data.frame(MatrizPerfilMutagenico)
colnames(MatrizPerfilMutagenico) <- PacientesLusc
rownames(MatrizPerfilMutagenico) <- c("PertenenciaCluster",HugoNames)
MatrizPerfilMutagenico[1,] <- GrupoK3_Lusc$nmf_subtypes


#Juntar todas las mutaciones en un unico data.frame 

Mutaciones <- rbind(TP53,STK11)
Mutaciones <- rbind(Mutaciones,PTEN)
Mutaciones <- rbind(Mutaciones,KEAP1)
Mutaciones <- rbind(Mutaciones,CDKN2A)

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

#Mutaciones de cada grupo: 
#1) TP53: 31
#2) STK11:8
#3) ATM: 3
#4) KEAP1: 8
#5) CDKN2A: 2


ComprobacionMatriz_TP53 <- sum(as.numeric(MatrizPerfilMutagenico[2,])) #31
ComprobacionMatriz_STK11 <- sum(as.numeric(MatrizPerfilMutagenico[3,])) #8
ComprobacionMatriz_PTEN <- sum(as.numeric(MatrizPerfilMutagenico[4,])) #3
ComprobacionMatriz_KEAP1 <- sum(as.numeric(MatrizPerfilMutagenico[5,])) #8
ComprobacionMatriz_CDKN2A <- sum(as.numeric(MatrizPerfilMutagenico[6,])) #2

saveRDS(MatrizPerfilMutagenico,"MatrizPerfilMutagenico_Lusc.rds")


# ······················ M A T R I Z  P E R F I L   M U T A G E N I C O (Generacion) ·················

# 1) Generacion de la matriz con el paquete de ggplot: ···············································

MatrizPerfilMutagenico <- readRDS("MatrizPerfilMutagenico_Lusc.rds")

#Hacemos la traspuesta de la matriz para tener los genes en las columnas y a los pacientes en las filas. 
Matriz <- data.frame(t(MatrizPerfilMutagenico))

#Eliminamos la pertenencia para poder hacer el pivot longer
Matriz_SinPertenencia <- Matriz[,-1]

#Cambiamos la organizacion de la matriz
Matriz_Longer <- Matriz_SinPertenencia %>% mutate(Row = rownames(Matriz_SinPertenencia)) %>%
  pivot_longer(-Row) 
#Cambiamos el nombre de la columna del nombre del paciente para luego poder hacer el inner Join y recuperar la pertenencia
colnames(Matriz_Longer)[1] <- "SampleID"

Matriz <- Matriz %>% mutate(SampleID=rownames(Matriz))
Matriz_Longer_Pertenencia <- inner_join(Matriz_Longer,Matriz,by="SampleID")
Matriz_Longer_Pertenencia <- Matriz_Longer_Pertenencia[,-c(5,6,7,8,9)]

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

#2) Generacion de grafico con pheatmap: ········································

library(pheatmap)

MatrizPerfilMutagenico <- readRDS("MatrizPerfilMutagenico_Lusc.rds")

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
Matriz_Longer_Pertenencia <- Matriz_Longer_Pertenencia[,-c(5,6,7,8,9)]



MatrizDatos <- Matriz[,c(2,3,4,5,6)]
MatrizDatos[,1] <- as.numeric(MatrizDatos[,1])
MatrizDatos[,2] <- as.numeric(MatrizDatos[,2])
MatrizDatos[,3] <- as.numeric(MatrizDatos[,3])
MatrizDatos[,4] <- as.numeric(MatrizDatos[,4])
MatrizDatos[,5] <- as.numeric(MatrizDatos[,5])

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






