# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-06-2021

#OBJETIVO:Es una prueba para intentar representar las lineas celulares de otra forma 
#en el que se ven mas claros si clusterizan o no con esos genes. 

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
#######################################################################################################################################

# FIRST APPROACH: 

#En este primer approach a partir de la matriz lo que voy ha hacer es calcular la media de expresion de los genes en cada una de las lineas
#y con respecto a esa media organizar las lineas celulares de mi matriz. 

Matriz <- readRDS("Matriz_Huellas_Lineas.rds")
Media <- colMeans(Matriz)
names(Media) <- colnames(Matriz)

Media <- sort(Media)

NuevaMatriz <- matrix(0,nrow=nrow(Matriz),ncol=ncol(Matriz))
colnames(NuevaMatriz) <- names(Media)
rownames(NuevaMatriz) <- rownames(Matriz)
i <- 1

for ( i in seq_along(Media)){
  
  LineaCelular <- names(Media)[i]
  
  CuentasLineaCelular <- Matriz[,which(colnames(Matriz)==LineaCelular)]
  
  NuevaMatriz[,i] <- CuentasLineaCelular
}

saveRDS(NuevaMatriz,"Matriz_Lineas_Ordenadas_Segun_Media.rds")

Matriz <- NuevaMatriz

Anotaciones <- readRDS("Anotaciones_Heatmap.rds")

rownames(Matriz) <- Anotaciones[,1]

H <- Matriz
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H, 
         color = Pallete1, 
         cluster_rows =F , 
         breaks = breaks,
         cluster_cols = F, 
         show_colnames = T,
         show_rownames = F,
         annotation_row=Anotaciones[2])


#Haciendo la traspuesta: 

H <- t(Matriz)

pheatmap(H, 
         color = Pallete1, 
         cluster_rows =T , 
         breaks = breaks,
         cluster_cols = F, 
         show_colnames = F,
         show_rownames = T,
         annotation_col = Anotaciones[2])

################################################################################################################################################
#En este caso vamos ha hacer una representacion en la que ponemos solo las lineas que me han clusterizado en el hclust. 

Matriz <- readRDS("Matriz_Huellas_Lineas.rds")
Anotaciones <- readRDS("Anotaciones_Heatmap.rds")

#Clusterizacion de las lienas: 
plot(hclust(dist(t(Matriz)))) #Clusteriza las filas de la matriz por lo tanto las lineas tienen que ir filas. 

#Clusterizacion de las huellas: 
plot(hclust(dist(Matriz)))

GrupoLineas <- c("SW900_LUNG","HCC15_LUNG",
                 "SKMES1_LUNG","CALU1_LUNG",
                 "LOUNH91_LUNG","NCIH1703_LUNG",
                 "SQ1_LUNG","EBC1_LUNG",
                 "KNS62_LUNG","EPLC272H_LUNG","NCIH2170_LUNG",
                 "NCIH226_LUNG","SW1573_LUNG")
NuevaMatriz <- matrix(0,ncol=length(GrupoLineas),nrow=nrow(Matriz))
colnames(NuevaMatriz) <- GrupoLineas
rownames(NuevaMatriz) <- rownames(Matriz)
i <- 1
for ( i in seq_along(GrupoLineas)){
  Linea <- GrupoLineas[i]
  NuevaMatriz[,i] <- Matriz[,which(colnames(Matriz)==Linea)]
}

rownames(NuevaMatriz) <- Anotaciones[,1]
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
#Esto es lineas en columnas 
pheatmap(NuevaMatriz, 
         color = Pallete1, 
         cluster_rows =F , 
         breaks = breaks,
         cluster_cols = T, 
         show_colnames = T,
         show_rownames = F,
         annotation_row = Anotaciones[2])

#Esto es lineas en filas: 

NuevaMatriz <- t(NuevaMatriz)
pheatmap(NuevaMatriz, 
         color = Pallete1, 
         cluster_rows =T , 
         breaks = breaks,
         cluster_cols = F, 
         show_colnames = F,
         show_rownames = T,
         annotation_col = Anotaciones[2])

#Tengo que encontrar la forma de comparar la expresion de las dos situaciones. Entre la huella y la expresion de dicha huella en las lineas





