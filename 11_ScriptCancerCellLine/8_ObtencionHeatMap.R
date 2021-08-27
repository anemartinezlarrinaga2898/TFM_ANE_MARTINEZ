# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 26-05-2021

#OBJETIVO:Script para obtener los heatmaps

################################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
Huellas <- readRDS("Huellas_LineasCelulares_Lusc.rds")

HuellasSinNan <- list()
for ( i in seq_along(Huellas)){
  Huella <- Huellas[[i]]
  
  Y <- na.omit(Huella)
  
  HuellasSinNan[[i]] <- Y
}

saveRDS(HuellasSinNan,"Huellas_LineasCelulares_Lusc_SinNan.rds")
################################################################################################################################################

Huellas <- readRDS("Huellas_LineasCelulares_Lusc_SinNan.rds")

################################################################################################################################################

#Primera visualizacion en la que represento los resultados tal cual, voy a intentar reducir "el ruido" un poco mas, intentado representar aquellos genes que tengan
#niveles de cuentas muy altos o niveles de cuentas muy bajos. Para no tener que repreentar todos los genes ya que si no tengo mucho ruido y relamente no se ve nada

# H E A T M A P :  G R U P O  1

#HeatMap donde predomina el negro
H1 <- t(Huellas[[1]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H1, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

#HeatMap donde predomina el verde
H1 <- t(Huellas[[1]])
Colors <- c("green","black","red")
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H1, color = Pallete1, cluster_rows =T ,cluster_cols = F, show_colnames = FALSE)

# H E A T M A P :  G R U P O  2

#HeatMap donde predomina el negro
H2 <- t(Huellas[[2]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H2, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

#HeatMap donde predomina el verde
H2 <- t(Huellas[[2]])
Colors <- c("green","black","red")
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H2, color = Pallete1, cluster_rows =T ,cluster_cols = F, show_colnames = FALSE)

# H E A T M A P :  G R U P O  3

#HeatMap donde predomina el negro
H3 <- t(Huellas[[3]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H3, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

#HeatMap donde predomina el verde
H3 <- t(Huellas[[3]])
Colors <- c("green","black","red")
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H3, color = Pallete1, cluster_rows =T ,cluster_cols = F, show_colnames = FALSE)

################################################################################################################################################

Huellas_Ordenadas_As <- list()
Huellas_Orden_Genes <- list()
i <- 1
for ( i in seq_along(Huellas)){
  Huella <- Huellas[[i]]
  Huella <- t(Huella)
  
  j <- 1
  
  Filas<- nrow(Huella)
  Columnas <- ncol(Huella)
  
  HuellasOrden <- matrix(0,nrow = Filas, ncol = Columnas)
  HuellasOrden <- as.data.frame(HuellasOrden)
  rownames(HuellasOrden) <- rownames(Huella)
  
  OrdenGenes <- matrix(0,nrow = Filas, ncol = Columnas)
  OrdenGenes <-  as.data.frame(OrdenGenes)
  rownames(OrdenGenes) <- rownames(Huella)
  
  for ( j in seq(1,Filas,1)){
    LineaCelular <- Huella[j,]
    names(LineaCelular) <- colnames(Huella)
    
    Ordenado <- LineaCelular[order(LineaCelular)]
    
    HuellasOrden[j,] <-Ordenado
    OrdenGenes[j,] <- names(Ordenado)
  }
  
  Huellas_Ordenadas_As[[i]] <- HuellasOrden
  Huellas_Orden_Genes[[i]] <- OrdenGenes
}
  
saveRDS(Huellas_Ordenadas_As,"Huellas_Ordenadas_Lusc.rds")
saveRDS(Huellas_Orden_Genes,"OrdenGenes_Lusc.rds")

################################################################################################################################################

Huellas <- readRDS("Huellas_Ordenadas_Lusc.rds")

H1 <- (Huellas[[1]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H1, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

H2 <- (Huellas[[2]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H2, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

H2 <- (Huellas[[2]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H2, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)

H3 <- (Huellas[[3]])
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H3, color = Pallete1, cluster_rows =T , breaks = breaks,cluster_cols = F, show_colnames = FALSE)
