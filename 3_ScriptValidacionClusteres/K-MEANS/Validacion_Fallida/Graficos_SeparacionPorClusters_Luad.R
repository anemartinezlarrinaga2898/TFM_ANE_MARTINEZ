#Separar los datos por los datos obtenidos: 

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/Objetos")

library(tidyverse)
library(NbClust)
library(cluster)
library(Rtsne)
library(tsne)

muestras <- readRDS('MuestrasLuadTraspuesta.rds')
silhouette <- readRDS("res_nb_silhouette_luad.rds")
dindex <- readRDS("res_nb_dindex_luad.rds")
dunn <- readRDS("res_nb_dunn_luad.rds")
gap <- readRDS("res_nb_gap_luad.rds")
tau <- readRDS("res_nb_tau_luad.rds")

particionSilhouete <- as.factor(silhouette$Best.partition)
particionDindex <- as.factor(dindex$Best.partition)
particionDunn <- as.factor(dunn$Best.partition)
particionGap <- as.factor(gap$Best.partition)
particionTau <- as.factor(tau$Best.partition)

tsne_coord <- Rtsne::Rtsne(as.matrix(muestras),perplexity = 10,check_duplicates=FALSE)
tsne_coord <- as.data.frame(tsne_coord$Y)

pacientes <- seq(1,70,1)

datosLuadTotales <- data.frame(pacientes=pacientes,
                               coordenadas=tsne_coord,
                               Silhouette=particionSilhouete,
                               Dunn=particionDunn,
                               Gap=particionGap,
                               Tau=particionTau)
                               
names(datosLuadTotales) <- c("Paciente","CoordenadaX","CoordenadaY", "Silhouette","Dunn","Gap","Tau")

graficoSilhouette <- ggplot(datosLuadTotales, aes(x=CoordenadaX,  y=CoordenadaY,col=Silhouette)) + geom_point()
graficoDunn <- ggplot(datosLuadTotales, aes(x=CoordenadaX,  y=CoordenadaY,col=Dunn)) + geom_point()
graficoGap<- ggplot(datosLuadTotales, aes(x=CoordenadaX,  y=CoordenadaY,col=Gap)) + geom_point()
graficoTau <- ggplot(datosLuadTotales, aes(x=CoordenadaX,  y=CoordenadaY,col=Tau)) + geom_point()

figura <- ggarrange(graficoDunn, graficoGap, graficoSilhouette,graficoTau,
                    labels = c("A", "B", "C","D"), ncol = 2, nrow = 2)

annotate_figure(figura,
                top = text_grob("Validacion Clusteres Kmeans LUAD", color = "Black", size = 14))
