directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/Objetos")

silhouette <- readRDS("res_nb_silhouette_10.rds")
dindex <- readRDS("res_nb_dindex_luad_10.rds")
dunn <- readRDS("res_nb_dunn_luad_10.rds")
gap <- readRDS("res_nb_gap_luad_10.rds")
tau <- readRDS("res_nb_tau_luad_10.rds")

indicesSilhouette <- matrix(0,10,9)
indicesDindex <- matrix(0,10,9)
indicesDunn <- matrix(0,10,9)
indicesGap <- matrix(0,10,9)
indicesTau <- matrix(0,10,9)

maximoSilhouette <- matrix(0,nrow=10,ncol=3)
maximoTau <- matrix(0,nrow=10,ncol=3)
maximoGap <- matrix(0,nrow=10,ncol=3)
maximoDunn <- matrix(0,nrow=10,ncol=3)

i <- 1
for (i in seq_along(silhouette)){
  maximo <- max(silhouette[[i]])
  maximoSilhouette[i,1] <- maximo
  indices <- which(silhouette[[i]]==maximo)+1
  maximoSilhouette[i,2] <- indices[1]
  maximoSilhouette[i,3] <- indices[2]
  
  maximo <- max(tau[[i]])
  maximoTau[i,1] <- maximo
  indices <- which(tau[[i]]==maximo)+1
  maximoTau[i,2] <- indices[1]
  maximoTau[i,3] <- indices[2]
  
  maximo <- max(gap[[i]])
  maximoGap[i,1] <- maximo
  indices <- which(gap[[i]]==maximo)+1
  maximoGap[i,2] <- indices[1]
  maximoGap[i,3] <- indices[2]
  
  maximo <- max(dunn[[i]])
  maximoDunn[i,1] <- maximo
  indices <- which(dunn[[i]]==maximo)+1
  maximoDunn[i,2] <- indices[1]
  maximoDunn[i,3] <- indices[2]
  
}

maximoDunn <- data.frame(maximoDunn)[,-3]
names(maximoDunn) <- c("ValorMaximoDunn","NumeroClusterDunn")

maximoSilhouette <- data.frame(maximoSilhouette)
names(maximoSilhouette) <- c("ValorMaximoSilhouette","NumeroCluster","NumeroClusterSilhouette")

maximoGap <- data.frame(maximoGap)[,-3]
names(maximoGap) <- c("ValorMaximoGap","NumeroClusterGap")

maximoTau <- data.frame(maximoTau)[,-3]
names(maximoTau) <- c("ValorMaximoTau","NumeroClusterTau")

indicesTotales <- cbind(maximoDunn,maximoGap,maximoSilhouette,maximoTau)



numeroClusteres <- factor(rep(2:10,each=10))

i <- 1
for (i in seq_along(silhouette)){
  indicesSilhouette[i,]<- (silhouette[[i]])
  indicesDindex[i,]<- (dindex[[i]])
  indicesDunn[i,]<- (dunn[[i]])
  indicesGap[i,]<- (gap[[i]])
  indicesTau[i,]<- (tau[[i]])
}

indicesSilhouettefilaUnica <- as.vector(indicesSilhouette)
indicesDindexfilaUnica <- as.vector(indicesDindex)
indicesDunnfilaUnica <- as.vector(indicesDunn)
indicesGapfilaUnica <- as.vector(indicesGap)
indicesTaufilaUnica <- as.vector(indicesTau)

SILHOUETTE <- data.frame(numeroClusteres=numeroClusteres,valores=indicesSilhouettefilaUnica)
DINDEX<- data.frame(numeroClusteres=numeroClusteres,valores=indicesDindexfilaUnica)
DUNN<- data.frame(numeroClusteres=numeroClusteres,valores=indicesDunnfilaUnica)
GAP<- data.frame(numeroClusteres=numeroClusteres,valores=indicesGapfilaUnica)
TAU <- data.frame(numeroClusteres=numeroClusteres,valores=indicesTaufilaUnica)

library(tidyverse)
GraficoSilhouette <- ggplot(data=SILHOUETTE,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_boxplot()+
  ylab("Indice Silhouette")+
  xlab("Numero de Clusteres")

GraficoDindex <- ggplot(data=DINDEX,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_boxplot()+
  ylab("Indice Dindex")+
  xlab("Numero de Clusteres")

GraficoDunn <- ggplot(data=DUNN,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_boxplot()+
  ylab("Indice Dunn")+
  xlab("Numero de Clusteres")

GraficoGap <- ggplot(data=GAP,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_boxplot()+
  ylab("Indice Gap")+
  xlab("Numero de Clusteres")

GraficoTau <- ggplot(data=TAU,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_boxplot()+
  ylab("Indice Tau")+
  xlab("Numero de Clusteres")

library(ggpubr)
figura <- ggarrange(GraficoSilhouette,GraficoDindex,GraficoDunn,GraficoGap,GraficoTau,
                    labels = c("A", "B", "C","D","E"), ncol = 3, nrow = 2)

annotate_figure(figura,
                top = text_grob("Validacion Clusteres Kmeans LUAD", color = "Black", size = 14))

