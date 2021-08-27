# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO:Generacion de graficos KMEANS

# A N A L I S I S  | R E S U L T A D O S: 

#TM: Todas las muestras
directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/4_ScriptsAnalisisResultados_ValidacionClusteres/")
library(ggpubr)
library(diceR)
library(tidyverse)

ResultadosLusc_TM <- readRDS("KMEANS/ResutladosKmean_TodasLasMuestras_Lusc.rds")
ResultadosLuad_TM <- readRDS("KMEANS/ResutladosKmean_TodasLasMuestras_Luad.rds")

muestrasLuad <- readRDS("Luad_Escalado.rds")
muestrasLuad <- t(muestrasLuad)
muestrasLusc <- readRDS("Lusc_Escalado.rds")
muestrasLusc <- t(muestrasLusc)


##############################################################################
# L U A D
K2_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[1]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K2_Luad_TM

K3_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[2]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K3_Luad_TM

K4_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[3]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K4_Luad_TM

K5_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[4]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K5_Luad_TM

K6_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[5]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K6_Luad_TM

K7_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[6]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K7_Luad_TM

K8_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[7]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K8_Luad_TM

K9_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[8]], 
                           data = muestrasLuad,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K9_Luad_TM

K10_Luad_TM <- fviz_cluster(ResultadosLuad_TM[[9]], 
                            data = muestrasLuad,
                            geom = "point",
                            ellipse.type = "convex", 
                            ggtheme = theme_bw()
)
K10_Luad_TM

figura <- ggarrange(K2_Luad_TM,K3_Luad_TM,K4_Luad_TM,K5_Luad_TM,K6_Luad_TM,K7_Luad_TM,K8_Luad_TM,K9_Luad_TM,K10_Luad_TM,
                    labels = c("K=2", "K=3", "K=4","K=5","K=6","K=7","K=8","K=9","K=10"), ncol = 5, nrow = 2)
annotate_figure(figura,
                top = text_grob("Luad, K-Means Clustering", color = "Black", size = 14))

  # Grafico de los grupos de Kmeans Luad
##############################################################################

#LUSC

K2_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[1]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)
K2_Lusc_TM

K3_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[2]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K3_Lusc_TM

K4_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[3]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K4_Lusc_TM

K5_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[4]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K5_Lusc_TM

K6_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[5]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K6_Lusc_TM

K7_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[6]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K7_Lusc_TM

K8_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[7]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K8_Lusc_TM

K9_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[8]], 
                           data = muestrasLusc,
                           geom = "point",
                           ellipse.type = "convex", 
                           ggtheme = theme_bw()
)

K9_Lusc_TM

K10_Lusc_TM <- fviz_cluster(ResultadosLusc_TM[[9]], 
                            data = muestrasLusc,
                            geom = "point",
                            ellipse.type = "convex", 
                            ggtheme = theme_bw()
)

K10_Lusc_TM

figura <- ggarrange(K2_Lusc_TM,K3_Lusc_TM,K4_Lusc_TM,K5_Lusc_TM,K6_Lusc_TM,K7_Lusc_TM,K8_Lusc_TM,K9_Lusc_TM,K10_Lusc_TM,
                    labels = c("K=2", "K=3", "K=4","K=5","K=6","K=7","K=8","K=9","K=10"), ncol = 5, nrow = 2)
annotate_figure(figura,
                top = text_grob("Lusc, K-Means Clustering", color = "Black", size = 14))
 # Grafico de los grupos de Kmeans Lusc
##############################################################################
directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/4_ScriptsAnalisisResultados_ValidacionClusteres/KMEANS/")

NumeroCluster <- seq(2,7,1)


I.Silhouette_TM_Luad <- readRDS("KMEANS/IndiceSilhouette_TodasLasMuestras_Luad.rds")
I.Silhouette_TM_Luad <- round(I.Silhouette_TM_Luad,digits=3)
I.Silhouette_TM_Lusc <- readRDS("KMEANS/IndiceSilhouette_TodasLasMuestras_Lusc.rds")
I.Silhouette_TM_Lusc <- round(I.Silhouette_TM_Lusc,digits=3)

I.Varianza_TM_Luad <- readRDS("KMEANS/IndiceVarianza_TodasLasMuestras_Luad.rds")
I.Varianza_TM_Luad <- round(I.Varianza_TM_Luad,digits=3)
I.Varianza_TM_Lusc <- readRDS("KMEANS/IndiceVarianza_TodasLasMuestras_Lusc.rds")
I.Varianza_TM_Lusc <- round(I.Varianza_TM_Lusc,digits=3)

I.Dissimilarity_TM_Luad <- readRDS("KMEANS/IndiceDissimilarity_TodasLasMuestras_Luad.rds")
I.Dissimilarity_TM_Luad <- round(I.Dissimilarity_TM_Luad,digits=3)
I.Dissimilarity_TM_Lusc <- readRDS("KMEANS/IndiceDissimilarity_TodasLasMuestras_Lusc.rds")
I.Dissimilarity_TM_Lusc <- round(I.Dissimilarity_TM_Lusc,digits=3)

IndicesLuad <- cbind(NumeroCluster,I.Silhouette_TM_Luad)
IndicesLuad <- cbind(IndicesLuad,I.Varianza_TM_Luad)
IndicesLuad <- cbind(IndicesLuad,I.Dissimilarity_TM_Luad)

IndicesLusc <- cbind(NumeroCluster,I.Silhouette_TM_Lusc)
IndicesLusc <- cbind(IndicesLusc,I.Varianza_TM_Lusc)
IndicesLusc <- cbind(IndicesLusc,I.Dissimilarity_TM_Lusc)

IndicesLuad <- data.frame(IndicesLuad)
colnames(IndicesLuad) <- c("NumeroCluster","Silhouette","Varianza","Dissimilarity")

IndicesLusc <- data.frame(IndicesLusc)
colnames(IndicesLusc) <- c("NumeroCluster","Silhouette","Varianza","Dissimilarity")

IndicesLuad_Graficos_S <- ggplot(data=IndicesLuad)+
  geom_line(aes(x=NumeroCluster,y=Silhouette),color="#26798E",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Silhouette),color="#26798E",size=2)+
  geom_text(aes(x=NumeroCluster,y=Silhouette,label=Silhouette),hjust=0.6, vjust=-0.5)

IndicesLuad_Graficos_V <- ggplot(data=IndicesLuad)+
  geom_line(aes(x=NumeroCluster,y=Varianza),color="#63CAA7",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Varianza),color="#63CAA7",size=2)+
  geom_text(aes(x=NumeroCluster,y=Varianza,label=Varianza),hjust=0.6, vjust=-0.5)

IndicesLuad_Graficos_D <- ggplot(data=IndicesLuad)+
  geom_line(aes(x=NumeroCluster,y=Dissimilarity),color="#FFC172",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Dissimilarity),color="#FFC172",size=2)+
  geom_text(aes(x=NumeroCluster,y=Dissimilarity,label=Dissimilarity),hjust=0.6, vjust=-0.5)

figuraLuad <- ggarrange(IndicesLuad_Graficos_S,IndicesLuad_Graficos_V,IndicesLuad_Graficos_D,
                        labels = c("A","B","C"), ncol = 3, nrow = 1)
annotate_figure(figuraLuad,
                top = text_grob("Loss Function of K-Means, Luad", color = "Black", size = 14))


IndicesLusc_Graficos_S <- ggplot(data=IndicesLusc)+
  geom_line(aes(x=NumeroCluster,y=Silhouette),color="#75BDE0",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Silhouette),color="#75BDE0",size=2)+
  geom_text(aes(x=NumeroCluster,y=Silhouette,label=Silhouette),hjust=0.6, vjust=-0.5)

IndicesLusc_Graficos_V <- ggplot(data=IndicesLusc)+
  geom_line(aes(x=NumeroCluster,y=Varianza),color="#F8D49B",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Varianza),color="#F8D49B",size=2)+
  geom_text(aes(x=NumeroCluster,y=Varianza,label=Varianza),hjust=0.6, vjust=-0.5)

IndicesLusc_Graficos_D <- ggplot(data=IndicesLusc)+
  geom_line(aes(x=NumeroCluster,y=Dissimilarity),color="#F89B9B",size=1.5,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Dissimilarity),color="#F89B9B",size=2)+
  geom_text(aes(x=NumeroCluster,y=Dissimilarity,label=Dissimilarity),hjust=0.6, vjust=-0.5)

figuraLusc <- ggarrange(IndicesLusc_Graficos_S,IndicesLusc_Graficos_V,IndicesLusc_Graficos_D,
                        labels = c("A","B","C"), ncol = 3, nrow = 1)
annotate_figure(figuraLusc,
                top = text_grob("Loss Function of K-Means, Lusc", color = "Black", size = 14))
 #Grafico de los tres indices, unica vuelta
############################################################################################################################################################    

#Analisis de los indices con partes de las muestras: 
directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/4_ScriptsAnalisisResultados_ValidacionClusteres/KMEANS")

#LUAD ···········································································

silhouette_Luad <- readRDS("IndiceSilhouette_ParteDeLasMuestras_Luad.rds")
dissimilarity_Luad <- readRDS("IndiceSilhouette_ParteDeLasMuestras_Luad.rds")
varianza_Luad <- readRDS("IndiceVarianza_ParteDeLasMuestras_Luad.rds")

indiceS_Luad <- matrix(0,100,6) #Indice Silhouette
indiceD_Luad <- matrix(0,100,6) #Indice Dissimilarity
indiceV_Luad <- matrix(0,100,6) #Indice Varianza

i <- 1
for (i in seq_along(silhouette_Luad)){
    indiceS_Luad[i,]<- (silhouette_Luad[[i]])
    indiceD_Luad[i,]<- (dissimilarity_Luad[[i]])
    indiceV_Luad[i,]<- (varianza_Luad[[i]])
  
}

numeroClusteres <- factor(rep(2:7,each=100))

indicesS_filaUnica_Luad <- as.vector(indiceS_Luad)
indicesD_filaUnica_Luad <- as.vector(indiceD_Luad)
indicesV_filaUnica_Luad <- as.vector(indiceV_Luad)

SILHOUETTE_LUAD <- data.frame(numeroClusteres=numeroClusteres,
                         valores=indicesS_filaUnica_Luad)

DISSIMILARITY_LUAD <- data.frame(numeroClusteres=numeroClusteres,
                         valores=indicesD_filaUnica_Luad)

VARIANZA_LUAD <- data.frame(numeroClusteres=numeroClusteres,
                         valores=indicesV_filaUnica_Luad)

library(tidyverse)

GraficoSilhouette_Luad <- ggplot(data=SILHOUETTE_LUAD,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Silhouette")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Silhouette para Luad, subsamples")

GraficoSilhouette_Luad

GraficoDissimilarity_Luad <- ggplot(data=DISSIMILARITY_LUAD,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Dissimilarity")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Dissimilarity para Luad, subsamples")

GraficoDissimilarity_Luad

GraficoVarianza_Luad <- ggplot(data=VARIANZA_LUAD,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Varianza")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Varianza para Luad, subsamples")

GraficoVarianza_Luad

####################################################################################################
#LUSC: 

directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/4_ScriptsAnalisisResultados_ValidacionClusteres/KMEANS")

silhouette_Lusc <- readRDS("IndiceSilhouette_ParteDeLasMuestras_Lusc.rds")
dissimilarity_Lusc <- readRDS("IndiceDissimilarity_ParteDeLasMuestras_Lusc.rds")
varianza_Lusc <- readRDS("IndiceVarianza_ParteDeLasMuestras_Lusc.rds")

indiceS_Lusc <- matrix(0,100,6) #Indice Silhouette
indiceD_Lusc <- matrix(0,100,6) #Indice Dissimilarity
indiceV_Lusc <- matrix(0,100,6) #Indice Varianza

i <- 1
for (i in seq_along(silhouette_Lusc)){
  indiceS_Lusc[i,]<- (silhouette_Lusc[[i]])
  indiceD_Lusc[i,]<- (dissimilarity_Lusc[[i]])
  indiceV_Lusc[i,]<- (varianza_Lusc[[i]])
}

numeroClusteres <- factor(rep(2:7,each=100))

indicesS_filaUnica_Lusc <- as.vector(indiceS_Lusc)
indicesD_filaUnica_Lusc <- as.vector(indiceD_Lusc)
indicesV_filaUnica_Lusc <- as.vector(indiceV_Lusc)

SILHOUETTE_LUSC <- data.frame(numeroClusteres=numeroClusteres,
                              valores=indicesS_filaUnica_Lusc)

DISSIMILARITY_LUSC <- data.frame(numeroClusteres=numeroClusteres,
                                 valores=indicesD_filaUnica_Lusc)

VARIANZA_LUSC <- data.frame(numeroClusteres=numeroClusteres,
                            valores=indicesV_filaUnica_Lusc)

library(tidyverse)

GraficoSilhouette_Lusc <- ggplot(data=SILHOUETTE_LUSC,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Silhouette")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Silhouette para Lusc, subsamples")

GraficoSilhouette_Lusc

GraficoDissimilarity_Lusc <- ggplot(data=DISSIMILARITY_LUSC,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Dissimilarity")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Dissimilarity para Lusc, subsamples")

GraficoDissimilarity_Lusc

GraficoVarianza_Lusc <- ggplot(data=VARIANZA_LUSC,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  ylab("Indice Varianza")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Varianza para Lusc, subsamples")

GraficoVarianza_Lusc

##############################################################################

