# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO:  Generacion de graficos para el algoritmo de NMF

#ANALISIS RESULTADOS NMF

library(tidyverse)
directoryI <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/4_ScriptsAnalisisResultados_ValidacionClusteres/NMF")
###########################################################################################################

Estimacion_Paquete_NMF_Luad <- readRDS("Resultados_NMF_Total_Luad.rds")
NumeroCluster <- seq(2,7,1)


##############################################################################################################

#ANALISIS TOTAL:  ·········································································································

ValorCofenetic_Total_Luad <- readRDS("ValorCofenetic_Total_Luad.rds")
ValorCofenetic_Total_Luad <- round(ValorCofenetic_Total_Luad,digits=3)
ValorCofenetic_Total_Luad <- data.frame(NumeroCluster=NumeroCluster,Valores=ValorCofenetic_Total_Luad)

Grafico_Luad_Cofenetic_Total <- ggplot(data=ValorCofenetic_Total_Luad,aes(x=NumeroCluster,y=Valores))+
  geom_line(color="#D0E6A5",size=2,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Valores),color="#D0E6A5",size=1.5)+
  geom_text(aes(x=NumeroCluster,y=Valores,label=Valores))+
  scale_x_continuous(breaks =seq(2,10,1))+
  labs(title="Cophenetic Correlation Luad")

Grafico_Luad_Cofenetic_Total

ValorCofenetic_Total_Lusc <- readRDS("ValorCofenetic_Total_Lusc.rds")
ValorCofenetic_Total_Lusc <- round(ValorCofenetic_Total_Lusc,digits=3)
ValorCofenetic_Total_Lusc <- data.frame(NumeroCluster=NumeroCluster,Valores=ValorCofenetic_Total_Lusc)

Grafico_Lusc_Cofenetic_Total <- ggplot(data=ValorCofenetic_Total_Lusc,aes(x=NumeroCluster,y=Valores))+
  geom_line(color="#CCABD8",size=2)+
  geom_point(aes(x=NumeroCluster,y=Valores),color="#CCABD8",size=1.5)+
  geom_text(aes(x=NumeroCluster,y=Valores,label=Valores))+
  scale_x_continuous(breaks =seq(2,10,1))+
  labs(title="Cophenetic Correlation Lusc")

Grafico_Lusc_Cofenetic_Total

####################################################################################################

#PARTE DE LAS MUESTRAS: 

numeroClusteres <- factor(rep(2:7,each=100))
cofeneticValue_Luad <- readRDS("ValorCofenetic_Parte_De_Las_Muestr_NMF_Luad.rds")
cofeneticValue_Luad <- as.matrix(cofeneticValue_Luad)

cofenetic_filaUnica_Luad <- as.vector(cofeneticValue_Luad)
COFENETIC_LUAD <- data.frame(numeroClusteres=numeroClusteres,
                              valores=cofenetic_filaUnica_Luad)

library(tidyverse)

GraficoCofenetic_Luad <- ggplot(data=COFENETIC_LUAD,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  ylab("Indice Cofenetic")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Cofenetic para Luad, subsamples")

GraficoCofenetic_Luad


numeroClusteres <- factor(rep(2:7,each=100))
cofeneticValue_Lusc <- readRDS("ValorCofenetic_Parte_De_Las_Muestr_NMF_Lusc.rds")
cofeneticValue_Lusc<- as.matrix(cofeneticValue_Lusc)

cofenetic_filaUnica_Lusc <- as.vector(cofeneticValue_Lusc)
COFENETIC_LUSC <- data.frame(numeroClusteres=numeroClusteres,
                             valores=cofenetic_filaUnica_Lusc)



GraficoCofenetic_Lusc <- ggplot(data=COFENETIC_LUSC,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  ylab("Indice Cofenetic")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Cofenetic para Lusc, subsamples")

GraficoCofenetic_Lusc

##########################################################################################################

# L U A D: 

MatricesConsenso_Luad <- readRDS("MatricesConsenso_Total_Luad.rds")

library(pheatmap)
library(ggplotify)
library(patchwork)

MatrizConsenso_Luad_K2 <- as.ggplot(pheatmap(MatricesConsenso_Luad[[1]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=2, Consensus Matrix"))


MatrizConsenso_Luad_K3 <- as.ggplot(pheatmap(MatricesConsenso_Luad[[2]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=3, Consensus Matrix"))

MatrizConsenso_Luad_K4 <- as.ggplot(pheatmap(MatricesConsenso_Luad[[3]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=4, Consensus Matrix"))

MatrizConsenso_Luad_K5 <- as.ggplot(pheatmap(MatricesConsenso_Luad[[4]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=5, Consensus Matrix"))

MatrizConsenso_Luad_K6 <- as.ggplot(pheatmap(MatricesConsenso_Luad[[5]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=6, Consensus Matrix"))

MatrizConsenso_Luad_K7 <-as.ggplot(pheatmap(MatricesConsenso_Luad[[6]],
                                   clustering_distance_rows="euclidean",
                                   clustering_distance_cols="euclidean",
                                   main="NMF Luad K=7, Consensus Matrix"))



figuraLuad_NMF_K2.K5 <- ggarrange(MatrizConsenso_Luad_K2,MatrizConsenso_Luad_K3,MatrizConsenso_Luad_K4,MatrizConsenso_Luad_K5,
                        labels = c("A","B","C","D"), ncol = 2, nrow = 2)
annotate_figure(figuraLuad_NMF_K2.K5,
                top = text_grob("NMF Consensus Clustering Luad ", color = "Black", size = 14))


figuraLuad_NMF_K6.K10 <- ggarrange(MatrizConsenso_Luad_K6,MatrizConsenso_Luad_K7,
                                  labels = c("A","B"), ncol = 2, nrow = 1)
annotate_figure(figuraLuad_NMF_K6.K10,
                top = text_grob("NMF Consensus Clustering Luad ", color = "Black", size = 14))


par(mfrow=c(2,2))
hist(MatricesConsenso_Luad[[1]], main="Comparacion K2,K3,K5",col="white")
hist(MatricesConsenso_Luad[[2]], add = TRUE, col = rgb(0, 1, 0, alpha = 0.3))
hist(MatricesConsenso_Luad[[4]], add = TRUE, col = rgb(1, 0, 0, alpha = 0.5))
hist(MatricesConsenso_Luad[[1]], main="K2-Luad-NMF, Frecuencias de la matriz consenso",col="white")
hist(MatricesConsenso_Luad[[2]], main="K3-Luad-NMF, Frecuencias de la matriz consenso",col=rgb(0, 1, 0, alpha = 0.3))
hist(MatricesConsenso_Luad[[4]], main="K5-Luad-NMF, Frecuencias de la matriz consenso",col=rgb(1, 0, 0, alpha = 0.3))


#LUSC: 


MatricesConsenso_Lusc <- readRDS("MatricesConsenso_Total_Lusc.rds")

MatrizConsenso_Lusc_K2 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[1]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=2, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))


MatrizConsenso_Lusc_K3 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[2]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=3, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))

MatrizConsenso_Lusc_K4 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[3]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=4, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))

MatrizConsenso_Lusc_K5 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[4]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=5, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))

MatrizConsenso_Lusc_K6 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[5]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=6, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))

MatrizConsenso_Lusc_K7 <- as.ggplot(pheatmap(MatricesConsenso_Lusc[[6]],
                                             clustering_distance_rows="euclidean",
                                             clustering_distance_cols="euclidean",
                                             main="NMF Lusc K=7, Consensus Matrix",
                                             show_colnames = F, show_rownames = F))


figuraLusc_NMF_K2.K5 <- ggarrange(MatrizConsenso_Lusc_K2,MatrizConsenso_Lusc_K3,MatrizConsenso_Lusc_K4,MatrizConsenso_Lusc_K5,
                                  labels = c("A","B","C","D"), ncol = 2, nrow = 2)
annotate_figure(figuraLusc_NMF_K2.K5,
                top = text_grob("NMF Consensus Clustering Lusc ", color = "Black", size = 14))


figuraLusc_NMF_K6.K10 <- ggarrange(MatrizConsenso_Lusc_K6,MatrizConsenso_Lusc_K7,
                                   labels = c("A","B"), ncol = 2, nrow = 1)
annotate_figure(figuraLusc_NMF_K6.K10,
                top = text_grob("NMF Consensus Clustering Lusc ", color = "Black", size = 14))


par(mfrow=c(2,2))
hist(MatricesConsenso_Luad[[1]], main="K2-Luad-NMF, Frecuencias de la matriz consenso",col="white")
hist(MatricesConsenso_Luad[[2]], main="K3-Luad-NMF, Frecuencias de la matriz consenso",col=rgb(0, 1, 0, alpha = 0.3))
hist(MatricesConsenso_Luad[[3]], main="K4-Luad-NMF, Frecuencias de la matriz consenso",col=rgb(0, 0, 1, alpha = 0.3))
hist(MatricesConsenso_Luad[[4]], main="K5-Luad-NMF, Frecuencias de la matriz consenso",col=rgb(1, 0, 0, alpha = 0.3))



#############################################################################

# M A T R I C E S  C O N S E N S O  L U A D: 

library(NMF)
resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

K2 <- resultadosNMFLuad[[1]]
K3 <- resultadosNMFLuad[[2]]
K4 <- resultadosNMFLuad[[3]]
K5 <- resultadosNMFLuad[[4]]
K6 <- resultadosNMFLuad[[5]]
K7 <- resultadosNMFLuad[[6]]

colnames(K2@consensus) <- sampleNames(K2)
rownames(K2@consensus) <- sampleNames(K2)
pdf("K2_LUAD.pdf")
consensusmap(K2, hclustfun="average")
dev.off()

colnames(K3@consensus) <- sampleNames(K3)
rownames(K3@consensus) <- sampleNames(K3)
pdf("K3_LUAD.pdf")
consensusmap(K3, hclustfun="average")
dev.off()

colnames(K4@consensus) <- sampleNames(K4)
rownames(K4@consensus) <- sampleNames(K4)
pdf("K4_LUAD.pdf")
consensusmap(K4, hclustfun="average")
dev.off()

colnames(K5@consensus) <- sampleNames(K5)
rownames(K5@consensus) <- sampleNames(K5)
pdf("K5_LUAD.pdf")
consensusmap(K5, hclustfun="average")
dev.off()

colnames(K6@consensus) <- sampleNames(K6)
rownames(K6@consensus) <- sampleNames(K6)
pdf("K6_LUAD.pdf")
consensusmap(K6, hclustfun="average")
dev.off()

colnames(K7@consensus) <- sampleNames(K7)
rownames(K7@consensus) <- sampleNames(K7)
pdf("K7_LUAD.pdf")
consensusmap(K7, hclustfun="average")
dev.off()

###########################################################################

# M A T R I C E S  C O N S E N S O  L U S C: 

library(NMF)
resultadosNMFLusc <- readRDS("Resultados_NMF_Total_Lusc.rds")

K2 <- resultadosNMFLusc[[1]]
K3 <- resultadosNMFLusc[[2]]
K4 <- resultadosNMFLusc[[3]]
K5 <- resultadosNMFLusc[[4]]
K6 <- resultadosNMFLusc[[5]]
K7 <- resultadosNMFLusc[[6]]

colnames(K2@consensus) <- sampleNames(K2)
rownames(K2@consensus) <- sampleNames(K2)
pdf("K2_LUSC.pdf")
consensusmap(K2, hclustfun="average")
dev.off()

colnames(K3@consensus) <- sampleNames(K3)
rownames(K3@consensus) <- sampleNames(K3)
pdf("K3_LUSC.pdf")
consensusmap(K3, hclustfun="average")
dev.off()

colnames(K4@consensus) <- sampleNames(K4)
rownames(K4@consensus) <- sampleNames(K4)
pdf("K4_LUSC.pdf")
consensusmap(K4, hclustfun="average")
dev.off()

colnames(K5@consensus) <- sampleNames(K5)
rownames(K5@consensus) <- sampleNames(K5)
pdf("K5_LUSC.pdf")
consensusmap(K5, hclustfun="average")
dev.off()

colnames(K6@consensus) <- sampleNames(K6)
rownames(K6@consensus) <- sampleNames(K6)
pdf("K6_LUSC.pdf")
consensusmap(K6, hclustfun="average")
dev.off()

colnames(K7@consensus) <- sampleNames(K7)
rownames(K7@consensus) <- sampleNames(K7)
pdf("K7_LUSC.pdf")
consensusmap(K7, hclustfun="average")
dev.off()

