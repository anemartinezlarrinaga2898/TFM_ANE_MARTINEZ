# Resultados una vuelta LUAD: 

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/Objetos")

#Primero leemos los resultados
silhouette <- readRDS("res_nb_silhouette_luad.rds")
dindex <- readRDS("res_nb_dindex_luad.rds")
dunn <- readRDS("res_nb_dunn_luad.rds")
gap <- readRDS("res_nb_gap_luad.rds")
tau <- readRDS("res_nb_tau_luad.rds")

inicio <- 2
final <- 10
paso <- 1
numerosClusteres <- data.frame(x=seq(inicio,final,paso))

#Obtenemos lo que vamos a graficar: 

valoresSL <- data.frame(silhouette[[1]])
names(valoresSL) <- ("Silhouette")
valoresDI <- data.frame(dindex)
names(valoresDI) <- ("Dindex")
valoresDU <- data.frame(dunn[[1]])
names(valoresDU) <- ("Dunn")
valoresGAP <- data.frame(gap[[1]])
names(valoresGAP) <- ("Gap")
valoresTAU <- data.frame(tau[[1]]) 
names(valoresTAU) <- ("Tau")

indices <- cbind(valoresSL,valoresDI,valoresDU,valoresGAP,valoresTAU)

informacion <- cbind(numerosClusteres,indices)
names(informacion)[[1]] <- "NumerosdeClusteres"
informacion <- data.frame(informacion)

library(tidyverse)

indicesSilhouetteGrafico <- ggplot(informacion,aes(x=NumerosdeClusteres,y=Silhouette))+
  geom_line(size=1,color="#F886A8")+
  geom_point(data=informacion,aes(x=NumerosdeClusteres[[1]],y=Silhouette[[1]]),size=2)+
  scale_x_continuous(breaks = round(seq(min(informacion$NumerosdeClusteres), max(informacion$NumerosdeClusteres), by = 1),1))+
  ylab("Indice Silhouette")+
  xlab("Numero de Clusteres")
  
indicesDindexeGrafico <- ggplot(informacion,aes(x=NumerosdeClusteres,y=Dindex))+
  geom_line(size=1,color="#FE8D6F")+
  geom_point(data=informacion,aes(x=NumerosdeClusteres[[1]],y=Dindex[[1]]),size=2)+
  scale_x_continuous(breaks = round(seq(min(informacion$NumerosdeClusteres), max(informacion$NumerosdeClusteres), by = 1),1))+
  ylab("Indice Dindex")+
  xlab("Numero de Clusteres")

indicesDunnGrafico <- ggplot(informacion,aes(x=NumerosdeClusteres,y=Dunn))+
  geom_line(size=1,color="#FDC453")+
  geom_point(data=informacion,aes(x=NumerosdeClusteres[[6]],y=Dunn[[6]]),size=2)+
  scale_x_continuous(breaks = round(seq(min(informacion$NumerosdeClusteres), max(informacion$NumerosdeClusteres), by = 1),1))+
  ylab("Indice Dunn")+
  xlab("Numero de Clusteres")

indicesGapGrafico <- ggplot(informacion,aes(x=NumerosdeClusteres,y=Gap))+
  geom_line(size=1,color="#DFDD6C")+
  geom_point(data=informacion,aes(x=NumerosdeClusteres[[1]],y=Gap[[1]]),size=2)+
  scale_x_continuous(breaks = round(seq(min(informacion$NumerosdeClusteres), max(informacion$NumerosdeClusteres), by = 1),1))+
  ylab("Indice Gap")+
  xlab("Numero de Clusteres")

indicesTauGrafico <- ggplot(informacion,aes(x=NumerosdeClusteres,y=Tau))+
  geom_line(size=1,color="#A0DDE0")+
  geom_point(data=informacion,aes(x=NumerosdeClusteres[[1]],y=Tau[[1]]),size=2)+
  scale_x_continuous(breaks = round(seq(min(informacion$NumerosdeClusteres), max(informacion$NumerosdeClusteres), by = 1),1))+
  ylab("Indice Tau")+
  xlab("Numero de Clusteres")

library(ggpubr)

figura <- ggarrange(indicesDindexeGrafico, indicesDunnGrafico, indicesGapGrafico,indicesSilhouetteGrafico,indicesTauGrafico,
          labels = c("A", "B", "C","D","E"), ncol = 3, nrow = 2)

annotate_figure(figura,
                top = text_grob("Validacion Clusteres Kmeans LUAD", color = "Black", size = 14))
