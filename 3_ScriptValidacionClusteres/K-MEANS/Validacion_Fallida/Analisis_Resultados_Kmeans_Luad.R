directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLuad/Objetos")

library(tidyverse)
library(NbClust)
library(cluster)
library(ggpubr)
library(factoextra)

resultadosKmeans <- readRDS("Resultados_Kmeans_ParaCadaK_luad.rds")
muestras <- data.frame(readRDS('MuestrasLuadTraspuesta.rds'))


kmeans_2 <- fviz_cluster(resultadosKmeans[[1]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_3 <- fviz_cluster(resultadosKmeans[[2]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_4 <- fviz_cluster(resultadosKmeans[[3]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_5 <- fviz_cluster(resultadosKmeans[[4]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_6 <- fviz_cluster(resultadosKmeans[[5]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_7 <- fviz_cluster(resultadosKmeans[[6]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_8 <- fviz_cluster(resultadosKmeans[[7]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_9 <- fviz_cluster(resultadosKmeans[[8]], data = muestras,
                         geom = "point",
                         ellipse.type = "convex", 
                         ggtheme = theme_bw()
)

kmeans_10 <- fviz_cluster(resultadosKmeans[[9]], data = muestras,
                          geom = "point",
                          ellipse.type = "convex", 
                          ggtheme = theme_bw()
)

figura <- ggarrange(kmeans_2,kmeans_3,kmeans_4,kmeans_5,kmeans_6,kmeans_7,kmeans_8,kmeans_9,kmeans_10,
                    labels = c("K=2", "K=3", "K=4","K=5","K=6","K=7","K=8","K=9","K=10"), ncol = 5, nrow = 2)
annotate_figure(figura,
                top = text_grob("Validacion Clusteres Kmeans LUAD", color = "Black", size = 14))


