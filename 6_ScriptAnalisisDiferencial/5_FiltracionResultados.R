# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-04-2021

#OBJETIVO: es para probar que tipo de filtracion es la mejor para obtener huellas mas pequeñas
###########################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)


###############################################################################################

#······························· L U A D ····································

# 1) C O M P A R A C I O N  1:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_1_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05)
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion1_Luad_DF.rds")

Filtracion0 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >0 | log2FoldChange < 0 )) 

GenesUp <- Filtracion0 %>% filter(log2FoldChange>0)
nrow(GenesUp)
GenesDown<- Filtracion0 %>% filter(log2FoldChange< 0)
nrow(GenesDown)

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion1_Luad.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Luad_Comparacion1.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)

# 2) C O M P A R A C I O N  2:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_2_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05)
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion2_Luad_DF.rds")

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion2_Luad.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Luad_Comparacion2.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)

# 3) C O M P A R A C I O N  3:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_3_Luad.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05)
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion3_Luad_DF.rds")

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion3_Luad.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Luad_Comparacion3.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)

###############################################################################################

#······························· L U S C ····································

# 1) C O M P A R A C I O N  1:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_1_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05) 
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion1_Lusc_DF.rds")

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion1_Lusc.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Lusc_Comparacion1.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)

# 2) C O M P A R A C I O N  2:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05)
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion2_Lusc_DF.rds")

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion2_Lusc.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Lusc_Comparacion2.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)

# 3) C O M P A R A C I O N  3:··············································

Resultados <- readRDS("ResultadosAnalasisExpresion_3_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
ResultadosDataFrame <- ResultadosDataFrame %>% filter(padj<0.05)
saveRDS(ResultadosDataFrame,"ResultadosDF/ResultadosComparacion3_Lusc_DF.rds")

Filtracion1 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1 | log2FoldChange < -1 )) 

GenesUp <- Filtracion1 %>% filter(log2FoldChange>1)
nrow(GenesUp)
GenesDown<- Filtracion1 %>% filter(log2FoldChange< -1)
nrow(GenesDown)

Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

write.csv(Filtracion2, file = "Comparacion3_Lusc.csv", row.names = TRUE)
saveRDS(Filtracion2,"ResultadosFiltrados_Lusc_Comparacion3.rds")

GenesUp <- Filtracion2 %>% filter(log2FoldChange>2)
nrow(GenesUp)
GenesDown<- Filtracion2 %>% filter(log2FoldChange< -2)
nrow(GenesDown)

Filtracion3 <-ResultadosDataFrame %>%
  filter(padj<0.05 & (log2FoldChange >1.5 | log2FoldChange < -1.5 ))

GenesUp <- Filtracion3 %>% filter(log2FoldChange>1.5)
nrow(GenesUp)
GenesDown<- Filtracion3 %>% filter(log2FoldChange< -1.5)
nrow(GenesDown)
