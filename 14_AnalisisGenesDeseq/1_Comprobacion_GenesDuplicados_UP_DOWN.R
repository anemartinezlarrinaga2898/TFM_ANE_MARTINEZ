# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 25-06-2021

#OBJETIVO:Analizar los genes que se repiten en mis huellas e intentar ver si se encuentra UP o DOWN

###########################################################################################v###################################
library(tidyverse)
library(data.table)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/14_AnalisisGenesDeseq/")

################################################################################################################################

# ································ L U A D ·················································

H1 <- readRDS("ResultadosFiltrados_Luad_Comparacion1.rds")
H2 <- readRDS("ResultadosFiltrados_Luad_Comparacion2.rds")
H3 <- readRDS("ResultadosFiltrados_Luad_Comparacion3.rds")

H1 <- H1 %>% mutate(Grupo="Huella 1",ENSEMBL=rownames(H1))
H2 <- H2 %>% mutate(Grupo="Huella 2",ENSEMBL=rownames(H2))
H3 <- H3 %>% mutate(Grupo="Huella 3",ENSEMBL=rownames(H3))

Secuencia <- seq(1,482,1)

H <- rbind(H1,H2)
H <- rbind(H,H3)

rownames(H) <- Secuencia

Duplicados_Ensembl <- H[duplicated(H$ENSEMBL),]
Duplicados_Hugo <- H[duplicated(H$SYMBOL),]

Duplicados <- Duplicados_Ensembl$ENSEMBL

NombresColumnas <- colnames(H)

LocalizacionDuplicados <- list()
i <- 1

for(i in seq_along(Duplicados)){
  GenDuplicados <- Duplicados[i]
  
  X<- H[which(H$ENSEMBL==GenDuplicados),]
  
  LocalizacionDuplicados[[i]] <- X
}

i <- 1

D2G <- list()
D3G <- list()

for(i in seq_along(LocalizacionDuplicados)){
  
   X<- LocalizacionDuplicados[[i]]
  
  Dimensiones <- dim(X)
  
  if(Dimensiones[1]==2){D2G[[i]] <- X}
  
  if(Dimensiones[1]==3){D3G[[i]] <- X}
  
}

D2G <- D2G[lengths(D2G) != 0]
D3G <- D3G[lengths(D3G) != 0]

D2GT <- rbind(D2G[[1]],D2G[[2]])

i <- 1

for(i in seq(3,length(D2G),1)){
  X <- D2G[[i]]
  D2GT1 <- rbind(D2GT,X)
  D2GT <- D2GT1
}

saveRDS(D2GT1,"Duplicados_2grupos_Luad.rds")

D3GT <- rbind(D3G[[1]],D3G[[2]])

i <- 1

for(i in seq(3,length(D3G),1)){
  X <- D3G[[i]]
  D3GT1 <- rbind(D3GT,X)
  D3GT <- D3GT1
}

saveRDS(D3GT1,"Duplicados_3grupos_Luad.rds")

################################################################################################################################

# ································ L U S C ·················································

H1 <- readRDS("ResultadosFiltrados_Lusc_Comparacion1.rds")
H2 <- readRDS("ResultadosFiltrados_Lusc_Comparacion2.rds")
H3 <- readRDS("ResultadosFiltrados_Lusc_Comparacion3.rds")

H1 <- H1 %>% mutate(Grupo="Huella 1",ENSEMBL=rownames(H1))
H2 <- H2 %>% mutate(Grupo="Huella 2",ENSEMBL=rownames(H2))
H3 <- H3 %>% mutate(Grupo="Huella 3",ENSEMBL=rownames(H3))


H <- rbind(H1,H2)
H <- rbind(H,H3)

Secuencia <- seq(1,758,1)

rownames(H) <- Secuencia

Duplicados_Ensembl <- H[duplicated(H$ENSEMBL),]
Duplicados_Hugo <- H[duplicated(H$SYMBOL),]

Duplicados <- Duplicados_Ensembl$ENSEMBL

NombresColumnas <- colnames(H)

LocalizacionDuplicados <- list()
i <- 1

for(i in seq_along(Duplicados)){
  GenDuplicados <- Duplicados[i]
  
  X<- H[which(H$ENSEMBL==GenDuplicados),]
  
  LocalizacionDuplicados[[i]] <- X
}

i <- 1

D2G <- list()
D3G <- list()

for(i in seq_along(LocalizacionDuplicados)){
  
  X<- LocalizacionDuplicados[[i]]
  
  Dimensiones <- dim(X)
  
  if(Dimensiones[1]==2){D2G[[i]] <- X}
  
  if(Dimensiones[1]==3){D3G[[i]] <- X}
  
}

D2G <- D2G[lengths(D2G) != 0]
D3G <- D3G[lengths(D3G) != 0]

D2GT <- rbind(D2G[[1]],D2G[[2]])

i <- 1

for(i in seq(3,length(D2G),1)){
  X <- D2G[[i]]
  D2GT1 <- rbind(D2GT,X)
  D2GT <- D2GT1
}

saveRDS(D2GT1,"Duplicados_2grupos_Lusc.rds")

D3GT <- rbind(D3G[[1]],D3G[[2]])

i <- 1

for(i in seq(3,length(D3G),1)){
  X <- D3G[[i]]
  D3GT1 <- rbind(D3GT,X)
  D3GT <- D3GT1
}

saveRDS(D3GT1,"Duplicados_3grupos_Lusc.rds")

