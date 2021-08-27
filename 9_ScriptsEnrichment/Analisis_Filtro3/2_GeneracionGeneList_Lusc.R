# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-06-2021

#OBJETIVO: Generacion de las geneList necesarias para hacer el GSEA. 

################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)

################################################################################################################################
#··································  L U A D ··········································

# GENE LIST de las compraciones TOTALES 

#RC1T:  Resultados comparacion 1 totales: 

RC1T <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion1_Lusc_DF.rds")
RC2T <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion2_Lusc_DF.rds")
RC3T <- readRDS("Analisis_Filtro3/DatosPartida/ResultadosComparacion3_Lusc_DF.rds")

# C O M P A R A C I O N  1: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC1T$log2FoldChange
NG <- rownames(RC1T)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLT1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLT1_Lusc.rds")

# C O M P A R A C I O N  2: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC2T$log2FoldChange
NG <- rownames(RC2T)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLT1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLT2_Lusc.rds")

# C O M P A R A C I O N  3: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC3T$log2FoldChange
NG <- rownames(RC3T)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLT1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLT3_Lusc.rds")

#···············································································

# GENE LIST de las compraciones FILTRADAS: 

RC1F <- readRDS("Analisis_Filtro3/ResultadosFiltrados/RFC1_Lusc.rds")
RC2F <- readRDS("Analisis_Filtro3/ResultadosFiltrados/RFC2_Lusc.rds")
RC3F <- readRDS("Analisis_Filtro3/ResultadosFiltrados/RFC3_Lusc.rds")

# C O M P A R A C I O N  1: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC1F$log2FoldChange
NG <- rownames(RC1F)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLT1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLF1_Lusc.rds")

# C O M P A R A C I O N  2: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC2F$log2FoldChange
NG <- rownames(RC2F)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLT1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLF2_Lusc.rds")

# C O M P A R A C I O N  3: ---------> T O T A L E S

#VLFC: Valores log fold change 
# NGT: Nombres Genes

VLFC <- RC3F$log2FoldChange
NG <- rownames(RC3F)

# GL es el data fram de gene List 
# gl es la lista en si

GL <- data.frame(GeneId=NG,LogFoldChange=VLFC)
gl <-  GL[,2]
names(gl) <-  as.character(GL[,1])
gl <-  sort(gl, decreasing = TRUE)

#GLF1: Gene list totales 1 

saveRDS(gl,"Analisis_Filtro3/GLF3_Lusc.rds")

