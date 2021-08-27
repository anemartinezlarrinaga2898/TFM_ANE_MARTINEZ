# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-06-2021

#OBJETIVO:Analisis del enrichment de LUAD y LUSC con los diferentes hallmarks con el filtro 1.5

##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

ObtencionCuentas <- function(x){
  Separacion <- str_replace_all(x,"/"," ")
  Cuentas <- length(unlist(strsplit(Separacion," ")))
  return(Cuentas)
}

H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)

C2 <-  msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene)

C5 <-  msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene)

C6 <-  msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, ensembl_gene)

C7 <-  msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, ensembl_gene)

##############################################################################################################


# ···················································· L U A D ············································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF1_Luad.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT1_Luad.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Luad_1_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Luad_1_Filtrado_Hallamarks.rds")

#··············································································································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF2_Luad.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT2_Luad.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Luad_2_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Luad_2_Filtrado_Hallamarks.rds")

#··············································································································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF3_Luad.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT3_Luad.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Luad_3_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Luad_3_Filtrado_Hallamarks.rds")

##############################################################################################################################

# ···················································· L U S C ············································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF1_Lusc.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT1_Lusc.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Lusc_1_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_1_Filtrado_Hallamarks.rds")

#··············································································································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF2_Lusc.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT2_Lusc.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Lusc_2_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_2_Filtrado_Hallamarks.rds")

#··············································································································

geneList_Filtrado <- readRDS("Analisis_Filtro3/Objetos_GSEA/FILTRADO/GLF3_Lusc.rds")
geneList_Total <- readRDS("Analisis_Filtro3/Objetos_GSEA/TOTAL/GLT3_Lusc.rds")
condiciones <- c("H","C2","C5","C6","C7")
geneList <- list(geneList_Filtrado,geneList_Total)

ValoresGsea <- list(H,C2,C5,C6,C7)
ResultadosGsea_Filtrado <- list()
ResultadosGsea_Total <- list()

i <- 1
j <- 1
for (j in seq_along(geneList)){
  listaGenes <- geneList[[j]]
  print(paste("j",j))
  if(i==5){i <- 1}
  for (i in seq_along(ValoresGsea)){
    TerminoGsea <- ValoresGsea[[i]]
    print(paste("i",i))
    gseaGO <- GSEA(listaGenes, 
                   TERM2GENE=TerminoGsea, 
                   minGSSize = 10,
                   pvalueCutoff = 0.05,
                   verbose=FALSE,
                   eps=0)
    if (j==1){ResultadosGsea_Filtrado[[i]] <-gseaGO }
    else(ResultadosGsea_Total[[i]] <-gseaGO)
  }
}

names(ResultadosGsea_Filtrado) <- condiciones
names(ResultadosGsea_Total) <- condiciones

saveRDS(ResultadosGsea_Total,"Analisis_Filtro3/Resultados_GSEA/Total/AnalisisGSEA_Lusc_3_Total_Hallamarks.rds")
saveRDS(ResultadosGsea_Filtrado,"Analisis_Filtro3/Resultados_GSEA/Filtrado/AnalisisGSEA_Lusc_3_Filtrado_Hallamarks.rds")



