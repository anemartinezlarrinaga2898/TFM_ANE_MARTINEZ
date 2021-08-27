# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 27-05-2021

#OBJETIVO:Script realizar las anotaciones de mis genes
################################################################################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(biomaRt)

Genes <- readRDS("EstructuraLonger_OrdenGenes.rds")
ensembel <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

################################################################################################################################################

i <- 1
GenesAnotados_G1 <- list()
GenesAnotados_G2 <- list()
GenesAnotados_G3 <- list()
NombresLinea <- unique(Genes[[1]]$NombreLinea)
SecuenciaLineaCelular <- seq(1,length(NombresLinea),1)
for(i in seq_along(Genes)){
  Grupo <- Genes[[i]]
  
  j <- 1
  for (j in seq(SecuenciaLineaCelular)){
    Linea <- NombresLinea[j]
    
    LineaGrupo <- Grupo %>% filter(NombreLinea==Linea)
    
    GenesGrupo <- LineaGrupo$Gen
    
    AnotacionesGenes <- getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id","description"), 
                              filters = 'ensembl_gene_id', 
                              values = GenesGrupo, 
                              mart = ensembel)
    
    if (i==1){
      GenesAnotados_G1[[j]] <- AnotacionesGenes
      names(GenesAnotados_G1)[[j]] <- Linea
    }
    
    if (i==2){
      GenesAnotados_G2[[j]] <- AnotacionesGenes
      names(GenesAnotados_G2)[[j]] <- Linea
    }
    
    if (i==3){
      GenesAnotados_G3[[j]] <- AnotacionesGenes
      names(GenesAnotados_G3)[[j]] <- Linea
    }
    
  }
  
}

saveRDS(GenesAnotados_G1,"GenesAnotados_G1.rds")
saveRDS(GenesAnotados_G2,"GenesAnotados_G2.rds")
saveRDS(GenesAnotados_G3,"GenesAnotados_G3.rds")


