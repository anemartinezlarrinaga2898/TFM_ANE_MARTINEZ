# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 27-04-2021

#OBJETIVO: es realizar el analisis de la expresion diferencial de mi grupo constante en LUAD, realmente no es un grupo como tal 
#ya que luego tiene varios integrantes extras pero estos no se mantienen contasnates. 
###########################################################################################

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/6_ScriptAnalisisDiferencial/")
library(tidyverse)

###########################################################################################
                      # G E N E R A C I O N  C O L D A T A: 
#1) Cargamos los datos y generamos las matrices: ·········································

GrupoConstante <- readRDS("GrupoConstante_Luad.rds")
Matriz <- readRDS("Luad_Filtrado.rds")

nombresMatriz <- colnames(Matriz)
nombresMatriz <- str_replace_all(nombresMatriz,"-",".")
colnames(Matriz) <- nombresMatriz

MatrizGrupoConstante <- Matriz[,GrupoConstante$item]
NombresMatrizGrupoConstnate <- colnames(MatrizGrupoConstante)

MatrizResto <- Matriz[,!colnames(Matriz) %in% colnames(MatrizGrupoConstante)]
NombresMatrizResto <- colnames(MatrizResto)

#2) Generamos el coldata: ·········································

GrupoConstante_Coldata <- data.frame(nombresPacientes=NombresMatrizGrupoConstnate,
                                     Grupo="Constante")
rownames(GrupoConstante_Coldata) <- NombresMatrizGrupoConstnate

GrupoResto_ColData <- data.frame(nombresPacientes=NombresMatrizResto,
                                 Grupo="NoConstante")
rownames(GrupoResto_ColData) <- NombresMatrizResto

ColData_GC <- rbind(GrupoConstante_Coldata,GrupoResto_ColData)

ColData_GC$Grupo <- as.factor(ColData_GC$Grupo)
ColData_GC$Grupo

ColData_GC$Grupo <- relevel(ColData_GC$Grupo, "NoConstante")
ColData_GC$Grupo

Matriz <- cbind(MatrizGrupoConstante,MatrizResto)

saveRDS(ColData_GC,"ColDataGrupoConstante.rds")
saveRDS(Matriz,"MatrizCuentasGrupoConstante.rds")

###################################################################################
                  # A N A L I S I S  D I F E R E N C I A L

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

Matriz <- Matriz
ColData <- ColData_GC

ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL

dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

DDS<- DESeq(dds)
DDS_Resultados <- results(DDS)
DDS_Resultados

DDS_Resultados$SYMBOL = mapIds(org.Hs.eg.db,
                    keys=rownames(DDS_Resultados),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

DDS_Resultados$GENENAME = mapIds(org.Hs.eg.db, 
                      keys=rownames(DDS_Resultados),
                      column="GENENAME",
                      keytype="ENSEMBL",
                      multiVals="first")

DDS_Resultados$PATH = mapIds(org.Hs.eg.db, 
                  keys=rownames(DDS_Resultados),
                  column="PATH",
                  keytype="ENSEMBL",
                  multiVals="first")

DDS_Resultados$ENTREZID = mapIds(org.Hs.eg.db, 
                      keys=rownames(DDS_Resultados),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

saveRDS(DDS,"ObjetoDeseq_Constante.rds")
saveRDS(DDS_Resultados,"Resultados_DeseqGrupoConstante.rds")

###################################################################################
# F I L T R A C I O N  R E S U L T A D O S: 

Resultados <- DDS_Resultados
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                     row.names=rownames(Resultados))
saveRDS(ResultadosDataFrame,"ResultadosGrupoConstante_DF.rds")

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
  filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))

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

