Matriz <- readRDS("Comparacion2_Lusc.rds")
ColData <- readRDS("Comparacion2ColData_Lusc.rds")
ENSEMBEL <- gsub("\\..*","",row.names(Matriz))
rownames(Matriz) <- ENSEMBEL
dds <- DESeqDataSetFromMatrix(countData = Matriz,
                              colData= ColData,
                              design = ~ Grupo)

#DDS <- DESeq(dds)

DDS <- readRDS("AnalisisExpresion_2_Lusc.rds")
Resultados <- results(DDS)

Resultados$SYMBOL = mapIds(org.Hs.eg.db,
                           keys=rownames(Resultados),
                           column="SYMBOL",
                           keytype="ENSEMBL",
                           multiVals="first")

Resultados <- readRDS("ResultadosAnalasisExpresion_2_Lusc.rds")
ResultadosDataFrame <- as.data.frame(dplyr::mutate(as.data.frame(Resultados)),
                                   row.names=rownames(Resultados))
library(dplyr)
Filtracion2 <-ResultadosDataFrame %>%
  dplyr::filter(padj<0.05 & (log2FoldChange >2 | log2FoldChange < -2 ))


NanEnResultados <- Resultados[which(is.na(Resultados$SYMBOL)),]
NanEnResultados <- ResultadosDataFrame[which(is.na(ResultadosDataFrame$SYMBOL)),]

BiocManager::install("EnsDb.Hsapiens.v75")
library("EnsDb.Hsapiens.v75")

Resultados$SYMBOL  <- mapIds(EnsDb.Hsapiens.v75, 
                     keys=rownames(Resultados), 
                     column="SYMBOL", 
                     keytype="GENEID", 
                     multiVals="first")
saveRDS(Resultados,"ResultadosAnalasisExpresion_2_Lusc.rds")


BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
library("BSgenome.Hsapiens.NCBI.GRCh38")

Resultados$SYMBOL  <- mapIds(BSgenome.Hsapiens.NCBI.GRCh38, 
                             keys=rownames(Resultados), 
                             column="SYMBOL", 
                             keytype="GENEID", 
                             multiVals="first")


