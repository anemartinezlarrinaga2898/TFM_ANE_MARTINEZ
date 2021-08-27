# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-06-2021

#OBJETIVO: Obtencion de la matriz con los genes en formato gene symbol. 
###########################################################################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/13_GeneRegulatoryNetwork/")
library(tidyverse)

library(biomaRt)

ensembel <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
anotaciones <- listAttributes(ensembel)
filtros <- listFilters(ensembel)

###########################################################################################################################################################
Luad <- readRDS("Luad_Filtrado.rds")
rownames(Luad) <- gsub("\\..*","",rownames(Luad))
Luad <- Luad[-c(17789,17788,17787),]


#Vamos a anotar los genes en vez de con el ensembl con el hugo symbol: 
NombresGenes <- data.frame(ensembl_gene_id=rownames(Luad))

# Anotamos los 17872 genes que tiene mi matriz de LUSC: 
AnotacionesGenes <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                          filters="ensembl_gene_id",
                          values=NombresGenes,
                          mart=ensembel)

Union_ConNaN <- full_join(NombresGenes,AnotacionesGenes,by="ensembl_gene_id") #Esta union tiene NaN pero tiene TODOS los genes

#Obtenemos los genes que tienen NaN
GenesConNan <-Union_ConNaN[which(is.na(Union_ConNaN$hgnc_symbol)),]

# Eliminamos los genes que tienen NaN
Union_SinNaN <- na.omit(Union_ConNaN) #Todos los genes pero SIN los que tiene NaN

i <- 1
LocalizacionGen <- c()
for(i in seq(nrow(GenesConNan))){
  NombreGenConNan <- GenesConNan[i,1]
  
  LocalizacionGen[i] <- which(rownames(Luad)==NombreGenConNan)
}

#Eliminamos los genes que tienen NaN en las anotaciones en la matriz original. 
Luad <- Luad[-LocalizacionGen,]

NombresGenesSinNan <- data.frame(ensembl_gene_id=rownames(Luad)) #Los nombres de los genes sin NaN

Union_FilasSinNan_Anotaciones <- full_join(NombresGenesSinNan,Union_SinNaN,by="ensembl_gene_id")

Unicos_Ensembl <- Union_FilasSinNan_Anotaciones[!duplicated(Union_FilasSinNan_Anotaciones$ensembl_gene_id),]
NombresGenesDuplicados_Ensembl <- Union_FilasSinNan_Anotaciones[duplicated(Union_FilasSinNan_Anotaciones$ensembl_gene_id),]

Unicos_Hugo <- Union_FilasSinNan_Anotaciones[!duplicated(Union_FilasSinNan_Anotaciones$hgnc_symbol),]
NombresGenesDuplicados_Hugo <- Union_FilasSinNan_Anotaciones[duplicated(Union_FilasSinNan_Anotaciones$hgnc_symbol),]

#Tengo duplicados en HUGO que tengo que eliminar lo que hace que baje aun mas el numero de genes con el que estoy trabajando. 

Unicos_Hugo_Ensembl <- Unicos_Hugo[!duplicated(Unicos_Hugo$hgnc_symbol),]
Comprobacion <- Unicos_Hugo_Ensembl[duplicated(Unicos_Hugo_Ensembl$hgnc_symbol),]

AnotacionesFinales <- Unicos_Hugo_Ensembl

i <- 1

LuscMatriz_SinDuplicados <- matrix(0,nrow=nrow(AnotacionesFinales),ncol=ncol(Luad))
colnames(LuscMatriz_SinDuplicados) <- colnames(Luad)

for(i in seq(nrow(AnotacionesFinales))){
  Gen <- AnotacionesFinales[i,1]
  
  GenMatriz <- Luad[which(rownames(Luad)==Gen),]
  
  LuscMatriz_SinDuplicados[i,] <- as.matrix(GenMatriz)
  
}

rownames(LuscMatriz_SinDuplicados) <- AnotacionesFinales[,2]

saveRDS(LuscMatriz_SinDuplicados,"MatrizLuad_Genes_HugoSymbol.rds")

