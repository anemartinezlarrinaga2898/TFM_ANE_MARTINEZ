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

# Cargamos la matriz de cuentas de LUSC 

Lusc <- readRDS("Lusc_Filtrado.rds")
Lusc <- Lusc[-c(17875,17874,17873),]
rownames(Lusc) <- gsub("\\..*","",rownames(Lusc))

#Vamos a anotar los genes en vez de con el ensembl con el hugo symbol: 
NombresGenes <- data.frame(ensembl_gene_id=rownames(Lusc))

# Anotamos los 17872 genes que tiene mi matriz de LUSC: 
AnotacionesGenes <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                          filters="ensembl_gene_id",
                          values=NombresGenes,
                          mart=ensembel)

#Unimos las anotaciones con el nombre de mis genes de la matriz, para asegurar que tienen el mismo orden: 
Union_ConNaN <- full_join(NombresGenes,AnotacionesGenes,by="ensembl_gene_id") #Esta union tiene NaN pero tiene TODOS los genes

#Obtenemos los genes que tienen NaN
GenesConNan <-Union_ConNaN[which(is.na(Union_ConNaN$hgnc_symbol)),]

# Eliminamos los genes que tienen NaN
Union_SinNaN <- na.omit(Union_ConNaN) #Todos los genes pero SIN los que tiene NaN

#Ahora vamos a eliminar los genes de la matriz que tienen NaN: localizamos los indices de los genes. 

i <- 1
LocalizacionGen <- c()
for(i in seq(nrow(GenesConNan))){
  NombreGenConNan <- GenesConNan[i,1]
  
  LocalizacionGen[i] <- which(rownames(Lusc)==NombreGenConNan)
}

#Eliminamos los genes que tienen NaN en las anotaciones en la matriz original. 
Lusc <- Lusc[-LocalizacionGen,] #En esta matriz he eliminado los NaN 

NombresGenesSinNan <- data.frame(ensembl_gene_id=rownames(Lusc)) #Los nombres de los genes sin NaN

Union_FilasSinNan_Anotaciones <- full_join(NombresGenesSinNan,Union_SinNaN,by="ensembl_gene_id")


#Ahora encontramos los genes que estan duplicados en el gene id: 

#Tengo algo que me esta dando duplicado: ensembel 
#En estos NO ESTA EL DUPLICADO, son los UNICOS
Unicos_Ensembl <- Union_FilasSinNan_Anotaciones[!duplicated(Union_FilasSinNan_Anotaciones$ensembl_gene_id),]
NombresGenesDuplicados_Ensembl <- Union_FilasSinNan_Anotaciones[duplicated(Union_FilasSinNan_Anotaciones$ensembl_gene_id),]

Unicos_Hugo <- Union_FilasSinNan_Anotaciones[!duplicated(Union_FilasSinNan_Anotaciones$hgnc_symbol),]
NombresGenesDuplicados_Hugo <- Union_FilasSinNan_Anotaciones[duplicated(Union_FilasSinNan_Anotaciones$hgnc_symbol),]

#Tengo duplicados en HUGO que tengo que eliminar lo que hace que baje aun mas el numero de genes con el que estoy trabajando. 

Unicos_Hugo_Ensembl <- Unicos_Hugo[!duplicated(Unicos_Hugo$hgnc_symbol),]
Comprobacion <- Unicos_Hugo_Ensembl[duplicated(Unicos_Hugo_Ensembl$hgnc_symbol),]

AnotacionesFinales <- Unicos_Hugo_Ensembl

i <- 1

LuscMatriz_SinDuplicados <- matrix(0,nrow=nrow(AnotacionesFinales),ncol=ncol(Lusc))
colnames(LuscMatriz_SinDuplicados) <- colnames(Lusc)

for(i in seq(nrow(AnotacionesFinales))){
  Gen <- AnotacionesFinales[i,1]
  
  GenMatriz <- Lusc[which(rownames(Lusc)==Gen),]
  
  LuscMatriz_SinDuplicados[i,] <- as.matrix(GenMatriz)
 
}

rownames(LuscMatriz_SinDuplicados) <- AnotacionesFinales[,2]

saveRDS(LuscMatriz_SinDuplicados,"MatrizLusc_Genes_HugoSymbol.rds")



