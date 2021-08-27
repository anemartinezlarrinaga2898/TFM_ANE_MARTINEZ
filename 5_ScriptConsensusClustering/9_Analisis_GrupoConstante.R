# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 16-04-2021


#OBJETIVO: analisis de los 12 pacientes que se agrupan de forma conjunta en LUAD en KMEANS tambien se agrupan de forma conjunta en NMF. 

######################################################################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(NMF)
library(tidyverse)

#####################################################################################################################################

# Primero genero los grupos de cada valor de K: 

resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

nmf_extract_group <- function(res, type="consensus", matchConseOrder=F){
  data <- NULL
  if(type=="consensus"){
    predict.consensus <- predict(res, what="consensus")
    silhouette.consensus <- silhouette(res, what="consensus")
    # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
    # that matches the original sampleNames(res) order
    # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
    # order of the samples in consensusmap(res). It is just for displaying
    # Therefore, the merged data frame sampleNames(res) + a.predict.consensus is the final
    # consensus results.
    data <- data.frame(Sample_ID=sampleNames(res),
                       nmf_subtypes = predict.consensus,
                       sil_width = signif(silhouette.consensus[, "sil_width"], 3))
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if(matchConseOrder){
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  }else if(type=="samples"){
    predict.samples <- predict(res, what="samples", prob=T)
    silhouette.samples <- silhouette(res, what="samples")
    data <- data.frame(Sample_ID=names(predict.samples$predict),
                       nmf_subtypes = predict.samples$predict,
                       sil_width = signif(silhouette.samples[, "sil_width"], 3),
                       prob = signif(predict.samples$prob, 3))
  }else{
    stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
  }
  return(data)
}

K2 <- resultadosNMFLuad[[1]]
K3 <- resultadosNMFLuad[[2]]
K4 <- resultadosNMFLuad[[3]]
K5 <- resultadosNMFLuad[[4]]
K6 <- resultadosNMFLuad[[5]]
K7 <- resultadosNMFLuad[[6]]

K2GLuad <- nmf_extract_group(K2)
IDK2<- str_replace_all(K2GLuad$Sample_ID,"-",".")
K2GLuad <- cbind(K2GLuad,IDK2)
K2GLuad <- K2GLuad[,-1]
colnames(K2GLuad)[3] <- "Sample_ID"

K3GLuad <- nmf_extract_group(K3)
IDK3<- str_replace_all(K3GLuad$Sample_ID,"-",".")
K3GLuad <- cbind(K3GLuad,IDK3)
K3GLuad <- K3GLuad[,-1]
colnames(K3GLuad)[3] <- "Sample_ID"

K4GLuad <- nmf_extract_group(K4)
IDK4<- str_replace_all(K4GLuad$Sample_ID,"-",".")
K4GLuad <- cbind(K4GLuad,IDK4)
K4GLuad <- K4GLuad[,-1]
colnames(K4GLuad)[3] <- "Sample_ID"

K5GLuad <- nmf_extract_group(K5)
IDK5<- str_replace_all(K5GLuad$Sample_ID,"-",".")
K5GLuad <- cbind(K5GLuad,IDK5)
K5GLuad <- K5GLuad[,-1]
colnames(K5GLuad)[3] <- "Sample_ID"

K6GLuad <- nmf_extract_group(K6)
IDK6<- str_replace_all(K6GLuad$Sample_ID,"-",".")
K6GLuad <- cbind(K6GLuad,IDK6)
K6GLuad <- K6GLuad[,-1]
colnames(K6GLuad)[3] <- "Sample_ID"

K7GLuad <- nmf_extract_group(K7)
IDK7<- str_replace_all(K7GLuad$Sample_ID,"-",".")
K7GLuad <- cbind(K7GLuad,IDK7)
K7GLuad <- K7GLuad[,-1]
colnames(K7GLuad)[3] <- "Sample_ID"

GrupoConstante <- readRDS("GrupoConstante_Luad.rds")
colnames(GrupoConstante) <- c("K","Cluster","Sample_ID","ItemConsensus")
GrupoConstante <- GrupoConstante[,-1]

K2_GC_Luad <- inner_join(GrupoConstante,K2GLuad,by="Sample_ID")
K2_GC_Luad <- K2_GC_Luad[,c(2,1,3,4,5)]
colnames(K2_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")

K3_GC_Luad <- inner_join(GrupoConstante,K3GLuad,by="Sample_ID")
K3_GC_Luad <- K3_GC_Luad[,c(2,1,3,4,5)]
colnames(K3_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")

K4_GC_Luad <- inner_join(GrupoConstante,K4GLuad,by="Sample_ID") #Aqui no me lo agrupa igual
K4_GC_Luad <- K4_GC_Luad[,c(2,1,3,4,5)]
colnames(K4_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")

K5_GC_Luad <- inner_join(GrupoConstante,K5GLuad,by="Sample_ID")
K5_GC_Luad <- K5_GC_Luad[,c(2,1,3,4,5)]
colnames(K5_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")

K6_GC_Luad <- inner_join(GrupoConstante,K6GLuad,by="Sample_ID")
K6_GC_Luad <- K6_GC_Luad[,c(2,1,3,4,5)]
colnames(K6_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")

K7_GC_Luad <- inner_join(GrupoConstante,K7GLuad,by="Sample_ID")
K7_GC_Luad <- K7_GC_Luad[,c(2,1,3,4,5)]
colnames(K7_GC_Luad) <- c("Sample_ID","Kmeans_Cluster","ConsensusKmeans","NMF_Cluster","Consesus_NMF")
