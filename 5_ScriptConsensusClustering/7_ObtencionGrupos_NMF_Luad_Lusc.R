# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 17-04-2021

#SCRIPT para obtener los grupos de NMF:

#####################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
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
library(NMF)
library(tidyverse)

#####################################################################################################################################

# ······················· L U A D ·············································

resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

K2_NMF_Luad <- resultadosNMFLuad[[1]]
K3_NMF_Luad <- resultadosNMFLuad[[2]]
K4_NMF_Luad <- resultadosNMFLuad[[3]]
K5_NMF_Luad <- resultadosNMFLuad[[4]]
K6_NMF_Luad <- resultadosNMFLuad[[5]]
K7_NMF_Luad <- resultadosNMFLuad[[6]]

K2_NMF_Luad <- nmf_extract_group(K2_NMF_Luad)
K3_NMF_Luad <- nmf_extract_group(K3_NMF_Luad)
K4_NMF_Luad <- nmf_extract_group(K4_NMF_Luad)
K5_NMF_Luad <- nmf_extract_group(K5_NMF_Luad)
K6_NMF_Luad <- nmf_extract_group(K6_NMF_Luad)
K7_NMF_Luad <- nmf_extract_group(K7_NMF_Luad)

saveRDS(K2_NMF_Luad,"K2_NMF_Luad.rds")
saveRDS(K3_NMF_Luad,"K3_NMF_Luad.rds")
saveRDS(K4_NMF_Luad,"K4_NMF_Luad.rds")
saveRDS(K5_NMF_Luad,"K5_NMF_Luad.rds")
saveRDS(K6_NMF_Luad,"K6_NMF_Luad.rds")
saveRDS(K7_NMF_Luad,"K7_NMF_Luad.rds")

K3.1_NMF_Luad <- K3_NMF_Luad %>% filter(nmf_subtypes==1)
K3.2_NMF_Luad <- K3_NMF_Luad %>% filter(nmf_subtypes==2)
K3.3_NMF_Luad <- K3_NMF_Luad %>% filter(nmf_subtypes==3)

saveRDS(K3.1_NMF_Luad,"K3.1_NMF_Luad.rds")
saveRDS(K3.2_NMF_Luad,"K3.2_NMF_Luad.rds")
saveRDS(K3.3_NMF_Luad,"K3.3_NMF_Luad.rds")

#####################################################################################################################################

resultadosNMFLusc <- readRDS("Resultados_NMF_Total_Lusc.rds")

K2_NMF_Lusc <- resultadosNMFLusc[[1]]
K3_NMF_Lusc <- resultadosNMFLusc[[2]]
K4_NMF_Lusc <- resultadosNMFLusc[[3]]
K5_NMF_Lusc <- resultadosNMFLusc[[4]]
K6_NMF_Lusc <- resultadosNMFLusc[[5]]
K7_NMF_Lusc <- resultadosNMFLusc[[6]]

K2_NMF_Lusc <- nmf_extract_group(K2_NMF_Lusc)
K3_NMF_Lusc <- nmf_extract_group(K3_NMF_Lusc)
K4_NMF_Lusc <- nmf_extract_group(K4_NMF_Lusc)
K5_NMF_Lusc <- nmf_extract_group(K5_NMF_Lusc)
K6_NMF_Lusc <- nmf_extract_group(K6_NMF_Lusc)
K7_NMF_Lusc <- nmf_extract_group(K7_NMF_Lusc)

saveRDS(K2_NMF_Lusc,"K2_NMF_Lusc.rds")
saveRDS(K3_NMF_Lusc,"K3_NMF_Lusc.rds")
saveRDS(K4_NMF_Lusc,"K4_NMF_Lusc.rds")
saveRDS(K5_NMF_Lusc,"K5_NMF_Lusc.rds")
saveRDS(K6_NMF_Lusc,"K6_NMF_Lusc.rds")
saveRDS(K7_NMF_Lusc,"K7_NMF_Lusc.rds")

K3.1_NMF_Lusc <- K3_NMF_Lusc %>% filter(nmf_subtypes==1)
K3.2_NMF_Lusc <- K3_NMF_Lusc %>% filter(nmf_subtypes==2)
K3.3_NMF_Lusc <- K3_NMF_Lusc %>% filter(nmf_subtypes==3)

saveRDS(K3.1_NMF_Lusc,"K3.1_NMF_Lusc.rds")
saveRDS(K3.2_NMF_Lusc,"K3.2_NMF_Lusc.rds")
saveRDS(K3.3_NMF_Lusc,"K3.3_NMF_Lusc.rds")

#####################################################################################################################################





