# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 01-04-2021

#OBJETIVO:script para ver lo que hay dentro del objetio resultado generado por NMF

##################################################################################################

directoy <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")
library(tidyverse)

################################################################################################

#LUAD: 

matricesConsensosLuad <- readRDS("MatricesConsenso_Total_Luad.rds")
matrizWLuad <- readRDS("MatrizW_Total_Luad.rds")
matrizHLuad <- readRDS("MatrizH_Total_Luad.rds")

resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

resultadosNMFLuad_K2 <- resultadosNMFLuad[[1]]

########################################################################################################
#Codigo copiado de: https://github.com/naikai/sake/blob/master/R/nmf_utils.R

  #Extraer features: no funciona
nmf_extract_feature <- function(res, rawdata=NULL, manual.num=0, method="default", math="mad", FScutoff=0.9){
  feature.score <- featureScore(res)
  predict.feature <- predict(res, what="features", prob=T)
  
  data.feature <- data.frame(Gene=names(feature.score),
                             featureScore=feature.score,
                             Group=predict.feature$predict,
                             prob=predict.feature$prob,
                             stringsAsFactors=FALSE)
  if(method=="total"){
    print("return all featureScores")
  }else{
    if(method=="default"){
      # extracted features for each group
      if (manual.num==0){
        extract.feature <- extractFeatures(res)
      }else if (manual.num>0 && manual.num<=length(featureNames(res))){
        extract.feature <- extractFeatures(res, manual.num)
      }else{
        stop("wrong number of (manual num) features ")
      }
      
      data.feature <- extract.feature %>%
        lapply(., function(x) data.feature[x, ]) %>%
        rbindlist %>%
        as.data.frame.matrix
    }else if(method=="rank"){
      if(is.null(rawdata)){
        stop("error: need to provide original expression data if method is 'rank' ")
      }
      # data.feature <- cbind(data.feature, math=apply(log2(rawdata+1), 1, math)) %>%
      data.feature <- cbind(data.feature, math=apply(rawdata, 1, math)) %>%
        filter(featureScore>=FScutoff) %>%
        arrange(Group, dplyr::desc(math), dplyr::desc(prob))
      if(manual.num>0){
        data.feature <- group_by(data.feature, Group) %>%
          top_n(manual.num)
        # top_n(manual.num, math)
        # filter(min_rank(desc(math))<=manual.num)
      }
    }
  }
  return(data.feature)
}

  #Extraer grupo del sample: 
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

  #Para dibujar los plots: no me funciona por el grid
nmf_plot <- function(res, type="consensus", subsetRow=TRUE, save.image=F, hclustfun="average", silorder=F, add_original_name=T){
  if(save.image)
    pdf(res.pdf, width=18, height=15)
  
  if(type=="result"){
    print(plot(res))
  }else{
    si <- silhouette(res, what=type)
    
    if(type=="features"){
      if(silorder){
        basismap(res, Rowv = si, subsetRow=subsetRow)
      }else{
        basismap(res, subsetRow = subsetRow)
      }
    }else if(type=="samples"){
      if(silorder){
        coefmap(res, Colv = si)
      }else{
        coefmap(res)
      }
    }else if(type=="consensus"){
      if(add_original_name){
        colnames(res@consensus) <- sampleNames(res)
        rownames(res@consensus) <- sampleNames(res)
      }
      consensusmap(res, hclustfun=hclustfun)
    }
  }
  if(save.image)
    dev.off()
}

  #Para dibujar el silhouette plot:  
nmf_silhouette_plot <- function(res, type="consensus", silorder=F){
  si <- silhouette(res, what=type)
  plot(si)
}

  #Para elegir el mejor plot: no me funciona
nmf_select_best_k <- function(res, type="consensus"){
  data <- NULL
  
  if(type=="consensus"){
    predict.consensus <- predict(res, what="consensus")
    data <- data.frame(Sample_ID=sampleNames(res),
                       nmf_subtypes = predict.consensus)
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if(matchConseOrder){
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  }else if(type=="samples"){
    predict.samples <- predict(res, what="samples", prob=T)
    data <- data.frame(Sample_ID=names(predict.samples$predict),
                       nmf_subtypes = predict.samples$predict,
                       prob = predict.samples$prob)
  }else{
    stop(paste("Wrong type:", type))
  }
  return(data)
}

########################################################################################################
  #Esta funcion me devuelve a que grupo pertenece cada sample!
SamplesNMFLuad <- nmf_extract_group(resultadosNMFLuad_K2)
FeaturesNMFLuad <- nmf_extract_feature(resultadosNMFLuad_K2)

#Graficos NMF Luad para K=2; 
  #Dibujar el consesus map con los nombres
colnames(resultadosNMFLuad_K2@consensus) <- sampleNames(resultadosNMFLuad_K2)
rownames(resultadosNMFLuad_K2@consensus) <- sampleNames(resultadosNMFLuad_K2)
pdf("PruebaConsensus.pdf")
consensusmap(resultadosNMFLuad_K2, hclustfun="average")
dev.off()

nmf_silhouette_plot(resultadosNMFLuad_K2)
