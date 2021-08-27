# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 06-04-2021


#OBJETIVO: representar los resultados del consesusClustering, obtenido la asginacion para cada cluster. 
  #1) ANALISIS DE LOS RESULTADOS DE NMF

########################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/5_ScriptConsensusClustering/")

#########################################################################################################

# L U A D: 

matricesConsensosLuad <- readRDS("MatricesConsenso_Total_Luad.rds")
resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

#········································································
#Codigo copiado de: https://github.com/naikai/sake/blob/master/R/nmf_utils.R
#········································································

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

# K=3 y K=5;

K3 <- resultadosNMFLuad[[2]]
K5 <- resultadosNMFLuad[[4]]

K3_Grupos_Luad <- nmf_extract_group(K3)
K5_Grupos_Luad <- nmf_extract_group(K5)

K3_Grupos_Luad_Ordenados <- K3_Grupos_Luad[order(K3_Grupos_Luad$nmf_subtypes),]
K5_Grupos_Luad_Ordenados <- K5_Grupos_Luad[order(K5_Grupos_Luad$nmf_subtypes),]

#########################################################################################################

# L U S C: 

matricesConsensosLusc <- readRDS("MatricesConsenso_Total_Lusc.rds")


resultadosNMFLusc <- readRDS("Resultados_NMF_Total_Lusc.rds")

# K=3 y K=5;

K3 <- resultadosNMFLusc[[2]]
K5 <- resultadosNMFLusc[[4]]

K3_Grupos_Lusc <- nmf_extract_group(K3)
K5_Grupos_Lusc <- nmf_extract_group(K5)

K3_Grupos_Lusc_Ordenados <- K3_Grupos_Lusc[order(K3_Grupos_Lusc$nmf_subtypes),]
K5_Grupos_Lusc_Ordenados <- K5_Grupos_Lusc[order(K5_Grupos_Lusc$nmf_subtypes),]

#########################################################################################################

#Comparar K3 y K5 de luad, ver si el grupo 3 hace split en el grupo 3 y 4 en el caso de k=5:

  #Primero dentro de todos los resultados de K3 y K5 filtramos por el numero de K=3 y K=4
Grupo_K3_Luad_K3 <- K3_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==3)
Grupo_K5_Luad_K3 <- K5_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==3)
Grupo_K5_Luad_K4 <- K5_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==4)

#En total en K3 hay 25 individuos mientras que en los grupos K5 (K5.3 y K5.4) tienen un total del 23

  #Organizamos el ID de los pacientes y creamos un data.frame para luego poder unir ambos. 
Grupo_K3_Luad_Pacientes <-Grupo_K3_Luad_K3$Sample_ID
Grupo_K3_Luad_Pacientes <- Grupo_K3_Luad_Pacientes[order(Grupo_K3_Luad_Pacientes)]
Grupo_K3_Luad_Pacientes_DF <- data.frame(NumeroPaciente=seq(1,25,1),ID=Grupo_K3_Luad_Pacientes)

Grupo_K5.3_Luad_Pacientes <-Grupo_K5_Luad_K3$Sample_ID
Grupo_K5.4_Luad_Pacientes <-Grupo_K5_Luad_K4$Sample_ID
  #En este caso tenemos que unirlos ambos porque tengo dos y luego lo ordeno. (Podriamos hbaer filtrado por igual a 3 y 4)
Grupos_K5.3.4_Luad_Pacientes <- c(Grupo_K5.3_Luad_Pacientes,Grupo_K5.4_Luad_Pacientes)
Grupos_K5.3.4_Luad_Pacientes <- Grupos_K5.3.4_Luad_Pacientes[order(Grupos_K5.3.4_Luad_Pacientes)]
Grupos_K5.3.4_Luad_Pacientes_DF <- data.frame(NumeroPaciente=seq(1,23,1),ID=Grupos_K5.3.4_Luad_Pacientes)

PacientesComunes <- inner_join(Grupo_K3_Luad_Pacientes_DF,Grupos_K5.3.4_Luad_Pacientes_DF,by="ID")
colnames(PacientesComunes) <- c("K3","ID","K5.3.4")

PacientesUnicosdeK3 <- setdiff(Grupo_K3_Luad_Pacientes,Grupos_K5.3.4_Luad_Pacientes)
PacientesUnicosdeK3 <- data.frame(NumeroPaciente=seq(1,length(PacientesUnicosdeK3),1),
                                  ID=PacientesUnicosdeK3)

PacientesUnicosdeK5 <- setdiff(Grupos_K5.3.4_Luad_Pacientes,Grupo_K3_Luad_Pacientes)
PacientesUnicosdeK5 <- data.frame(NumeroPaciente=seq(1,length(PacientesUnicosdeK5),1),
                                  ID=PacientesUnicosdeK5)

#########################################################################################################

#GRUPO 1 en K=3:

Grupo_K3_Luad_K1 <- K3_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==1)
Grupo_K3_Luad_Pacientes_K1 <-Grupo_K3_Luad_K1$Sample_ID

Grupo_K3_Luad_Pacientes_K1 <- Grupo_K3_Luad_Pacientes_K1[order(Grupo_K3_Luad_Pacientes_K1)]
Grupo_K3_Luad_Pacientes_K1_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K3_Luad_Pacientes_K1),1),
                                            Sample_ID=Grupo_K3_Luad_Pacientes_K1)

#GRUPO 1 Y 2 en K=5

Grupo_K5_Luad_K1.2 <- K5_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==1|nmf_subtypes==2)
Grupo_K5_Luad_Pacientes_K1.2 <-Grupo_K5_Luad_K1.2$Sample_ID

Grupo_K5_Luad_Pacientes_K1.2 <- Grupo_K5_Luad_Pacientes_K1.2[order(Grupo_K5_Luad_Pacientes_K1.2)]
Grupo_K5_Luad_Pacientes_K1.2_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K5_Luad_Pacientes_K1.2),1),
                                            Sample_ID=Grupo_K5_Luad_Pacientes_K1.2)


#Vamos a buscar los pacientes del grupo 1 en k=3 donde se encuentran en k=5?
Grupo5_Luad <- K5_Grupos_Luad_Ordenados[,c(1,2)]
Pacientes_K3.1_K5 <- inner_join(Grupo_K3_Luad_Pacientes_K1_DF,K5_Grupos_Luad_Ordenados,by="Sample_ID")

#Comprobacion: 
Grupo5_Luad_Comprobacion <- inner_join(Grupo_K3_Luad_Pacientes_K1_DF,Grupo_K5_Luad_Pacientes_K1.2_DF,by="Sample_ID")
colnames(Grupo5_Luad_Comprobacion) <- c("K3.1","ID","K5.1.2")


#########################################################################################################

#GRUPO 2: en K=3

Grupo_K3_Luad_K2 <- K3_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==2)
Grupo_K3_Luad_Pacientes_K2 <-Grupo_K3_Luad_K2$Sample_ID

Grupo_K3_Luad_Pacientes_K2 <- Grupo_K3_Luad_Pacientes_K2[order(Grupo_K3_Luad_Pacientes_K2)]
Grupo_K3_Luad_Pacientes_K2_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K3_Luad_Pacientes_K2),1),
                                            Sample_ID=Grupo_K3_Luad_Pacientes_K2)
Grupo_K3_Luad_Pacientes_K2_DF

#GRUPO 3 y 4: en K=5;

Grupo_K5_Luad_K3.4 <- K5_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==3|nmf_subtypes==4)
Grupo_K5_Luad_Pacientes_K3.4 <-Grupo_K5_Luad_K3.4$Sample_ID

Grupo_K5_Luad_Pacientes_K3.4 <- Grupo_K5_Luad_Pacientes_K3.4[order(Grupo_K5_Luad_Pacientes_K3.4)]
Grupo_K5_Luad_Pacientes_K3.4_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K5_Luad_Pacientes_K3.4),1),
                                              Sample_ID=Grupo_K5_Luad_Pacientes_K3.4)

Grupo5_Luad <- K5_Grupos_Luad_Ordenados[,c(1,2)]
Pacientes_K3.2_K5.3.4 <- inner_join(Grupo_K3_Luad_Pacientes_K2_DF,K5_Grupos_Luad_Ordenados,by="Sample_ID")

#Paciente que estan en los clusteres 3 y 4 de K=5: 
PacientesUnicosdeK5.3.4 <- setdiff(Grupo_K5_Luad_Pacientes_K3.4,Grupo_K3_Luad_Pacientes_K2) #Estos 5 pacientes no estan en el grupo 2 del K=3

#########################################################################################################

#GRUPO 3 en K=3:

Grupo_K3_Luad_K3 <- K3_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==3)
Grupo_K3_Luad_Pacientes_K3 <-Grupo_K3_Luad_K3$Sample_ID

Grupo_K3_Luad_Pacientes_K3 <- Grupo_K3_Luad_Pacientes_K3[order(Grupo_K3_Luad_Pacientes_K3)]
Grupo_K3_Luad_Pacientes_K3_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K3_Luad_Pacientes_K3),1),
                                            Sample_ID=Grupo_K3_Luad_Pacientes_K3)


#GRUPO 5  en K=5;

Grupo_K5_Luad_K5 <- K5_Grupos_Luad_Ordenados %>% filter(nmf_subtypes==5)
Grupo_K5_Luad_Pacientes_K5 <-Grupo_K5_Luad_K5$Sample_ID

Grupo_K5_Luad_Pacientes_K5 <- Grupo_K5_Luad_Pacientes_K5[order(Grupo_K5_Luad_Pacientes_K5)]
Grupo_K5_Luad_Pacientes_K5_DF <- data.frame(NumeroPaciente=seq(1,length(Grupo_K5_Luad_Pacientes_K5),1),
                                              Sample_ID=Grupo_K5_Luad_Pacientes_K5)

#Los 5 pacientes que estaban en, K=3 en el anterior punto estan en K=5?

PacientesUnicosdeK5.3.4 <- setdiff(Grupo_K5_Luad_Pacientes_K3.4,Grupo_K3_Luad_Pacientes_K2)
PacientesUnicosdeK5.3.4_DF <- data.frame(NumeroPaciente=seq(1,length(PacientesUnicosdeK5.3.4),1),
                                         Sample_ID=PacientesUnicosdeK5.3.4)

PacientesFaltan <- inner_join(Grupo_K3_Luad_Pacientes_K3_DF,PacientesUnicosdeK5.3.4_DF,by="Sample_ID")





