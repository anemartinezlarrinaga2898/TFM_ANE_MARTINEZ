# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO: El objetivo de este script es mediante diferentes estrategias llegar a la validacion final del numero de clusteres: 

#####################################################################################################################################

#Antes de empezar a trabajar es necesario: 
#Establecemos un seed para que los resultados puedan ser reproducibles: 
library(ClusterR)
library(NMF)
set.seed(123456789)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/3_ScriptValidacionClusteres/")
muestras <- readRDS('Luad_Filtrado.rds')
#muestras <- readRDS('Lusc_Filtrado.rds')

#Tenemos que hacer la traspuesta porque kmeans espera a los samples en las filas y a los features en las columnas. 
muestras <- t(muestras)
dimensiones <- dim(muestras)

NumeroCluster_inicial <- 2
NumeroCluster_final <- 7
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)

##############################################################################################################################

IndiceSilhouette_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                            max_clusters = clusteresTotales, #Numero de clusteres
                                            criterion ="silhouette", #Criterio de la loss function
                                            num_init=100, #Numero de veces que se va a inicializar
                                            initializer="kmeans++", #Metodo de inicializacion
                                            plot_clusters=TRUE) #Si luego quiero que me plotee el resultado

saveRDS(IndiceSilhouette_TodasLasMuestras,"IndiceSilhouette_TodasLasMuestras_Luad.rds")     
print("IndiceSilhouette_TodasLasMuestras LUAD CALCULADO y GUARDADO")

IndiceVarianza_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                            max_clusters = clusteresTotales, #Numero de clusteres
                                            criterion ="variance_explained", #Criterio de la loss function
                                            num_init=100, #Numero de veces que se va a inicializar
                                            initializer="kmeans++", #Metodo de inicializacion
                                            plot_clusters=TRUE) #Si luego quiero que me plotee el resultado

saveRDS(IndiceVarianza_TodasLasMuestras,"IndiceVarianza_TodasLasMuestras_Luad.rds")     
print("IndiceVarianza_TodasLasMuestras LUAD CALCULADO y GUARDADO")

IndiceDissimilarity_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                          max_clusters = clusteresTotales, #Numero de clusteres
                                          criterion ="dissimilarity", #Criterio de la loss function
                                          num_init=100, #Numero de veces que se va a inicializar
                                          initializer="kmeans++", #Metodo de inicializacion
                                          plot_clusters=TRUE) #Si luego quiero que me plotee el resultado

saveRDS(IndiceDissimilarity_TodasLasMuestras,"IndiceDissimilarity_TodasLasMuestras_Luad.rds")     
print("IndiceDissimilarity_TodasLasMuestras LUAD CALCULADO y GUARDADO")

resultadosKmean_TodasLasMuestras <- list()

i <- 1
for (i in seq_along(clusteresTotales)){
  print(i)
  k <- clusteresTotales[i]
  kmeans.resultado <- kmeans(x=muestras,centers = k,nstart=100)
  resultadosKmean_TodasLasMuestras[[i]] <- kmeans.resultado
}

saveRDS(resultadosKmean_TodasLasMuestras,"ResutladosKmean_TodasLasMuestras_Luad.rds")
print("Calculo de TODAS las MUESTRAS REALIZADOS LUAD")

#####################################################################################################################################

print("Trabajo con parte de las muestras")

IndiceSilhouette_ParteDeLasMuestras <- list()
IndiceVarianza_ParteDeLasMuestras <- list()
IndiceDissimilarity_ParteDeLasMuestras <- list()

i <- 1
for (i in seq(1,100,1)){
  print(i)
  subsample <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ]
  IndiceSilhouette_ParteDeLasMuestras[i] <- Optimal_Clusters_KMeans(data=subsample,
                                                               max_clusters = clusteresTotales, 
                                                               criterion ="silhouette", 
                                                               num_init=100,
                                                               initializer="kmeans++", 
                                                               plot_clusters=TRUE) 
  IndiceVarianza_ParteDeLasMuestras <- Optimal_Clusters_KMeans(data=subsample,
                                                             max_clusters = clusteresTotales, 
                                                             criterion ="variance_explained", 
                                                             num_init=100, 
                                                             initializer="kmeans++", 
                                                           plot_clusters=TRUE)
  IndiceDissimilarity_ParteDeLasMuestras <- Optimal_Clusters_KMeans(data=subsample,
                                                                  max_clusters = clusteresTotales, #Numero de clusteres
                                                                  criterion ="dissimilarity", #Criterio de la loss function
                                                                  num_init=100, #Numero de veces que se va a inicializar
                                                                  initializer="kmeans++", #Metodo de inicializacion
                                                                  plot_clusters=TRUE)
  
}

saveRDS(IndiceSilhouette_ParteDeLasMuestras,"IndiceSilhouette_ParteDeLasMuestras_Luad.rds")
saveRDS(IndiceVarianza_ParteDeLasMuestras,"IndiceVarianza_ParteDeLasMuestras_Luad.rds")
saveRDS(IndiceDissimilarity_ParteDeLasMuestras,"IndiceDissimilarity_ParteDeLasMuestras_Luad.rds")





