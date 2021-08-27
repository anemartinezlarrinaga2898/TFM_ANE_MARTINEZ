set.seed(1234)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/AnalisisLusc/")
muestras <- readRDS('MuestrasLusc.rds')
muestras <- t(muestras)
dimensiones <- dim(muestras)

#####################################################################################################################################

print("Trabajo con todas las muestras")

library(ClusterR)

NumeroCluster_inicial <- 2
NumeroCluster_final <- 10
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)

IndiceSilhouette_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                                             max_clusters = clusteresTotales, #Numero de clusteres
                                                             criterion ="silhouette", #Criterio de la loss function
                                                             num_init=100, #Numero de veces que se va a inicializar
                                                             initializer="kmeans++", #Metodo de inicializacion
                                                             plot_clusters=TRUE) #Si luego quiero que me plotee el resultado

saveRDS(IndiceSilhouette_TodasLasMuestras,"IndiceSilhouette_TodasLasMuestras_Lusc.rds")     
print("IndiceSilhouette_TodasLasMuestras CALCULADO y GUARDADO LUSC")

IndiceVarianza_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                                           max_clusters = clusteresTotales, #Numero de clusteres
                                                           criterion ="variance_explained", #Criterio de la loss function
                                                           num_init=100, #Numero de veces que se va a inicializar
                                                           initializer="kmeans++", #Metodo de inicializacion
                                                           plot_clusters=TRUE) #Si luego quiero que me plotee el resultado

saveRDS(IndiceVarianza_TodasLasMuestras,"IndiceVarianza_TodasLasMuestras_Lusc.rds")     
print("IndiceVarianza_TodasLasMuestras CALCULADO y GUARDADO LUSC")

IndiceDissimilarity_TodasLasMuestras <- Optimal_Clusters_KMeans(data=muestras,
                                                                max_clusters = clusteresTotales, #Numero de clusteres
                                                                criterion ="dissimilarity", #Criterio de la loss function
                                                                num_init=100, #Numero de veces que se va a inicializar
                                                                initializer="kmeans++", #Metodo de inicializacion
                                                                plot_clusters=TRUE) #Si luego quiero que me plotee el resultado
saveRDS(IndiceDissimilarity_TodasLasMuestras,"IndiceDissimilarity_TodasLasMuestras_Lusc.rds")     
print("IndiceDissimilarity_TodasLasMuestras CALCULADO y GUARDADO LUSC")

resultadosKmean_TodasLasMuestras <- list()
i <- 1
for (i in seq_along(clusteresTotales)){
  print(i)
  k <- clusteresTotales[i]
  kmeans.resultado <- kmeans(x=muestras,centers = k,nstart=100)
  resultadosKmean[[i]] <- kmeans.resultado
}

saveRDS(resultadosKmean_TodasLasMuestras,"ResutladosKmean_TodasLasMuestras_Lusc.rds")
print("Calculo de TODAS las MUESTRAS REALIZADOS LUSC")

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

saveRDS(IndiceSilhouette_ParteDeLasMuestras,"IndiceSilhouette_ParteDeLasMuestras_Lusc.rds")
saveRDS(IndiceVarianza_ParteDeLasMuestras,"IndiceVarianza_ParteDeLasMuestras_Lusc.rds")
saveRDS(IndiceDissimilarity_ParteDeLasMuestras,"IndiceDissimilarity_ParteDeLasMuestras_Lusc.rds")

print("Indices para parte de las muestras calculados LUSC")

resultadosKmeans_ParteDeLasMuestras <- list()
MatricesConsenso <- list()

i <- 1

for (i in seq_along(clusteresTotales)){
  k <- clusteresTotales[i]
  MatricesDeConectividad <- list()
  for (j in seq(1,100,1)){
    subsample <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ]
    kmeans.output <- kmeans(subsample,nstart=100,centers = k)
    resultadosKmeans_ParteDeLasMuestras[[j]] <- kmeans.output
    conectividad <- connectivity(kmeans.output)
    MatricesDeConectividad[[j]] <- conectividad
    
  }
  consenso <- consensus(MatricesDeConectividad)
  MatricesConsenso[[i]] <- consenso
}

saveRDS(MatricesConsenso,"MatricesConsenso_Lusc.rds")
saveRDS(resultadosKmeans_ParteDeLasMuestras,"resultadosKmeans_ParteDelasMuestras_Lusc.rds")

print("Calculos hechos y guardados, LUSC")


