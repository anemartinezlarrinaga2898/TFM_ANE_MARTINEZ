muestras <- readRDS('MuestrasLuadTraspuesta.rds')

library(NbClust)

#CALCULO DE UN POOL DE INDICES: 

#Primero creamos las listas vacias que necesitamos para actuar
res_nb_silhouette <- list() #Aqui se guardan los resultados que obtener del resultado
res_nb_dindex <- list()
res_nb_dunn <- list()
res_nb_gap <- list()
res_nb_tau <- list()


for(iter in seq(1,10)){
  print(iter)
  #Porque las muestras es la traspuesta entonces tengo en filas pacientes y en columnas los genes. 
  #Cogemos solo el 70% de las muestras y como lo tengo en filas hacemos una indexacion. 
  #Estamos cogiendo todas las columnas y parte de las filas
  subsample <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ]
  
  #Calculamos cada uno de los indices: 
  system.time(res_nb_silhouette_tmp <- NbClust(subsample, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='silhouette'))
  system.time(res_nb_dindex_tmp     <- NbClust(subsample, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='dindex'))
  system.time(res_nb_dunn_tmp       <- NbClust(subsample, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='dunn'))
  system.time(res_nb_gap_tmp        <- NbClust(subsample, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='gap'))
  system.time(res_nb_tau_tmp        <- NbClust(subsample, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='tau'))
  
  #Guardamos cada una de las interacciones en su correspondiente lista
  res_nb_silhouette[iter] <- res_nb_silhouette_tmp
  
  res_nb_dindex[iter]     <- res_nb_dindex_tmp
  
  res_nb_dunn[iter]       <- res_nb_dunn_tmp
  
  res_nb_gap[iter]        <- res_nb_gap_tmp
  
  res_nb_tau[iter]        <- res_nb_tau_tmp
  
}


saveRDS(res_nb_silhouette, 'res_nb_silhouette_10.rds') #AQUI APARECERA UN 100 o un 500 dependiendo de las veces que realicemos el bucle
saveRDS(res_nb_dindex, 'res_nb_dindex_luad_10.rds')
saveRDS(res_nb_dunn, 'res_nb_dunn_luad_10.rds')
saveRDS(res_nb_gap, 'res_nb_gap_luad_10.rds')
saveRDS(res_nb_tau, 'res_nb_tau_luad_10.rds')
