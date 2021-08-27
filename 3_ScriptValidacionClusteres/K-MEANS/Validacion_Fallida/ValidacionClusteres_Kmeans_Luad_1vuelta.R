muestras <- readRDS('MuestrasLuadTraspuesta.rds')

library(NbClust)

#CALCULO DE UN POOL DE INDICES:
  
  #Calculamos cada uno de los indices: 
system.time(res_nb_silhouette <- NbClust(muestras, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='silhouette'))
print("Calculado Silhouette_luad")

system.time(res_nb_dindex     <- NbClust(muestras, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='dindex'))
print("Calculado Dindex_luad")

system.time(res_nb_dunn       <- NbClust(muestras, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='dunn'))
print("Calculado Dunn_luad")

system.time(res_nb_gap        <- NbClust(muestras, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='gap'))
print("Calculado gap_luad")

system.time(res_nb_tau        <- NbClust(muestras, distance = 'euclidean', method = 'kmeans', min.nc=2, max.nc=10, index='tau'))
print("Calculado tau_luad")

  
saveRDS(res_nb_silhouette, 'res_nb_silhouette_luad.rds') 
saveRDS(res_nb_dindex, 'res_nb_dindex_luad.rds')
saveRDS(res_nb_dunn, 'res_nb_dunn_luad.rds')
saveRDS(res_nb_gap, 'res_nb_gap_luad.rds')
saveRDS(res_nb_tau, 'res_nb_tau_luad.rds')

print("Calculos terminados luad")