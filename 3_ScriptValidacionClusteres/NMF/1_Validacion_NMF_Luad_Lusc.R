# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 23-03-2021

#OBJETIVO:
#SCRIPT para la validacion de los clusteres generados por NMF: 
#Esto es un script para validar el correlation cophenetic value de mis muestras en ambos grupos. 
#Para validar voy a trabajar con las muestras filtradas. 

set.seed(123456789)
library(NMF)
######################################################################################################

# V A L I D A C I O N  L U A D 

muestras <- readRDS("LuadFiltrado.rds")
muestras <- data.frame(muestras)

NumeroCluster_inicial <- 2
NumeroCluster_final <- 10
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)

#Estimacion realizada con el paquete de NMF: 
# muestras <- as.matrix(muestras)
# EstimationOfRank_Luad <- nmfEstimateRank(muestras,range=clusteresTotales,nrun=100,method="brunet",seed="random")
# saveRDS(EstimationOfRank_Luad,"Estimacion_NMF_Luad_Total.rds")
# 
# #Calculo del cophor value: 
# i <- 1
# valorCorrelacionCofenetci_Total <- vector(mode="logical",length =length(clusteresTotales))
# matricesConsenso <- list()
# Resultados <- list()
# MatrizH <- list()
# MatrizW <- list()
# 
# for (i in seq_along(clusteresTotales)){
#   cluster <- clusteresTotales[i]
#   print(paste("Numero de Cluster",cluster,"en NMF para Luad en Todas las Muestras"))
#   Resultados[[i]] <-  nmf(muestras,cluster,method="brunet",seed="random",nrun=100)
#   valorCorrelacionCofenetci_Total[[i]] <- cophcor(Resultados[[i]])
#   matricesConsenso[[i]] <- consensus(Resultados[[i]])
#   MatrizH[[i]] <- coef(Resultados[[i]])
#   MatrizW[[i]] <- basis(Resultados[[i]])
# }
# 
# saveRDS(valorCorrelacionCofenetci_Total,"ValorCofenetic_Luad_Total.rds")
# saveRDS(matricesConsenso,"MatricesConsenso_Luad_Total.rds")
# saveRDS(MatrizW,"MatrizW_Luad_Total.rds")
# saveRDS(MatrizH,"MatrizH_Luad_Total.rds")

# i <- 1 #Son las filas
# j <- 1 #Son las columnas
# valorCorrelacionCofenetci_ParteMuestras <- matrix(0,nrow=length(clusteresTotales),ncol=100)
# valorCorrelacionCofenetci_ParteMuestras <- data.frame(valorCorrelacionCofenetci_ParteMuestras)
# rownames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("cluster",2:10)
# colnames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("NumeroVuelta",1:ncol(valorCorrelacionCofenetci_ParteMuestras))
# Resultados <- list()
# 
# for (i in seq_along(clusteresTotales)){
#   cluster <- clusteresTotales[i]
#   print(paste("Cluster",cluster,"NMF, Luad con Parte de las muestras"))
#   for (j in seq(1,100,1)){
#     print(paste("NumeroVuelta",j,"Luad, trabajo parte de las muestras"))
#     subsample <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ]
#     Resultados[[j]] <- nmf(subsample,cluster,method="brunet",seed="random",nrun=100)
#     valorCorrelacionCofenetci_ParteMuestras[i,j] <- cophcor(Resultados[[j]]) #Guardamos en el cluster I en la posicion j
#   }
# }
# 
# saveRDS(valorCorrelacionCofenetci_ParteMuestras,"ValorCofenetic_Luad_CrossValidation.rds")
# saveRDS(Resultados,"Resultados_Luad_Total.rds")

print("T R A B A J O  L U A D T E R M I N A D O")

######################################################################################################

# V A L I D A C I O N    L U S C  

muestras <- readRDS("LuscFiltrado.rds")
muestras <- data.frame(muestras)

# muestras <- as.matrix(muestras)
# EstimationOfRank_Lusc <- nmfEstimateRank(muestras,range=clusteresTotales,nrun=1000,method="brunet",seed="random")
# saveRDS(EstimationOfRank_Lusc,"Estimacion_NMF_Lusc_Total.rds")

#Calculo del cophor value: 
i <- 1
valorCorrelacionCofenetci_Total <- vector(mode="logical",length =length(clusteresTotales))
matricesConsenso <- list()
Resultados <- list()
MatrizH <- list()
MatrizW <- list()

for (i in seq_along(clusteresTotales)){
  cluster <- clusteresTotales[i]
  print(paste("Numero de Cluster",cluster,"en NMF para Lusc en Todas las Muestras"))
  Resultados[[i]] <-  nmf(muestras,cluster,method="brunet",seed="random",nrun=100)
  valorCorrelacionCofenetci_Total[[i]] <- cophcor(Resultados[[i]])
  matricesConsenso[[i]] <- consensus(Resultados[[i]])
  MatrizH[[i]] <- coef(Resultados[[i]])
  MatrizW[[i]] <- basis(Resultados[[i]])
}

saveRDS(valorCorrelacionCofenetci_Total,"ValorCofenetic_Lusc_Total_ii.rds")
saveRDS(matricesConsenso,"MatricesConsenso_Lusc_Total_ii.rds")
saveRDS(MatrizW,"MatrizW_Lusc_Total_ii.rds")
saveRDS(MatrizH,"MatrizH_Lusc_Total_ii.rds")

i <- 1 #Son las filas
j <- 1 #Son las columnas
valorCorrelacionCofenetci_ParteMuestras <- matrix(0,nrow=length(clusteresTotales),ncol=100)
valorCorrelacionCofenetci_ParteMuestras <- data.frame(valorCorrelacionCofenetci_ParteMuestras)
rownames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("cluster",2:10)
colnames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("NumeroVuelta",1:ncol(valorCorrelacionCofenetci_ParteMuestras))
Resultados <- list()

for (i in seq_along(clusteresTotales)){
  cluster <- clusteresTotales[i]
  print(paste("Cluster",cluster,"NMF, Lusc con Parte de las muestras"))
  for (j in seq(1,100,1)){
    subsample <- muestras[sample(nrow(muestras), round(nrow(muestras) * 0.7) ), ]
    Resultados[[j]] <- nmf(subsample,cluster,method="brunet",seed="random",nrun=100)
    valorCorrelacionCofenetci_ParteMuestras[i,j] <- cophcor(Resultados[[j]]) #Guardamos en el cluster i en la posicion j
  }
}

saveRDS(valorCorrelacionCofenetci_ParteMuestras,"ValorCofenetic_Lusc_CrossValidation_ii.rds")
saveRDS(Resultados,"Resultados_Lusc_Total.rds")
print("T R A B A J O  L U S C T E R M I N A D O")

######################################################################################################

#CONSENSUS CLUSTERING 

######################################################################################################

#LUAD: 

library(diceR)

muestras <- readRDS("LuadFiltrado.rds")
muestras <- data.frame(muestras)


resultadosConsensus_Luad <- consensus_cluster(data=muestras,
                                              nk = 2:10,
                                              p.item = 0.7,
                                              reps = 1000,
                                              algorithms = "nmf",
                                              nmf.method="brunet",
                                              distance = "euclidean",
                                              scale = FALSE,
                                              seed.data=1262118388.71279,
)

saveRDS(resultadosConsensus_Luad,"ResultadosConsensus_Luad.rds")
print("ConsensusClustering Luad Done")

#LUSC: 

muestras <- readRDS("LuscFiltrado.rds")
muestras <- data.frame(muestras)

resultadosConsensus_Lusc <- consensus_cluster(data=muestras,
                                              nk = 2:10,
                                              p.item = 0.7,
                                              reps = 1000,
                                              algorithms = "nmf",
                                              nmf.method="brunet",
                                              distance = "euclidean",
                                              scale = FALSE,
                                              seed.data=1262118388.71279,
)

saveRDS(resultadosConsensus_Luad,"ResultadosConsensus_Lusc.rds")
print("ConsensusClustering Lusc Done")

#######################################################################################################
