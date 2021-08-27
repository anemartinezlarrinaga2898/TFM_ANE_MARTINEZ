# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 14-06-2021

#OBJETIVO:Nuevo script para obtener los heatmaps. En este caso es necesario que las lineas celulas esten en las columnas y los genes en filas. 


################################################################################################################################################

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/11_ScriptCancerCellLine/")
library(tidyverse)
library(pheatmap)
#Cargamos los datos de mis huellas: 
Huellas <- readRDS("Huellas_LineasCelulares_Lusc_SinNan.rds")
Condiciones <- c("Grupo 1","Grupo 2","Grupo 3")


#Vamos ha hacer un bucle para añadir una columna en la que ponga cada gen de que grupo es. 
#Vamos a anotar las huellas, y vamos a crear un DF con las anotaciones para poder ponerlo en el heatmap. 
DataFramesAnotaciones <- list()
i <- 1
for ( i in seq_along(Huellas)){
  #Seleccionamos la huella correspondiente
  Huella <- data.frame(Huellas[[i]])
  #A la huella seleccionada le añadimos la condicion con la que vamos a trabajar
  HuellaNueva <- Huella %>% mutate(Grupo=Condiciones[i])
  
  #Generamos un data.frame con las anotaciones
  DataFramesAnotaciones[[i]] <- data.frame(Genes=rownames(HuellaNueva),Grupo=HuellaNueva$Grupo)
}

#Unimos las cuentas en un unica matriz; en el orden, GRUPO 1, GRUPO 2 , GRUPO 3. 
MatrizUnica <- rbind(Huellas[[1]],Huellas[[2]])
MatrizUnica <- rbind(MatrizUnica,Huellas[[3]])

saveRDS(MatrizUnica,"Matriz_Huellas_Lineas.rds")
saveRDS(DataFramesAnotaciones,"Anotaciones_Heatmap.rds")
#Las anotaciones estaban por separado , cada una en un elemento de una lista, lo que hacemos es coger, 
#cada elemento y los unimos uno debajo del otro. 

AnotacionesUnica <- rbind(DataFramesAnotaciones[[1]],DataFramesAnotaciones[[2]])
AnotacionesUnica <- rbind(AnotacionesUnica,DataFramesAnotaciones[[3]])

#Como no podemos añadirle el nombre de los genes porque las huellas se repite ( lo cual es interesante )
#Añadimo una secuencia de numero y cada gen de cada grupo tendra un numero concreto .

SecuenciaCodigo <- seq(1,nrow(MatrizUnica),1)
AnotacionesUnica <- cbind(AnotacionesUnica,SecuenciaCodigo)
AnotacionesUnica <- AnotacionesUnica[,c(3,2,1)] #En este data.frame es donde tengo cada numero a que gen corresponde. 
rownames(AnotacionesUnica) <- AnotacionesUnica[,1]

#Para comprobar que el orden es el mismo
rownames(MatrizUnica)==AnotacionesUnica$Genes
rownames(MatrizUnica) <- AnotacionesUnica[,1]

saveRDS(AnotacionesUnica,"Anotaciones_Heatmap.rds")
#Matriz unica ahora mismo esta: 
  # Lineas celulares en columnas 
  # Genes en filas. 

#Por lo tanto el heatmap: 
  # 1) Las anotaciones van en las filas 
  # 2) Quiero los nombres de las columnas. 


H <- MatrizUnica
Colors <- c("green","black","red")
breaks <- seq(-2,2,by=0.1)
Pallete1 <- colorRampPalette(Colors)(length(breaks))
pheatmap(H, 
         color = Pallete1, 
         cluster_rows =F , 
         clustering_distance_cols = "euclidean",
         clustering_method="ward.D",
         breaks = breaks,
         cluster_cols = T, 
         show_colnames = T,
         show_rownames = F,
         annotation_row=AnotacionesUnica[2])

plot(hclust(dist(t(H))))

############################################################################################################

#Voy a intentar ha hacerlo con otra funcion: 

library/

