#Comprobacion de huellos: 

#Vamos a comprar cuantos elementos son unicos de cada huella con el filtro2: 

direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/12_Huellas/")

# luad 


Huella1 <- readRDS("Huella_Ensembel_Luad_1.rds")
Huella2 <- readRDS("Huella_Ensembel_Luad_2.rds")
Huella3 <- readRDS("Huella_Ensembel_Luad_3.rds")

ElementosUnicosHuella1_Huella2 <- setdiff(Huella1,Huella2)
# 74 son unicos 
# 38 estan repetidos 

ElementosUnicosHuella1_Huella3 <- setdiff(Huella1,Huella3)
# 93 son unicos 
# 18 estan repetidos 

ElementosUnicosHuella2_Huella1 <- setdiff(Huella2,Huella1)
# 231 son unicos 
# 38 estan repetidos 

ElementosUnicosHuella2_Huella3 <- setdiff(Huella2,Huella3)
# 202 son unicos 
# 67 repetidos 

ElementosUnicosHuella3_Huella2 <- setdiff(Huella3,Huella2)
# 34 son unicos 
# 67 repetidos 

ElementosUnicosHuella3_Huella1 <- setdiff(Huella3,Huella1)
#82 son unicos 
#19 estan repetidos. 