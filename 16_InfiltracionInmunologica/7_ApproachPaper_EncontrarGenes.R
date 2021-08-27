# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-07-2021

#OBJETIVO:Encontrar los genes de las huellas

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
library(data.table)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")
############################################################################################################################################################################################################################################################################################

HuellasSistemaInmune <- readRDS("HuellasSistemaInmune_DF.rds")

Huella1 <- readRDS("Huella_1_LUAD_Anotada.rds")
Huella2 <- readRDS("Huella_2_LUAD_Anotada.rds")
Huella3 <- readRDS("Huella_3_LUAD_Anotada.rds")

Huellas <- list(Huella1,Huella2,Huella3)
names(Huellas) <- c("Huella1","Huella2","Huella3")

i <- 1

HuellasInmune1 <- copy(Huella1)
HuellasInmune1[,] <- 0
HuellasInmune1$TipoCelular <- 0
HuellasInmune1 <- HuellasInmune1[1:189,]
rownames(HuellasInmune1) <- seq(1,nrow(HuellasSistemaInmune),1)

HuellasInmune2 <- copy(Huella2)
HuellasInmune2[,] <- 0
HuellasInmune2$TipoCelular <- 0
HuellasInmune2 <- HuellasInmune3[1:189,]
rownames(HuellasInmune2) <- seq(1,nrow(HuellasSistemaInmune),1)

HuellasInmune3 <- copy(Huella3)
HuellasInmune3[,] <- 0
HuellasInmune3$TipoCelular <- 0
HuellasInmune3 <- HuellasInmune3[1:189,]
rownames(HuellasInmune3) <- seq(1,nrow(HuellasSistemaInmune),1)



for(i in seq_along(Huellas)){
  Huella <- Huellas[[i]]
  
    k <- 1
    
    for(k in seq(1,nrow(HuellasSistemaInmune),1)){
    
      X <- HuellasSistemaInmune$hgnc_symbol[k]
    
      GenHuella_Inmune <- Huella[which(Huella$hgnc_symbol==X),]
      
      if(nrow(GenHuella_Inmune)==0){
        if(i==1){
          HuellasInmune1[k,1:3] <- 0
          HuellasInmune1[k,4] <- 0
        }
        
        if(i==2){
          HuellasInmune2[k,1:3] <- 0
          HuellasInmune2[k,4] <-0
        }
        
        if(i==3){
          HuellasInmune3[k,1:3] <- 0
          HuellasInmune3[k,4] <- 0
        }
      }
      
      if(nrow(GenHuella_Inmune)!=0){
        if(i==1){
          HuellasInmune1[k,1:3] <- GenHuella_Inmune
          HuellasInmune1[k,4] <- HuellasSistemaInmune$TipoCelular[k]
        }
        
        if(i==2){
          HuellasInmune2[k,1:3] <- GenHuella_Inmune
          HuellasInmune2[k,4] <- HuellasSistemaInmune$TipoCelular[k]
        }
        
        if(i==3){
          HuellasInmune3[k,1:3] <- GenHuella_Inmune
          HuellasInmune3[k,4] <- HuellasSistemaInmune$TipoCelular[k]
        }
        
      
        
      }
      
    }
    
}

#COMPROBACION CON INNER JOIN: 
colnames(HuellasSistemaInmune)[1] <- colnames(Huella1)[3]
Compracion1 <- inner_join(Huella1,HuellasSistemaInmune,by="hgnc_symbol")
Compracion2 <- inner_join(Huella2,HuellasSistemaInmune,by="hgnc_symbol")
Compracion3 <- inner_join(Huella3,HuellasSistemaInmune,by="hgnc_symbol")


#########################################################################################################################################################

HuellasSistemaInmune <- readRDS("HuellasSistemaInmune_DF.rds")

Huella1 <- readRDS("Huella_1_LUSC_Anotada.rds")
Huella2 <- readRDS("Huella_2_LUSC_Anotada.rds")
Huella3 <- readRDS("Huella_3_LUSC_Anotada.rds")

Huellas <- list(Huella1,Huella2,Huella3)
names(Huellas) <- c("Huella1","Huella2","Huella3")

i <- 1

HuellasInmune1 <- copy(Huella1)
HuellasInmune1[,] <- 0
HuellasInmune1$TipoCelular <- 0
HuellasInmune1 <- HuellasInmune3[1:313,]
rownames(HuellasInmune1) <- seq(1,nrow(HuellasSistemaInmune),1)

HuellasInmune2 <- copy(Huella2)
HuellasInmune2[,] <- 0
HuellasInmune2$TipoCelular <- 0
HuellasInmune2 <- HuellasInmune3[1:313,]
rownames(HuellasInmune2) <- seq(1,nrow(HuellasSistemaInmune),1)

HuellasInmune3 <- copy(Huella3)
HuellasInmune3[,] <- 0
HuellasInmune3$TipoCelular <- 0
HuellasInmune3 <- HuellasInmune3[1:313,]
rownames(HuellasInmune3) <- seq(1,nrow(HuellasSistemaInmune),1)



for(i in seq_along(Huellas)){
  Huella <- Huellas[[i]]
  
  k <- 1
  
  for(k in seq(1,nrow(HuellasSistemaInmune),1)){
    
    X <- HuellasSistemaInmune$GenesInmune[k]
    
    GenHuella_Inmune <- Huella[which(Huella$hgnc_symbol==X),]
    
    if(nrow(GenHuella_Inmune)==0){
      if(i==1){
        HuellasInmune1[k,1:3] <- 0
        HuellasInmune1[k,4] <- 0
      }
      
      if(i==2){
        HuellasInmune2[k,1:3] <- 0
        HuellasInmune2[k,4] <-0
      }
      
      if(i==3){
        HuellasInmune3[k,1:3] <- 0
        HuellasInmune3[k,4] <- 0
      }
    }
    
    if(nrow(GenHuella_Inmune)!=0){
      if(i==1){
        HuellasInmune1[k,1:3] <- GenHuella_Inmune
        HuellasInmune1[k,4] <- HuellasSistemaInmune$TipoCelular[k]
      }
      
      if(i==2){
        HuellasInmune2[k,1:3] <- GenHuella_Inmune
        HuellasInmune2[k,4] <- HuellasSistemaInmune$TipoCelular[k]
      }
      
      if(i==3){
        HuellasInmune3[k,1:3] <- GenHuella_Inmune
        HuellasInmune3[k,4] <- HuellasSistemaInmune$TipoCelular[k]
      }
      
      
      
    }
    
  }
  
}

#COMPROBACION CON INNER JOIN: 
colnames(HuellasSistemaInmune)[1] <- colnames(Huella1)[3]
Compracion1 <- inner_join(Huella1,HuellasSistemaInmune,by="hgnc_symbol")
Compracion2 <- inner_join(Huella2,HuellasSistemaInmune,by="hgnc_symbol")
Compracion3 <- inner_join(Huella3,HuellasSistemaInmune,by="hgnc_symbol")



