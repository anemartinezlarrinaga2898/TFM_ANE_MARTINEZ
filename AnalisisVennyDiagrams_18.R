directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/15_AnalisisEnrichment/")

GenesUnicosLuad <- readRDS("Comparacion18/Anotaciones_NivelesExpresion_TG1.rds")
GenesUnicosLuad <- GenesUnicosLuad[,-1]
GenesUnicosLusc <- readRDS("Comparacion18/Anotaciones_NivelesExpresion_TG2.rds")
GenesUnicosLusc <- GenesUnicosLusc[,-1]

Huellaslusc <- readRDS("Huellas_Lusc_Anotadas.rds")
Huellaslusc <- Huellaslusc[,-1]
Huellaslusc <- Huellaslusc %>% filter(ValoresExpresion< (-2) | ValoresExpresion > (2))

Huellasluad <- readRDS("HuellasLuad_Anotadas.rds")
Huellasluad <- Huellasluad[,-2]
Huellasluad <- Huellasluad %>% filter(Valores< (-2) | Valores > (2))

# ¿A que grupo pertenece cada gen unico?

GenUnico_Grupo_Luad <- inner_join(GenesUnicosLuad,Huellasluad,by="hgnc_symbol")
GenUnico_Grupo_Lusc <- inner_join(GenesUnicosLusc,Huellaslusc,by="hgnc_symbol")

G1 <- GenUnico_Grupo_Luad %>% filter(Grupo=="Grupo1")
G1 <- G1[!duplicated(G1$hgnc_symbol),]
G2 <- GenUnico_Grupo_Luad %>% filter(Grupo=="Grupo2")
G2 <- G2[!duplicated(G2$hgnc_symbol),]
G3 <- GenUnico_Grupo_Luad %>% filter(Grupo=="Grupo3")
G3 <- G3[!duplicated(G3$hgnc_symbol),]

G1 <- GenUnico_Grupo_Lusc %>% filter(Grupo=="Grupo1")
G1 <- G1[!duplicated(G1$hgnc_symbol),]
G2 <- GenUnico_Grupo_Lusc %>% filter(Grupo=="Grupo2")
G2 <- G2[!duplicated(G2$hgnc_symbol),]
G3 <- GenUnico_Grupo_Lusc %>% filter(Grupo=="Grupo3")
G3 <- G3[!duplicated(G3$hgnc_symbol),]

##############################################################################################################

#¿Cuantos de la comparacion 18 coinciden con la comparacion 17, dentro de lusc?

Lusc18 <- readRDS("Comparacion18/Anotaciones_NivelesExpresion_TG1.rds")
G1 <- readRDS("Comparacion17/G1.rds")
G2 <- readRDS("Comparacion17/G2.rds")
G3 <- readRDS("Comparacion17/G3.rds")

G1_18_17 <- inner_join(G1,GenUnico_Grupo_Lusc,by="hgnc_symbol")
G1_18_17 <- G1_18_17[,-c(5,6,8)]
G1_18_17 <- G1_18_17[!duplicated(G1_18_17$ensembl_gene_id),]

G2_18_17 <- inner_join(G2,GenUnico_Grupo_Lusc,by="hgnc_symbol")
G2_18_17 <- G2_18_17[,-c(5,6,8)]
G2_18_17 <- G2_18_17[!duplicated(G2_18_17$ensembl_gene_id),]

G3_18_17 <- inner_join(G3,GenUnico_Grupo_Lusc,by="hgnc_symbol")
G3_18_17 <- G3_18_17[,-c(5,6,8)]
G3_18_17 <- G3_18_17[!duplicated(G3_18_17$ensembl_gene_id),]

##############################################################################################################

#¿Cuantos de la comparacion 18 coinciden con la comparacion 16, dentro de luad?

G1 <- readRDS("Comparacion16/G1.rds")
G2 <- readRDS("Comparacion16/G2.rds")
G3 <- readRDS("Comparacion16/G3.rds")

G1_18_17 <- inner_join(G1,GenUnico_Grupo_Luad,by="hgnc_symbol")
G1_18_17 <- G1_18_17[,-c(5,6,8)]
G1_18_17 <- G1_18_17[!duplicated(G1_18_17$ensembl_gene_id),]

G2_18_17 <- inner_join(G2,GenUnico_Grupo_Luad,by="hgnc_symbol")
G2_18_17 <- G2_18_17[,-c(5,6,8)]
G2_18_17 <- G2_18_17[!duplicated(G2_18_17$ensembl_gene_id),]

G3_18_17 <- inner_join(G3,GenUnico_Grupo_Luad,by="hgnc_symbol")
G3_18_17 <- G3_18_17[,-c(5,6,8)]
G3_18_17 <- G3_18_17[!duplicated(G3_18_17$ensembl_gene_id),]

