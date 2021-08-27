# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 19-07-2021

# OBJETIVO:Encontrar los factores de transcripcion que se encuentran DE
########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(VennDiagram)
directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/17_FactoresTranscripcion/")
################################################################################################################################

#································ L U A D ···················································

FTLuad <- readRDS("FT_TOTALES_LUAD.rds")

H1 <- FTLuad %>% filter(padj<0.05) %>% filter(Grupo=="Grupo1")
H1 <- H1[!duplicated(H1$Symbol),]
rownames(H1) <- H1$Symbol
H2 <- FTLuad %>% filter(padj<0.05) %>% filter(Grupo=="Grupo2")
H2 <- H2[!duplicated(H2$Symbol),]
rownames(H2) <- H2$Symbol
H3 <- FTLuad %>% filter(padj<0.05) %>% filter(Grupo=="Grupo3")
H3 <- H3[!duplicated(H3$Symbol),]
rownames(H3) <- H3$Symbol

NH1 <- H1$Symbol
NH2 <- H2$Symbol
NH3 <- H3$Symbol

venn.diagram(
  x = list(NH1, NH2,NH3),
  category.names = c("Group 1" , "Group 2","Group 3"),
  filename = "Huellas_Luad.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#A7DA8F", '#F4DC55',"#55CBF4"),
  fill = c(alpha("#A7DA8F",0.5), alpha('#F4DC55',0.5),alpha('#55CBF4',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

OverLap_TG1_TG2 <- calculate.overlap(x=list("Grupo1"=NH1,
                                            "Grupo 2"=NH2,
                                            "Grupo 3"=NH3))

N1_Unicos <- OverLap_TG1_TG2[[5]]
N2_Unicos <- OverLap_TG1_TG2[[6]]
N3_Unicos <- OverLap_TG1_TG2[[7]]

N1_Unicos <- H1[N1_Unicos,]
N2_Unicos <- H1[N2_Unicos,]
N3_Unicos <- H1[N3_Unicos,]

H1_H3 <- rbind(H1,H3)
H1_H3 <- H1_H3[order(H1_H3$Symbol),]

saveRDS(N1_Unicos,"Grupo1_LUAD_FT_Unicos.rds")
saveRDS(N2_Unicos,"Grupo2_LUAD_FT_Unicos.rds")
saveRDS(N3_Unicos,"Grupo3_LUAD_FT_Unicos.rds")

################################################################################################################################

#································ L U S C ···················································

FTLusc <- readRDS("FT_TOTALES_LUSC.rds")

H1 <- FTLusc %>% filter(padj<0.05) %>% filter(Grupo=="Grupo1")
H1 <- H1[!duplicated(H1$Symbol),]
rownames(H1) <- H1$Symbol
H2 <- FTLusc %>% filter(padj<0.05) %>% filter(Grupo=="Grupo2")
H2 <- H2[!duplicated(H2$Symbol),]
rownames(H2) <- H2$Symbol
H3 <- FTLusc %>% filter(padj<0.05) %>% filter(Grupo=="Grupo3")
H3 <- H3[!duplicated(H3$Symbol),]
rownames(H3) <- H3$Symbol

NH1 <- H1$Symbol
NH2 <- H2$Symbol
NH3 <- H3$Symbol

venn.diagram(
  x = list(NH1, NH2,NH3),
  category.names = c("Group 1" , "Group 2","Group 3"),
  filename = "Huellas_Lusc.png",
  output=TRUE,
  lwd = 2,
  lty = 1,
  col=c("#6988E1", '#F3679C',"#B78FDA"),
  fill = c(alpha("#6988E1",0.5), alpha('#F3679C',0.5),alpha('#B78FDA',0.5)),
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
)

OverLap_TG1_TG2 <- calculate.overlap(x=list("Grupo1"=NH1,
                                            "Grupo 2"=NH2,
                                            "Grupo 3"=NH3))

N1_Unicos <- OverLap_TG1_TG2[[5]]
N2_Unicos <- OverLap_TG1_TG2[[6]]
N3_Unicos <- OverLap_TG1_TG2[[7]]

N1_Unicos <- H1[N1_Unicos,]
N2_Unicos <- H1[N2_Unicos,]
N3_Unicos <- H1[N3_Unicos,]

saveRDS(N1_Unicos,"Grupo1_LUSC_FT_Unicos.rds")
saveRDS(N2_Unicos,"Grupo2_LUSC_FT_Unicos.rds")
saveRDS(N3_Unicos,"Grupo3_LUSC_FT_Unicos.rds")

H1_H3 <- rbind(H2,H3)
H1_H3 <- H1_H3[order(H1_H3$Symbol),]


