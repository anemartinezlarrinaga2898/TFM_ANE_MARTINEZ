# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-06-2021

#OBJETIVO:ANALISIS BP 
#########################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
#########################################################################################################

# ······································L U A D ····················································

# ·················C O M P A R A C I O N  1

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Luad.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF1_Luad.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF1U_Luad.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF1D_Luad.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")



GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Luad_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Luad_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Luad_Filtrado_DOWN.rds")

# ·················C O M P A R A C I O N  2

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Luad.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF2_Luad.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF2U_Luad.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF2D_Luad.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Luad_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Luad_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Luad_Filtrado_DOWN.rds")

GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

# ·················C O M P A R A C I O N  3

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Luad.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF3_Luad.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF3U_Luad.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF3D_Luad.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Luad_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Luad_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Luad_Filtrado_DOWN.rds")

GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
pdf("prubea.pdf")
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

#########################################################################################################

# ······································L U S C ····················································

# ·················C O M P A R A C I O N  1

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Lusc.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF1_Lusc.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF1U_Lusc.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF1D_Lusc.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Lusc_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Lusc_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo1_Lusc_Filtrado_DOWN.rds")

GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

# ·················C O M P A R A C I O N  2

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Lusc.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF2_Lusc.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF2U_Lusc.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF2D_Lusc.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Lusc_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Lusc_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo2_Lusc_Filtrado_DOWN.rds")


GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

# ·················C O M P A R A C I O N  3

GT <- readRDS("Analisis_Filtro3/Objetos_ORA/TOTALES/GFT_Lusc.rds")
GF <- readRDS("Analisis_Filtro3/Objetos_ORA/FILTRADOS/GF3_Lusc.rds")
GU <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF3U_Lusc.rds")
GD <- readRDS("Analisis_Filtro3/Objetos_ORA/UP_DOWN/GF3D_Lusc.rds")

EF <- enrichGO(gene=GF,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

EU <- enrichGO(gene=GU,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

ED <- enrichGO(gene=GD,
               universe = GT,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               pAdjustMethod = "BH",
               ont="BP")

saveRDS(EF,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Lusc_Filtrado_Filtrado.rds")
saveRDS(EU,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Lusc_Filtrado_UP.rds")
saveRDS(ED,"Analisis_Filtro3/Resultados/ORA/ORA_Grupo3_Lusc_Filtrado_DOWN.rds")


GF <- barplot(EF, showCategory = 10, title="TOTAL")
GF
GU <- barplot(EU,showCategory=10, title = "UP")
GU
GD <- barplot(ED,showCategory=10, title="DOWN")
GD

