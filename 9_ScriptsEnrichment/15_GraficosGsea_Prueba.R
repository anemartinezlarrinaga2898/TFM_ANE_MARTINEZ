# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 18-05-2021

#OBJETIVO:Graficos de gseaPLOT de los 4 primeros tanto de los filtrados como del total


##############################################################################################################
direcotry <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/9_ScriptsEnrichment/")
library(tidyverse)
library(ggplot2) 
library(enrichplot)
library(ggpubr)
##############################################################################################################

# ···································· L U A D················································

# C O M P A R A C I O N   1: ··························································································

RFiltrados <- readRDS("AnalisisGSEA_BP/GenesFiltrados/AnalisisGSEA_Luad_1_Filtrados.rds")
RTotales<- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_1_Total.rds")

GeneSet_1 <- enrichplot::gseaplot2(RFiltrados,pvalue_table = TRUE ,geneSetID = "GO:0009966",title="Regulation of signal transduction")
GeneSet_2<- enrichplot::gseaplot2(RFiltrados,pvalue_table = TRUE,geneSetID = "GO:0071704",title="Organic substance metabolic process")
GeneSet_3 <- enrichplot::gseaplot2(RFiltrados,pvalue_table=TRUE,geneSetID = "GO:0044238",title="Primary metabolic process")
GeneSet_4 <- enrichplot::gseaplot2(RFiltrados,pvalue_table = TRUE,geneSetID = "GO:0009967",title="Positive regulation of signal transduction")

GeneSet_1
GeneSet_2
GeneSet_3
GeneSet_4

GeneSet_1 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE ,geneSetID = "GO:0002682",title="Regulation of immune system process")
GeneSet_2<- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0046649",title="Lymphocyte activation")
GeneSet_3 <- enrichplot::gseaplot2(RTotales,pvalue_table=TRUE,geneSetID = "GO:0045321",title="Positive regulation of immune system process")
GeneSet_4 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0006952",title="Leukocyte activation")

GeneSet_1
GeneSet_2
GeneSet_3
GeneSet_4

# C O M P A R A C I O N   2: ··························································································

RTotales<- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_2_Total.rds")

GeneSet_1 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE ,geneSetID = "GO:0006952",title="Defense Response")
GeneSet_2<- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0045321",title="Leukocyte activation")
GeneSet_3 <- enrichplot::gseaplot2(RTotales,pvalue_table=TRUE,geneSetID = "GO:0002684",title="Positive regulation of immune system process")
GeneSet_4 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0006954",title="Inflammatory response")

GeneSet_1
GeneSet_2
GeneSet_3
GeneSet_4

# C O M P A R A C I O N   3: ··························································································

RTotales<- readRDS("AnalisisGSEA_BP/GenesTotales/AnalisisGSEA_Luad_3_Total.rds")

GeneSet_1 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE ,geneSetID = "GO:0006952",title="Defense Response")
GeneSet_2<- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0045321",title="Leukocyte activation")
GeneSet_3 <- enrichplot::gseaplot2(RTotales,pvalue_table=TRUE,geneSetID = "GO:0002684",title="Positive regulation of immune system process")
GeneSet_4 <- enrichplot::gseaplot2(RTotales,pvalue_table = TRUE,geneSetID = "GO:0006954",title="Inflammatory response")

GeneSet_1
GeneSet_2
GeneSet_3
GeneSet_4
