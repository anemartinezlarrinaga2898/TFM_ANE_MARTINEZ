# TFM, Master en Metodos Computacionales UNAV
#Autor: Ane Martinez Larrinaga 
#Tutores: Silvestre Vicent Cambra y Mikel Hernaez
#Fecha: 02-07-2021

#OBJETIVO:

###########################################################################################################################################################################################################################################################################################v###################################
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

directory <- setwd("/Users/anemartinezlarrinaga/OneDrive/2_MASTER_Computacional/5_TFM/CODIGO/SCRIPTS_Y_DATOS/16_InfiltracionInmunologica/")

################################################################################################################################

#SISTEMA ADPATATIVO: 

CelulasB <- c("MS4A1,TCL1A ,HLA-DOB, PNOC,KIAA0125,CD19,CR2,FCRL2,BLK,IGHG1,COCH,OSBPL10,IGHA1,TNFRSF17,ABCB4,BLNK,GLDC,MEF2C,IGHM,FAM30A,SPIB,BCL11A,GNG7,IGKC,CD72,MICAL3,BCL11A,BACH2,IGL@,CCR9,QRSL1,DTNB,HLA-DQA1,SCN3A,QRSL1,SLC15A2")
CelulasB <- as.vector(str_split(CelulasB,pattern = ",",simplify = TRUE))

CelulasT <- c("PRKCQ,CD3D,CD3G,CD28,LCK,TRAT1,PRKCQ,BCL11B,CD2,LCK,TRBC1,ITM2A,SH2D1A,CD6,CD96,NCALD,GIMAP5,CD3E,SKAP1,TRA@")
CelulasT <- as.vector(str_split(CelulasT,pattern = ",",simplify = TRUE))

CelulasTFoliculares <- c("CHI3L2, CXCL13, MYO74, CHGB, ICAI, HEY1, CDK5R1, ST8SIA1, PDCD1, BLR1, KIAA1234, PVALB, TSHR, C18orf1, HEY1, TOX,BLR1, SLC7A10, SMAD1, POMT1, PASK, MKL2, PTPN13, PASK, KCNK5, ZNF764, MAF, MYO6,SIRPG, THADA, MAGEH1, B3GAT1, SH3TC1, HISTI1H3K, STK39")
CelulasTFoliculares <- as.vector(str_split(CelulasTFoliculares,pattern = ",",simplify = TRUE))

CelulasTReguladoras <-"FOXP3"

CD8 <- c("CD8B,CD8A,CD8B,PF4,PRR5,SF1,LIME1,DNAJB1,ARHGAP8,GZMM,SLC16A7,SFRS7,APBA2,C4orf15,LEPROTL1,ZFP36L2,GADD45A,ZFP36L2,MYST3,ZEB1,ZNF609,C12orf47,THUMPD1,VAMP2,ZNF91,ZNF22,TMC6,DNAJB1,FLT3LG,CDKN2AIP,TSC22D3,TBCC,RBM3,ABT1,C19orf6,CAMLG,PPP1R2,AES,KLF9,PRF1,TRD@")
CD8 <- as.vector(str_split(CD8,pattern = ",",simplify = TRUE))

CD4 <- c("KLRD1,KLRF1,GNLY,CTSW,KLRB1,KLRK1,NKG7,GZMH,KLRD1,SIGIRR,ZBTB16,RUNX3,APOL3,RORA,APBA2,SIGIRR,WHDC1L1,DUSP2,GZMA")
CD4 <- as.vector(str_split(CD4,pattern = ",",simplify = TRUE))

NK <- c("LOC643313,GAGE2,ZNF747,XCL1,XCL2,AF107846,SLC30A5,MCM3AP,TBXA2R,CDC5L,LOC7300096,FUT5,FGF18,MRC2,SPN,PSMD4,PRX,FZR1,ZNF205,ZNF528,BCL2,FZR1,PDLMI4,LDB3,ADARB1,SMKE1,TCTN2,TINAGL1,IGFBP5,ALDH1B1,NCR1")
NK <- as.vector(str_split(NK,pattern = ",",simplify = TRUE))

DC <- c("CD209,CCL17,HSD11B1,CCL13,CCL22,PPFIBP2,NPR1")
DC <- as.vector(str_split(DC,pattern = ",",simplify = TRUE))

aDC <- c("CCL1,EBI3,INDO,LAMP3,OAS3")
aDC <- as.vector(str_split(aDC,pattern = ",",simplify = TRUE))

Eosinofilos <- c("IL5RA,KCNH2,TKTL1,IL5RA,EMR1,KCNH2,CCR3,ACACB,THBS1,GALC,TKTL1,RNU2CLC,THBS1,HIST1H1C,CYSLTR2,HRH4,RNASE2,CAT, LRP5L,SYNJ1,THBS4GPR44,KBTBD11,HES1,ABHD2,TIPARP,SMPD3,MYO15B,TGIF1,RRP12,ACACBIGSF2,HES1,RCOR3,EPN2,C9orf156,SIAH1,ACACB")
Eosinofilos <- as.vector(str_split(Eosinofilos,pattern = ",",simplify = TRUE))

Macrofagos <- c("MARCO,CXCL5,SCG5, ,SULT1C2,MSR1,CTSK,PTGDS,COLEC12,GPC4,MSR1,PCOLCE2,CHIT1,PTGDS,KAL1,CLEC5A,GPC4,ME1,DNASE2B,CCL7,FN1,CD163,GM2A,SCARB2,BCAT1,RAI14,MSR1,COL8A2,CD163,APOE,CHI3L1,ATG7,CD84,FDX1,MS4A4A,SGMS1,EMP1,CYBB,CD68")
Macrofagos <- as.vector(str_split(Macrofagos,pattern = ",",simplify = TRUE))

Neutrofilos <- c("CSF3R,CYP4F3,VNN3,FPRL1,KCNJ15,MME,IL8RA,IL8RB,MME,FCGR3B,DYSF,KCNJ15,FCAR ,CEACAM3,FPRL,HIST1H2BC,HPSE,FLJ11151,CREB5,S100A12,FCGR3B,TNFRSF10C,SLC22A4,KIAA0329,SLC25A37,BST1,FCAR,CEACAM3,CRISPLD2,TNFRSF10C,G0S2,SIGLEC5,CD93,MGAM,ALPL,FPR1,CD93,PDE4B,LILRB2")
Neutrofilos <- as.vector(str_split(Macrofagos,pattern = ",",simplify = TRUE))

PDL1 <- "CD27A"
CTLA4 <- "CTLA4"
PDCD1 <- "PDCD1"
PDCD2 <- "PDCD1LG2"

HuellasSistemaInmune <- list(CelulasB,CelulasT,CelulasTFoliculares,CelulasTReguladoras,
                             CD8,CD4,NK,DC,aDC,Eosinofilos,Macrofagos,Neutrofilos,
                             PDL1,CTLA4,PDCD1,PDCD2)
names(HuellasSistemaInmune) <- c("B Cells","T Cells ","T Folicular","T Regulator","CD8","CD4","Natural Killer","Dendritic Cell","Active Dendritic Cells","Eosinophils","Macrophages","Neutrophils",
                                 "PDL1","CTLA4","PDCD1","PDCD2")
saveRDS(HuellasSistemaInmune,"HuellasSistemaInmune.rds")

################################################################################################################################

HuellasSistemaInmune <- readRDS("HuellasSistemaInmune.rds")

HuellasDF <- data.frame()

i <- 1
for(i in seq_along(HuellasSistemaInmune)){
  Celula <- HuellasSistemaInmune[[i]]
  
  Celula_DF <- data.frame(GenesInmune=Celula)
  Celula_DF$TipoCelular <- names(HuellasSistemaInmune)[i]
  
  HuellasCelulasDF <- rbind(HuellasDF,Celula_DF)
  
  HuellasDF <- HuellasCelulasDF
  
  
}

saveRDS(HuellasCelulasDF,"HuellasSistemaInmune_DF.rds")



