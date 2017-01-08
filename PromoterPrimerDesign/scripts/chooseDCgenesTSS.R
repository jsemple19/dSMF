#2016-12-23
#chooseDCgenesTSS.R
# find genes that are both in list of dosage compensated genes and have a well defined TSS shared
# between Kruesi and Chen datasets
setwd("/media/jenny/670FC52648FA85C4/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")

#read in lists of DC genes
L3_X_dc<-readRDS("../../DCgenes/scripts/L3_X_dc.Rds")
L3_X_dc_Julie<-readRDS("../../DCgenes/scripts/L3_X_dc_Julie.Rds")
L3same_X_dcAny<-readRDS("../../DCgenes/scripts/L3same_X_dcAny.Rds")
Adsame_X_dcAny<-readRDS("../../DCgenes/scripts/Adsame_X_dcAny.Rds")
L3_X_dpy21<-readRDS("../../DCgenes/scripts/L3_X_dpy21.Rds")
L3notSame_X_dcAny<-readRDS("../../DCgenes/scripts/L3notSame_X_dcAny.Rds")
L3notSame_X_dpy26tev<-readRDS("../../DCgenes/scripts/L3notSame_X_dpy26tev.Rds")
L3_X_ndc<-readRDS("../../DCgenes/scripts/L3_X_ndc.Rds")

#read in lists of genes with defined TSS
TSS<-readRDS("../../TSS/scripts/confirmedTSS_Chen&Kreus_WS235.Rds")
names(TSS)<-mcols(TSS)$WBGeneID

i<-which(L3_X_dc$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 42 genes
L3_X_dc_TSS<-L3_X_dc[i,]
L3_X_dc_badTSS<-L3_X_dc[-i,]
saveRDS(L3_X_dc_TSS,file="L3_X_dc_TSS.Rds")
#saveRDS(L3_X_dc_badTSS,file="L3_X_dc_badTSS.Rds")


i<-which(L3_X_dc_Julie$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 23 genes
L3_X_dc_Julie_TSS<-L3_X_dc_Julie[i,]
saveRDS(L3_X_dc_Julie_TSS,file="L3_X_dc_Julie_TSS.Rds")


i<-which(L3same_X_dcAny$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 7 genes
L3same_X_dcAny_TSS<-L3same_X_dcAny[i,]
saveRDS(L3same_X_dcAny_TSS,file="L3same_X_dcAny_TSS.Rds")

i<-which(Adsame_X_dcAny$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 0 genes
Adsame_X_dcAny_TSS<-Adsame_X_dcAny[i,]
saveRDS(Adsame_X_dcAny_TSS,file="Adsame_X_dcAny_TSS.Rds")

i<-which(L3_X_dpy21$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 89 genes for -75 to 25
L3_X_dpy21_TSS<-L3_X_dpy21[i,]
#L3_X_dpy21_badTSS<-L3_X_dpy21[-i,]
saveRDS(L3_X_dpy21_TSS,file="L3_X_dpy21_TSS.Rds")
#saveRDS(L3_X_dpy21_badTSS,file="L3_X_dpy21_badTSS.Rds")

i<-which(L3notSame_X_dcAny$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 2 genes
L3notSame_X_dcAny_TSS<-L3notSame_X_dcAny[i,]
saveRDS(L3notSame_X_dcAny_TSS,file="L3notSame_X_dcAny_TSS.Rds")

i<-which(L3notSame_X_dpy26tev$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 7 genes
L3notSame_X_dpy26tev_TSS<-L3notSame_X_dpy26tev[i,]
saveRDS(L3notSame_X_dpy26tev_TSS,file="L3notSame_X_dpy26tev_TSS.Rds")

i<-which(L3_X_ndc$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 213 genes for -75 to 25
L3_X_ndc_TSS<-L3_X_ndc[i,]
#L3_X_ndc_badTSS<-L3_X_ndc[-i,]
saveRDS(L3_X_ndc_TSS,file="L3_X_ndc_TSS.Rds")
#saveRDS(L3_X_ndc_badTSS,file="L3_X_ndc_badTSS.Rds")


chosen<-list("L3_dc_TSS"=L3_X_dc_TSS$Gene_WB_ID[,drop=TRUE],
     #"L3_dc_badTSS"=L3_X_dc_badTSS$Gene_WB_ID[,drop=TRUE],
     "L3_dc_Julie_TSS"=L3_X_dc_Julie_TSS$Gene_WB_ID[,drop=TRUE],
     "L3_dcAny_TSS"=L3same_X_dcAny_TSS$Gene_WB_ID[,drop=TRUE],
     "Adsame_dcAny_TSS"=Adsame_X_dcAny_TSS$Gene_WB_ID[,drop=TRUE],
     #"L3_dcAny_badTSS"=L3same_X_dcAny_badTSS$Gene_WB_ID[,drop=TRUE],
     "L3_dpy21_TSS"=L3_X_dpy21_TSS$Gene_WB_ID[,drop=TRUE],
     #"L3_dpy21_badTSS"=L3_X_dpy21_badTSS$Gene_WB_ID[,drop=TRUE],
     "L3notSame_dcAny_TSS"=L3notSame_X_dcAny_TSS$Gene_WB_ID[,drop=TRUE],
     "L3notSame_dpy26tev_TSS"=L3notSame_X_dpy26tev_TSS$Gene_WB_ID[,drop=TRUE],
     "L3_ndc_TSS"=L3_X_ndc_TSS$Gene_WB_ID[,drop=TRUE])


saveRDS(chosen,"chosenTSS.Rds")



#### do design for autosomal promoters too
#read in lists of DC genes
L3_A_dc<-readRDS("../../DCgenes/scripts/L3_A_dc.Rds")
# 134 genes
L3_A_ndc<-readRDS("../../DCgenes/scripts/L3_A_ndc.Rds")
# some duplicated lines to remove (!?)
# WBGene00008555, WBGene00020649
L3_A_ndc<-L3_A_ndc[!duplicated(L3_A_ndc),]
# 12725 genes

i<-which(L3_A_dc$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 34 genes
L3_A_dc_TSS<-L3_A_dc[i,]
saveRDS(L3_A_dc_TSS,file="L3_A_dc_TSS.Rds")


i<-which(L3_A_ndc$Gene_WB_ID %in% mcols(TSS)$WBGeneID )
length(i)
# 1974 genes
L3_A_ndc_TSS<-L3_A_ndc[i,]

saveRDS(L3_A_ndc_TSS,file="L3_A_ndc_TSS.Rds")


chosenA<-list("L3_dc_A_TSS"=L3_A_dc_TSS$Gene_WB_ID[,drop=TRUE],
             "L3_ndc_A_TSS"=L3_A_ndc_TSS$Gene_WB_ID[,drop=TRUE])


saveRDS(chosenA,"chosenTSS_A.Rds")
