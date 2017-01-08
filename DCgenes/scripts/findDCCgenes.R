# findDCCgenes.R
# 2016-12-22
# produces lists of genes that are either dosage compensated or not, according to several criteria
setwd("/media/jenny/670FC52648FA85C4/Documents/MeisterLab/dSMF/DCgenes/scripts")

library(data.table)
# using suplemental table 3 from Kramer ... Ercan et al. PlotGen (2015) 
# "Developmental Dynamics of DCC and H4K20me1 in C. elegans"
# this table contains log2 fold changes and adjusted p values from DEseq2
# most interested in 
# L3 stage hermaphrodites vs mixed (A: no change, X: upregulated herm)
# or L3 stage dpy-27 mutant vs control (A: no change, X: upregulated mutant)

#load DC gene list from Julie's RNAseq data
JulieDCC<-read.table("/media/jenny/670FC52648FA85C4/Documents/MeisterLab/otherPeopleProjects/Julie_RNAseq/DCCgenes_JulieRNAseq.xls",
                     stringsAsFactors=FALSE)

expn<-read.csv("../dataSets/Kramer2015/Kramer2015_S3_File_DEseq2results.csv",
               header=TRUE, sep=",",stringsAsFactors=FALSE)

#remove columns for (dpy21), set1 and set4 deletions
#remove<-grep("dpy21",names(expn))
#expn<-expn[,-remove]
remove<-grep("set1",names(expn))
expn<-expn[,-remove]
remove<-grep("set4",names(expn))
expn<-expn[,-remove]

####### get list of which genes on X vs autosomes
library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

annotFile<-"/SharedDocuments/MeisterLab/GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3.gz"
genomeVer="WS250"

txdb<-makeTxDbFromGFF(annotFile,format="auto",dataSource="WormBase",organism="Caenorhabditis elegans")
genesWS250<-genes(txdb)

XchrGenes<-genesWS250[seqnames(genesWS250)=="X",]
AchrGenes<-genesWS250[seqnames(genesWS250) %in% c("I","II","III","IV","V"),]
############3


i<-which(expn$Gene_WB_ID %in% XchrGenes$gene_id)
expn_X<-expn[i,]
expn_A<-expn[-i,]

## published DCC genes:
pubDCgenes<-c("WBGene00003514",# myo-2
              "WBGene00000149", # uvt-4 / apl-1
              "WBGene00003003", # lin-14
              "WBGene00023497", # lin-15B
              "WBGene00009621", # fipr-21
              "WBGene00009813") # haly-1
    
pubDCgenes %in% expn_X$Gene_WB_ID
# all found in expn_X

pubNDCgenes<-c("WBGene00006316", # sup-7 (tRNA)
               "WBGene00006327", # sup-21 (tRNA)
               "WBGene00015791") # C15C7.5
#not found in expn_X (not protein coding)
pubNDCgenes %in% expn_X$Gene_WB_ID
#last one
pubNDCgenes %in% expn_A$Gene_WB_ID
#none

L3_hermVmixd_A<-expn_A[expn_A[,"hermaphrodite_mixed_sex_L3_log2_padj"]>0.05 & 
                         abs(expn_A[,"hermaphrodite_mixed_sex_L3_log2_fold_change"])<log2(2),]
L3_hermVmixd_A<-L3_hermVmixd_A[complete.cases(L3_hermVmixd_A[,c("hermaphrodite_mixed_sex_L3_log2_padj",
                                              "hermaphrodite_mixed_sex_L3_log2_fold_change")]),]
#14433 genes

## further subset autosomal genes to ones that do or don't chagne upon depletion of dpy-27
L3_A_dc<-L3_hermVmixd_A[complete.cases(L3_hermVmixd_A[,"dpy27_RNAi_L3_padj"]),]
L3_A_dc<-L3_A_dc[L3_A_dc[,"dpy27_RNAi_L3_padj"]<=0.05,]
#134 genes affected by DC
L3_A_ndc<-L3_hermVmixd_A[complete.cases(L3_hermVmixd_A[,c("dpy27_RNAi_L3_padj","dpy27_RNAi_L3_log2_fold_change")]),]
L3_A_ndc<-L3_A_ndc[L3_A_ndc[,"dpy27_RNAi_L3_padj"]>0.05 &
                     abs(L3_A_ndc[,"dpy27_RNAi_L3_log2_fold_change"])<log2(1.4),]
#12731 autosomal genes not affected by DC deletion

saveRDS(L3_A_dc,file="L3_A_dc.Rds")
saveRDS(L3_A_ndc,file="L3_A_ndc.Rds")

## subset to genes that do not change significantly between hermaphrodites and mixed sex population
L3_hermVmixd_X<-expn_X[expn_X[,"hermaphrodite_mixed_sex_L3_log2_padj"]>0.05 & 
                         abs(expn_X[,"hermaphrodite_mixed_sex_L3_log2_fold_change"])<log2(2),]
L3_hermVmixd_X<-L3_hermVmixd_X[complete.cases(L3_hermVmixd_X[,c("hermaphrodite_mixed_sex_L3_log2_padj",
                                                                "hermaphrodite_mixed_sex_L3_log2_fold_change")]),]
#2354 genes
pubDCgenes %in% L3_hermVmixd_X$Gene_WB_ID
#myo-2 removed
pubNDCgenes %in% L3_hermVmixd_X$Gene_WB_ID
# prtn genes there

sum(is.na(match(JulieDCC[,1],expn_X$Gene_WB_ID)))
# 3 genes not expressed (let-4,ncRNAx2)  #WBGene00002282

#checking for Julie's genes
sum(is.na(match(as.vector(JulieDCC[,1]),L3_hermVmixd_X$Gene_WB_ID)))
#43 genes excluded => 129 genes from Julie's list equally expressed

## further subset to genes on X that change significantly upon deletion of dpy-27
L3_X_dc<-L3_hermVmixd_X[complete.cases(L3_hermVmixd_X[,"dpy27_RNAi_L3_padj"]),]
L3_X_dc<-L3_X_dc[L3_X_dc[,"dpy27_RNAi_L3_padj"]<=0.05,]
# 162 genes

pubDCgenes %in% L3_X_dc$Gene_WB_ID
# 5th remains


genesIncluded<-L3_X_dc$Gene_WB_ID
# get a subset for genes from Julie that are DC
sum(is.na(match(as.vector(JulieDCC[,1]),L3_X_dc$Gene_WB_ID)))
#142 genes excluded => could add another 99 genes to dc list
#take genes that are not in final list but ARE amoung genes equally expressed in hermaphrodites
# and mixed populations
takeNew<-JulieDCC[!(JulieDCC[,1] %in% L3_X_dc$Gene_WB_ID) & 
  (JulieDCC[,1] %in% L3_hermVmixd_X$Gene_WB_ID),1]

L3_X_dc_Julie<-L3_hermVmixd_X[L3_hermVmixd_X[,"Gene_WB_ID"] %in% takeNew,]
#99 genes
saveRDS(L3_X_dc_Julie,file="L3_X_dc_Julie.Rds")



saveRDS(L3_X_dc,file="L3_X_dc.Rds")
#saveRDS(L3_X_dc_fc,file="L3_X_dc_fc.Rds") # too few "survive" to add this filter


## also genes with same expression at adult stage
Ad_hermVmixd_X<-expn_X[complete.cases(expn_X[,c("hermaphrodite_mixed_sex_young_adult_padj",
                                                "hermaphrodite_mixed_sex_young_adult_log2_fold_change")]),]

Ad_hermVmixd_X<-Ad_hermVmixd_X[Ad_hermVmixd_X[,"hermaphrodite_mixed_sex_young_adult_padj"]>0.05 & 
                         abs(Ad_hermVmixd_X[,"hermaphrodite_mixed_sex_young_adult_log2_fold_change"])<log2(2),]
#1703 genes
pubDCgenes %in% Ad_hermVmixd_X$Gene_WB_ID
# all but 4th


########################## additional genes selected here below #################
#already included L3_X_dc in genesIncluded variable above
## update list of genes already chosen for primer design with Julie data (above):
genesIncluded<-union(genesIncluded,L3_X_dc_Julie$Gene_WB_ID)

#### add further subsets of genes. Choosing first geens that are equally expressed in  
# wt herm and mixed population, and then subset those that are dc at any stage

L3same_X_dcAny<-L3_hermVmixd_X[complete.cases(L3_hermVmixd_X[,c("dpy27_RNAi_mixed_embryos_padj",
                                                                "dpy27_mutant_L1_padj")]),]
L3same_X_dcAny<-L3same_X_dcAny[L3same_X_dcAny[,"dpy27_RNAi_mixed_embryos_padj"]<=0.05 |
                                 L3same_X_dcAny[,"dpy27_mutant_L1_padj"]<=0.05,]
# 42 genes
pubDCgenes %in% L3same_X_dcAny$Gene_WB_ID
# only last two

#remove genes already in L3 list
i<-which(L3same_X_dcAny$Gene_WB_ID %in% genesIncluded)
L3same_X_dcAny<-L3same_X_dcAny[-i,]
#28 genes

#update genesIncluded
genesIncluded<-union(genesIncluded,L3same_X_dcAny$Gene_WB_ID)

saveRDS(L3same_X_dcAny,file="L3same_X_dcAny.Rds")



#### also add genes same in adults and dc at any stage
Adsame_X_dcAny<-Ad_hermVmixd_X[complete.cases(Ad_hermVmixd_X[,c("dpy27_RNAi_L3_padj",
                                                                "dpy27_RNAi_mixed_embryos_padj",
                                                                "dpy27_mutant_L1_padj")]),]
#this deletes too many cases. but anything else is more complicated to code.
Adsame_X_dcAny<-Adsame_X_dcAny[c(Adsame_X_dcAny[,"dpy27_RNAi_L3_padj"]<=0.05 |
                                Adsame_X_dcAny[,"dpy27_RNAi_mixed_embryos_padj"]<=0.05 |
                                Adsame_X_dcAny[,"dpy27_mutant_L1_padj"]<=0.05),]
#100 genes
pubDCgenes %in% Adsame_X_dcAny$Gene_WB_ID
#last two remain

#remove genes already in L3 list
i<-which(Adsame_X_dcAny$Gene_WB_ID %in% genesIncluded)
Adsame_X_dcAny<-Adsame_X_dcAny[-i,]
# 8 genes
saveRDS(Adsame_X_dcAny,file="Adsame_X_dcAny.Rds")

#update genesIncluded
genesIncluded<-union(genesIncluded,Adsame_X_dcAny$Gene_WB_ID)


### export list of L3same for which dpy21 evidence of dc exists
L3_X_dpy21<-L3_hermVmixd_X[complete.cases(L3_hermVmixd_X[,c("dpy21_mutant_L3_padj",
                                          "dpy21_mutant_L3_log2_fold_change")]),]
L3_X_dpy21<-L3_X_dpy21[(L3_X_dpy21[,"dpy21_mutant_L3_padj"]<=0.05 &
                           L3_X_dpy21[,"dpy21_mutant_L3_log2_fold_change"]>log2(1.5)),]
# >1.5 fold change: 587 genes 
pubDCgenes %in% L3_X_dpy21$Gene_WB_ID
#2,3,5

#remove genes already in included list
i<-which(L3_X_dpy21$Gene_WB_ID %in% genesIncluded)
L3_X_dpy21<-L3_X_dpy21[-i,]
#446 genes

saveRDS(L3_X_dpy21,file="L3_X_dpy21.Rds")

#update genesIncluded
genesIncluded<-union(genesIncluded,L3_X_dpy21$Gene_WB_ID)
pubDCgenes %in% genesIncluded
#2,3,5,6


### export list of genes for which there was no evidence that expn is same
# in males and hermaphrodites, but at least there is DC evidence from dpy-27 deletion.
L3notSame_X_dcAny<-expn_X[complete.cases(expn_X[,c("dpy27_RNAi_mixed_embryos_padj",
                                                   "dpy27_mutant_L1_padj",
                                                   "dpy27_RNAi_L3_padj")]),]
L3notSame_X_dcAny<-L3notSame_X_dcAny[L3notSame_X_dcAny[,"dpy27_RNAi_mixed_embryos_padj"]<=0.05 |
                                       L3notSame_X_dcAny[,"dpy27_mutant_L1_padj"]<=0.05 |
                                       L3notSame_X_dcAny[,"dpy27_RNAi_L3_padj"]<=0.05,]
#201 genes
pubDCgenes %in% L3notSame_X_dcAny$Gene_WB_ID
#last two remain

i<-which(L3notSame_X_dcAny$Gene_WB_ID %in% genesIncluded)
L3notSame_X_dcAny<-L3notSame_X_dcAny[-i,]
#11 genes

saveRDS(L3notSame_X_dcAny,file="L3notSame_X_dcAny.Rds")

#update genesIncluded
genesIncluded<-union(genesIncluded,L3notSame_X_dcAny$Gene_WB_ID)
pubDCgenes %in% genesIncluded


### export list of remaining JulieDCC genes for which there was no evidence that expn is same
# in males and hermaphrodites, but at least there is DCC evidence.
i<-setdiff(JulieDCC[,1], genesIncluded) 

L3notSame_X_dpy26tev<-expn_X[expn_X$Gene_WB_ID %in% i,]
#38 genes
saveRDS(L3notSame_X_dpy26tev,file="L3notSame_X_dpy26tev.Rds")

#update genesIncluded
genesIncluded<-union(genesIncluded,L3notSame_X_dpy26tev$Gene_WB_ID)
pubDCgenes %in% genesIncluded
#1,2,3,5,6

### finally export all other chr X genes not in previous lists
i<-which(L3_hermVmixd_X$Gene_WB_ID %in% genesIncluded)
L3_X_ndc<-L3_hermVmixd_X[-i,]
#1619 genes
pubDCgenes %in% L3_X_ndc$Gene_WB_ID
#4

saveRDS(L3_X_ndc,file="L3_X_ndc.Rds")

#update genesIncluded
genesIncluded<-union(genesIncluded,L3_X_ndc$Gene_WB_ID)

