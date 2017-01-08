#2016-11-16
# convert supplemental tables from Kreus(2013) and Chen(2013). Column titles were changed in the tables
# to remove spaces. Removed empty columns in Kreus data. Used convertingGeneNames.R on Chen data to add 
# WBGeneXXXXXX gene names to the excel table. 
# tables were then exported as .csv for import into this script.
# saved WS235 files in bed gff and Rds formats (bed has no metadata, gff converts numbers to characters)

setwd("/SharedDocuments/MeisterLab/dSMF/TSS/scripts/")
library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
source("./getTSS_functions.R")

######### create TSS file from WB annotation - do once #######################
annotFile<-"/SharedDocuments/MeisterLab/GenomeVer/annotations/c_elegans.PRJNA13758.WS240.annotations.gff3.gz"
genomeVer="WS240"
library(GenomicFeatures)
txdb<-makeTxDbFromGFF(annotFile,format="auto",dataSource="WormBase",organism="Caenorhabditis elegans")
promReg <- promoters(txdb, upstream=0,downstream=1,columns=c("gene_id","tx_id"))
#txpts<-transcripts(txdb)
gnId <- sapply(mcols(promReg)$gene_id, paste, collapse=",")
promRegSel <- promReg[ match(unique(gnId), gnId) ]
names(promRegSel) <- gsub("Gene:","",unique(gnId))
export(promRegSel,paste0("TSS_1perGene_",genomeVer,".bed"),format="bed")
export(promRegSel,paste0("TSS_1perGene_",genomeVer,".gff3"),format="gff3")
export(promReg,paste0("TSS_1perTxpt_",genomeVer,".bed"),format="bed")
export(promReg,paste0("TSS_1perTxpt_",genomeVer,".gff3"),format="gff3")


########## create TSS file from Kreus et al. 2013 ########################3
KreusTSSfile<-"../dataSets/Kreus2013/TSS_Kreus_2013_elife_suplfig1-data2.csv"
#Note Kreus data uses WS230 (ce10)
library(genomation)
library(GenomicRanges)
KreusTSS<-readGeneric(KreusTSSfile,chr=1,start=5,end=6,
            keep.all.metadata=TRUE,header=TRUE,sep=",")
strand(KreusTSS)<-ifelse(mcols(KreusTSS)$Strand==1,"+","-")
#seqlevels(KreusTSS)<- sub('chr','',seqlevels(KreusTSS))
library(rtracklayer)
# use liftover in R directly on GRanges object:
# ce10toce11<-import.chain("/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain")
# KreusTSS_WS235<-liftOver(KreusTSS,ce10toce11)
# KreusTSS_WS235<-unlist(KreusTSS_WS235)
# this gives a granges list with disjoint granges for some
# genes. When then unlisted to create a granges object, you 
# get more than one range per gene. Command line liftover does not do that

# do liftover using a system call to the command line and re-import
export(KreusTSS,"../dataSets/Kreus2013/KreusTSS_WS230.bed",format="bed")
liftOverProg<-"/SharedDocuments/MeisterLab/GenomeVer/liftOver" # where your lift over program is installed
chainFile<-"/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain" #appropriate chain file
dataDir<-"../dataSets/Kreus2013/" # folder where your datasets are
data2lift<-paste0(dataDir,"KreusTSS_WS230.bed") # bed file to liftover
liftedData<-paste0(dataDir,"KreusTSS_WS235.bed") # name of output bedfile
# systerm call to liftover program: 
system(paste(liftOverProg,data2lift,chainFile,liftedData, paste0(dataDir,"unlifted.bed"))) #file for unlifted seq (just in case it doesn't work)

# import new Granges object with new coordinates
KreusTSS_WS235<-import("../dataSets/Kreus2013/KreusTSS_WS235.bed",format="bed")
# attach metadata to new genome coordinates (NOTE: metadata is still using old coordinates!!!)
mcols(KreusTSS_WS235)<-mcols(KreusTSS[,-3])

# list columns that contain coordinate information that needs to be updated 
names(mcols(KreusTSS_WS235))
# make list of column names with coordinate information
cols2Convert<-c("DCCmutantEmbryo_TSS","WildTypeEmbryo_TSS","WildTypeStarvedL1_TSS",
                "WildTypeL3_TSS","Allen2011_TSAPosMostSL1reads_WS230")
# lift over coordinates for metadata columns 
KreusTSS_WS235<-liftOverMcols(KreusTSS_WS235,cols2Convert)
export(KreusTSS_WS235,"../dataSets/Kreus2013/KreusTSS_WS235.gff3",format="gff",version="3")
saveRDS(KreusTSS_WS235,"../dataSets/Kreus2013/KreusTSS_WS235.Rds")


####################### do same of gene expression files ###############3
### embryos wt and DCC
KreusGeneExprFile<-"../dataSets/Kreus2013/GeneExprEmbDCC_Kreus_2013_elife_suplfig4-data1.csv"
KreusGeneExpr_embDCC<-readGeneric(KreusGeneExprFile,chr=1,start=5,end=6,
                      keep.all.metadata=TRUE,header=TRUE,sep=",")
strand(KreusGeneExpr_embDCC)<-ifelse(mcols(KreusGeneExpr_embDCC)$Strand==1,"+","-")

export(KreusGeneExpr_embDCC,"../dataSets/Kreus2013/KreusGeneExpr_embDCC_WS230.bed",format="bed")
liftOverProg<-"/SharedDocuments/MeisterLab/GenomeVer/liftOver" # where your lift over program is installed
chainFile<-"/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain" #appropriate chain file
dataDir<-"../dataSets/Kreus2013/" # folder where your datasets are
data2lift<-paste0(dataDir,"KreusGeneExpr_embDCC_WS230.bed") # bed file to liftover
liftedData<-paste0(dataDir,"KreusGeneExpr_embDCC_WS235.bed") # name of output bedfile
# systerm call to liftover program: 
system(paste(liftOverProg,data2lift,chainFile,liftedData, paste0(dataDir,"unlifted.bed"))) #file for unlifted seq (just in case it doesn't work)

# import new Granges object with new coordinates
KreusGeneExpr_embDCC_WS235<-import("../dataSets/Kreus2013/KreusGeneExpr_embDCC_WS235.bed",format="bed")
# attach metadata to new genome coordinates (NOTE: metadata is still using old coordinates!!!)
mcols(KreusGeneExpr_embDCC_WS235)<-mcols(KreusGeneExpr_embDCC[,-3])

export(KreusGeneExpr_embDCC_WS235,"../dataSets/Kreus2013/KreusGeneExpr_embDCC_WS235.gff3",format="gff3")
saveRDS(KreusGeneExpr_embDCC_WS235,"../dataSets/Kreus2013/KreusGeneExpr_embDCC_WS235.Rds")


### starved L1s
KreusGeneExprFile<-"../dataSets/Kreus2013/GeneExprStarvedL1s_Kreus_2013_elife_suplfig4-data1.csv"
KreusGeneExpr_starvedL1s<-readGeneric(KreusGeneExprFile,chr=1,start=5,end=6,
                                  keep.all.metadata=TRUE,header=TRUE,sep=",")
strand(KreusGeneExpr_starvedL1s)<-ifelse(mcols(KreusGeneExpr_starvedL1s)$Strand==1,"+","-")

export(KreusGeneExpr_starvedL1s,"../dataSets/Kreus2013/KreusGeneExpr_starvedL1s_WS230.bed",format="bed")
liftOverProg<-"/SharedDocuments/MeisterLab/GenomeVer/liftOver" # where your lift over program is installed
chainFile<-"/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain" #appropriate chain file
dataDir<-"../dataSets/Kreus2013/" # folder where your datasets are
data2lift<-paste0(dataDir,"KreusGeneExpr_starvedL1s_WS230.bed") # bed file to liftover
liftedData<-paste0(dataDir,"KreusGeneExpr_starvedL1s_WS235.bed") # name of output bedfile
# systerm call to liftover program: 
system(paste(liftOverProg,data2lift,chainFile,liftedData, paste0(dataDir,"unlifted.bed"))) #file for unlifted seq (just in case it doesn't work)

# import new Granges object with new coordinates
KreusGeneExpr_starvedL1s_WS235<-import("../dataSets/Kreus2013/KreusGeneExpr_starvedL1s_WS235.bed",format="bed")
# attach metadata to new genome coordinates (NOTE: metadata is still using old coordinates!!!)
mcols(KreusGeneExpr_starvedL1s_WS235)<-mcols(KreusGeneExpr_starvedL1s[,-3])

export(KreusGeneExpr_starvedL1s_WS235,"../dataSets/Kreus2013/KreusGeneExpr_starvedL1s_WS235.gff3",format="gff3")
saveRDS(KreusGeneExpr_starvedL1s_WS235,"../dataSets/Kreus2013/KreusGeneExpr_starvedL1s_WS235.Rds")


### L3s
KreusGeneExprFile<-"../dataSets/Kreus2013/GeneExprL3s_Kreus_2013_elife_suplfig4-data1.csv"
KreusGeneExpr_L3s<-readGeneric(KreusGeneExprFile,chr=1,start=5,end=6,
                                      keep.all.metadata=TRUE,header=TRUE,sep=",")
strand(KreusGeneExpr_L3s)<-ifelse(mcols(KreusGeneExpr_L3s)$Strand==1,"+","-")

export(KreusGeneExpr_L3s,"../dataSets/Kreus2013/KreusGeneExpr_L3s_WS230.bed",format="bed")
liftOverProg<-"/SharedDocuments/MeisterLab/GenomeVer/liftOver" # where your lift over program is installed
chainFile<-"/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain" #appropriate chain file
dataDir<-"../dataSets/Kreus2013/" # folder where your datasets are
data2lift<-paste0(dataDir,"KreusGeneExpr_L3s_WS230.bed") # bed file to liftover
liftedData<-paste0(dataDir,"KreusGeneExpr_L3s_WS235.bed") # name of output bedfile
# systerm call to liftover program: 
system(paste(liftOverProg,data2lift,chainFile,liftedData, paste0(dataDir,"unlifted.bed"))) #file for unlifted seq (just in case it doesn't work)

# import new Granges object with new coordinates
KreusGeneExpr_L3s_WS235<-import("../dataSets/Kreus2013/KreusGeneExpr_L3s_WS235.bed",format="bed")
# attach metadata to new genome coordinates (NOTE: metadata is still using old coordinates!!!)
mcols(KreusGeneExpr_L3s_WS235)<-mcols(KreusGeneExpr_L3s[,-3])

export(KreusGeneExpr_L3s_WS235,"../dataSets/Kreus2013/KreusGeneExpr_L3s_WS235.gff3",format="gff3")
saveRDS(KreusGeneExpr_L3s_WS235,"../dataSets/Kreus2013/KreusGeneExpr_L3s_WS235.Rds")


############# import Chen 2013 data ###############
ChenTSSfile<-"../dataSets/Chen2013/TSS_Chen2013_Supp_TableS2.csv"
#Note Chen data uses WS220 (ce10)
library(genomation)
library(GenomicRanges)
ChenTSS<-readGeneric(ChenTSSfile,chr=3,start=4,end=5,strand=6,
                      keep.all.metadata=TRUE,header=TRUE,sep=",")

seqlevels(ChenTSS)<- paste0("chr",seqlevels(ChenTSS))

hist(width(ChenTSS),breaks=100)
quantile(width(ChenTSS),c(0.90,0.95,0.99))
#90% 95% 99% 
#86 114 182
sum(width(ChenTSS)<=100)/length(width(ChenTSS))
#0.9293469


library(rtracklayer)

# note: one region did not liftover:
# #Partially deleted in new
# chrIII	3162620	3162660	.	0	+

GR2del<-GRanges(seqnames=c("chrIII"),IRanges(start=3162620,end=3162660),strand="+")
findOverlaps(GR2del,ChenTSS)
ChenTSS<-ChenTSS[-35991,]

# do liftover using a system call to the command line and re-import
export(ChenTSS,"../dataSets/Chen2013/ChenTSS_WS230.bed",format="bed")
liftOverProg<-"/SharedDocuments/MeisterLab/GenomeVer/liftOver" # where your lift over program is installed
chainFile<-"/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain" #appropriate chain file
dataDir<-"../dataSets/Chen2013/" # folder where your datasets are
data2lift<-paste0(dataDir,"ChenTSS_WS230.bed") # bed file to liftover
liftedData<-paste0(dataDir,"ChenTSS_WS235.bed") # name of output bedfile
# systerm call to liftover program: 
system(paste(liftOverProg,data2lift,chainFile,liftedData, paste0(dataDir,"unlifted.bed"))) #file for unlifted seq (just in case it doesn't work)

# import new Granges object with new coordinates
ChenTSS_WS235<-import("../dataSets/Chen2013/ChenTSS_WS235.bed",format="bed")
# attach metadata to new genome coordinates (NOTE: metadata is still using old coordinates!!!)
mcols(ChenTSS_WS235)<-mcols(ChenTSS)

# list columns that contain coordinate information that needs to be updated 
names(mcols(ChenTSS_WS235))
# make list of column names with coordinate information
cols2Convert<-c("modePosition")
# lift over coordinates for metadata columns 
ChenTSS_WS235<-liftOverMcols(ChenTSS_WS235,cols2Convert)

export(ChenTSS_WS235,"../dataSets/Chen2013/ChenTSS_WS235.gff3",format="gff",version="3")
saveRDS(ChenTSS_WS235,"../dataSets/Chen2013/ChenTSS_WS235.Rds")
