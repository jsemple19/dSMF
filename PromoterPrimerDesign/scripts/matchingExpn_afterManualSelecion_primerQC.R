# 2017-05-08 modified to include better TSS choices
# 2016-12-20
# matchingExpn.R
# takes list of genes for which primers were designed on X chr and assigns them to quantiles of
# L3 gene expression in order to match to autosomal genes

setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")


##############3 after doing manual selection: check for similar expression patterns
finalChosenA<-read.csv("./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                       header=TRUE,stringsAsFactors=FALSE)
finalChosenA<-finalChosenA[finalChosenA$finalSelected=="y","Gene_WB_ID"]

######## recreating "expressed"  and chosenX object
chosenX<-read.csv("./XdcPrimers/DM_DE_promoter_test_primers_WS250_subset.csv",header=TRUE,stringsAsFactors=FALSE)
chosenX<-chosenX[chosenX$finalSelected=="y","Gene_WB_ID"]

#use combined expression filtered by Peter
wormG<-read.csv("../../expressionMatching/Gene_expression_Gerstein_combined.csv",stringsAsFactors=FALSE)
sum(is.na(wormG$Gene))
#convert gene names to wb gene ids
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
conversionTable<-convertGeneNames(wormG$Gene,inputType="seqID")
#add gene name annotation to wormG
i<-match(wormG$Gene,conversionTable$seqID)
wormG<-cbind(conversionTable[i,],wormG)
#remove lines where gene names are NA (often deprcated/merged genes)
i<-which(rowSums(is.na(wormG[,1:3]))>0)
wormG<-wormG[-i,]
#remove genes not expressed in L3
expressed<-wormG[wormG$L3_N2_L3.1>0,]


expnChosen<-rbind(cbind("chr"=rep("chrX",48),expressed[expressed$WBgeneID %in% chosenX,]),
                  cbind("chr"=rep("Autosomes",48),expressed[expressed$WBgeneID %in% finalChosenA,]))

library(ggplot2)
library(reshape2)
#ggplot(expnChosen[,c(1,3,7,8,10:17)],aes(x="chr",y="PolII_EE_FE"))+geom_boxplot()

pdf(file="./primerQC.pdf", width=11,height=8,paper="a4r")

mm=melt(expnChosen[,c(1,3,7,8,10:17)],id=c("WBgeneID","chr"))
names(mm)<-c("WBgeneID","chr","dataset","expression")
ggplot(mm,aes(x=chr,y=log2(expression),fill=chr)) + geom_boxplot() + facet_wrap(~dataset,nrow=2)


############### compare other attirbutes of primers
#read in full chosenX table
chosenX<-read.csv("./XdcPrimers/DM_DE_promoter_test_primers_WS250_subset.csv",header=TRUE,stringsAsFactors=FALSE)
chosenX<-chosenX[chosenX$finalSelected=="y",]

#names(chosenX)<-c(names(chosenX)[1],"fragID",names(chosenX)[2:18])
#names(chosenX)[5]<-"Amplicon"
#library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
# library(genomation)
# library(GenomicRanges)
# library(rtracklayer)
txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% chosenX$Gene_WB_ID,]
i<-match(chosenX$Gene_WB_ID,mcols(genes)$gene_id)
chosenX$strand<-as.vector(strand(genes[i,]))

names(chosenX)[6]<-"orientation"

finalChosenA<-read.csv("./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                       header=TRUE,stringsAsFactors=FALSE)
finalChosenA<-finalChosenA[finalChosenA$finalSelected=="y",]
names(finalChosenA)[6]<-"NameFromBEDfile"
#names(finalChosenA)[2]<-"fragID"
#names(finalChosenA)[5]<-"PrimerID"
names(finalChosenA)[7]<-"Amplicon"
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% finalChosenA$Gene_WB_ID,]
i<-match(finalChosenA$Gene_WB_ID,mcols(genes)$gene_id)
finalChosenA$strand<-as.vector(strand(genes[i,]))



primerData<-rbind(cbind("chr"=rep("chrX",48),"strand"=chosenX$strand,chosenX[,4:19]),
                  cbind("chr"=rep("Autosomes",48),"strand"=finalChosenA$strand,finalChosenA[,c(6,7,5,8:20)]))

primerData<-cbind(primerData,"location"=unlist(lapply(sapply(primerData$Amplicon,strsplit,":"),'[[',1)))
table(primerData$location)
#chrI  chrII chrIII  chrIV   chrV   chrX
#7      9      6     10     16     48

######### WRONG!!! orientation is which strand of DNA is amplified from DNA after BS
table(primerData$orientation,primerData$chr)
#chrX Autosomes
#fw   18        22
#rc   30        26

#names(mm)<-c("wormbase_gene","chr","dataset","expression")
#ggplot(mm,aes(x=chr,y=log2(expression),fill=chr)) + geom_boxplot() + facet_wrap(~dataset,nrow=2)
mm=melt(primerData[,c(1,3,8:18)],id=c("NameFromBEDfile","chr"))
ggplot(mm,aes(x=chr,y=value,fill=chr)) + geom_boxplot() +
    facet_wrap(~variable,scales="free")

mm=melt(primerData[,c(1,3,9:18)],id=c("NameFromBEDfile","chr"))
ggplot(mm,aes(value,fill=chr)) + geom_histogram() +
  facet_wrap(~chr+variable,scales="free",nrow=2)

dev.off()

table(primerData$FwC.covered+primerData$RvC.covered)
#0  1  2  3  4  5  6  7  8  9 10
#6  8 18 12 22 11  7  5  4  1  2

primerTable<-data.frame("primerNames"=c(paste0("X",1:48,"_f"),paste0("A",1:48,"_f"), paste0("X",1:48,"_r"), paste0("A",1:48,"_r")),
      "seq"=c(primerData$Fwseq,primerData$Rvseq))

write.csv(primerTable,"primerSeqs2order.csv",row.names=FALSE,quote=FALSE)




