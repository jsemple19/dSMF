# 2017-05-08 modified to include better TSS choices
# 2016-12-20
# matchingExpn.R
# takes list of genes for which primers were designed on X chr and assigns them to quantiles of
# L3 gene expression in order to match to autosomal genes

setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")

#read in genes selected from ChrX
chosenX<-read.csv("./XdcPrimers/DM_DE_promoter_test_primers_WS250_subset.csv",header=TRUE,stringsAsFactors=FALSE)
chosenX<-chosenX[chosenX$finalSelected=="y","Gene_WB_ID"]

# read in original expression dataset to check it
#wormG<-read.csv("../../expressionMatching/worm_gene.csv")
#wormG<-wormG[,grep("L3",names(wormG))]
#cor(wormG$L3,wormG$N2_L3-1)
#cor(wormG$L3,wormG$L3_N2_L3-1)
#cor(wormG$N2_L3-1,wormG$L3_N2_L3-1)
#correlation strongest with N2_L3-1 and L3_N2_L3-1 probably becuase most of sequencing reads came
#from those runs? L3_N2_L3-1 is probably a weighted average of the two experiments.

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


#remove nonexpressed genes
expressed<-wormG[wormG$L3_N2_L3.1>0,]
#determine the cutoffs for gene expression to define deciles in this dataset
qt<-quantile(expressed$L3_N2_L3.1,probs=seq(0,1,0.1))
ranks<-sapply(expressed$L3_N2_L3.1,function(x) {max(which(x>qt))})

### now look at the distribution of my chosen genes on X among the deciles:
#filter table for genes with primers on X
chosenXexpn<-wormG[match(chosenX, wormG$wormbase_gene),]
#check which decile they fall into
ranksX<-sapply(chosenXexpn$L3_N2_L3.1,function(x) {max(which(x>qt))})
#look at the distribution
table(ranksX)
######old:
#decile:          3  4  5  6  7  8  9 10
#number of genes: 4  5  4 12  7  5  3  8
######new:
# decile:         3  4  5  6  7  8  9 10
#number of genes: 2  4  5  8  8  3  2 16

x<-barplot(table(ranksX),xlab="quantile of gene expression",ylab="number of genes",ylim=c(0,17),
           main="Distribution in expression quantiles of chosen X chromosome genes")
text(x,table(ranksX)+0.5,labels=table(ranksX))

### now see which deciles the autosomal genes fall into

#read in genes for which primers were designed on autosomes
chosenA<-read.table("./AdcPrimers/DM_DE_promoter_test_primers_WS250.tab",header=TRUE,stringsAsFactors=FALSE)

#subset the expression table to the autosomal genes with primers
i<-match(chosenA$Gene_WB_ID, expressed$wormbase_gene)
chosenAexpn<-expressed[i,]
#check which decile they fall into
ranksA<-sapply(chosenAexpn$L3_N2_L3.1,function(x) {max(which(x>qt))})
#look at the distribution
table(ranksA)
############old:
#decile:            -Inf  1    2    3    4    5    6    7    8    9   10
#number of genes:     7   13   36   40   67   90  124  144  199  192  236
###########new:
#decile:           -Inf    1    2    3    4    5    6    7    8    9   10
#number of genes:    14   29   51   48   93  111  151  172  233  257  352

x<-barplot(table(ranksA)[2:11],xlab="quantile of gene expression",ylab="number of genes",
           ylim=c(0,375),main="Distribution in expression quantiles of autosomal genes with primers")

text(x,table(ranksA)[2:11]+11,labels=table(ranksA)[2:11])


#save the ranks into the primer design table
chosenA<-cbind("decile"=ranksA,chosenA)

#### deal with problem that output .tab and .bed have different numbers of rows (amplicons
# that are rejected??):
bedA<-read.delim("./AdcPrimers/test_primers_bed_WS250.bed",stringsAsFactors=FALSE,header=FALSE)
#extract wb gene id from primer id column
bedA<-data.frame("Gene_WB_ID"=sapply(strsplit(bedA$V4,".",fixed=T), '[[',2),
                 "location"=paste0(bedA$V1,":",bedA$V2,"-",bedA$V3),
                 bedA, stringsAsFactors=FALSE)
bedA<-bedA[bedA$Gene_WB_ID %in% chosenA$Gene_WB_ID,]
bedA<-bedA[order(match(bedA$Gene_WB_ID,chosenA$Gene_WB_ID)),]
#add bedfile data for location (and names to double check sorting):
chosenA$nameInBedFile<-bedA$V4
chosenA$location<-bedA$location
chosenA$finalSelected<-""
#reorder columns
chosenA<-chosenA[,c(19,1:3,17:18,4:16)]

#write in which ones have already been selected (just need to adjust numbers in categories
#not rechoose them all)
oldChosenA<-read.csv("./AdcPrimers_old/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                                   header=TRUE,stringsAsFactors=FALSE)
oldChosenA<-oldChosenA[oldChosenA$FinalSelected=="y","Gene_WB_ID"]
chosenA$finalSelected[chosenA$Gene_WB_ID %in% oldChosenA]<-"y"

write.csv(chosenA,"./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks.csv")


##############3 after doing manual selection: check for similar expression patterns
finalChosenA<-read.csv("./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                       header=TRUE,stringsAsFactors=FALSE)
finalChosenA<-finalChosenA[finalChosenA$finalSelected=="y","Gene_WB_ID"]

expnChosen<-rbind(cbind("chr"=rep("chrX",48),expressed[expressed$wormbase_gene %in% chosenX,]),
                  cbind("chr"=rep("Autosomes",48),expressed[expressed$wormbase_gene %in% finalChosenA,]))

library(ggplot2)
library(reshape)
#ggplot(expnChosen[,c(1,3,7,8,10:17)],aes(x="chr",y="PolII_EE_FE"))+geom_boxplot()

pdf(file="./primerQC.pdf", width=11,height=8,paper="a4r")

mm=melt(expnChosen[,c(1,3,7,8,10:17)],id=c("wormbase_gene","chr"))
names(mm)<-c("wormbase_gene","chr","dataset","expression")
ggplot(mm,aes(x=chr,y=log2(expression),fill=chr)) + geom_boxplot() + facet_wrap(~dataset,nrow=2)


############### compare other attirbutes of primers
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
mm=melt(primerData[,c(1,2,7:17)],id=c("NameFromBEDfile","chr"))
ggplot(mm,aes(x=chr,y=value,fill=chr)) + geom_boxplot() +
    facet_wrap(~variable,scales="free")

mm=melt(primerData[,c(1,2,9:17)],id=c("NameFromBEDfile","chr"))
ggplot(mm,aes(value,fill=chr)) + geom_histogram() +
  facet_wrap(~chr+variable,scales="free",nrow=2)

dev.off()

table(primerData$FwC.covered+primerData$RvC.covered)
#0  1  2  3  4  5  6  7  8  9 10
#6  8 18 12 22 11  7  5  4  1  2

primerTable<-data.frame("primerNames"=c(paste0("X",1:48,"_f"),paste0("A",1:48,"_f"), paste0("X",1:48,"_r"), paste0("A",1:48,"_r")),
      "seq"=c(primerData$Fwseq,primerData$Rvseq))

write.csv(primerTable,"primerSeqs2order.csv",row.names=FALSE,quote=FALSE)


#### get length of Amplicons to check PCR results
AmpliconLengths<-data.frame("primerNames"=c(paste0("X",1:48,"_f"),paste0("A",1:48,"_f"), paste0("X",1:48,"_r"), paste0("A",1:48,"_r")),
                               "seq"=c(primerData$Fwseq,primerData$Rvseq),"amplicon"=c(primerData$Amplicon,primerData$Amplicon),
                               "rowNum"=rep(rep(c("A","B","C","D","E","F","G","H"), each=12),2),"colNum"=rep(rep(1:12,8),2))
startEnd<-sapply(strsplit(as.character(AmpliconLengths$amplicon),":",fixed=T), '[[',2)
start<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',1))
end<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',2))
AmpliconLengths["width"]<-end-start

AmpliconLengths<-AmpliconLengths[c(1:96),]
AmpliconLengths
#(8+10+9+6+9+11)/96


interleave<-function(v1,v2) {
   ord1<-2*(1:length(v1))-1
   ord2<-2*(1:length(v2))
   c(v1,v2)[order(c(ord1,ord2))]
}

printColPair_plate<-function(data,col1,col2) {
   r1=data[data$colNum==col1,]
   r2=data[data$colNum==col2,]
   out<-interleave(r1$width[order(r1$rowNum,decreasing=T)],r2$width[order(r2$rowNum,decreasing=T)])
   names(out)<-rep(r1$rowNum[order(r1$rowNum,decreasing=T)],each=2)
   return(out)
}

printColPair_plate(AmpliconLengths,1,2)
printColPair_plate(AmpliconLengths,3,4)
printColPair_plate(AmpliconLengths,5,6)
printColPair_plate(AmpliconLengths,7,8)
printColPair_plate(AmpliconLengths,9,10)
printColPair_plate(AmpliconLengths,11,12)

printColPair<-function(data,col1,col2) {
   r1=data[data$colNum==col1,]
   r2=data[data$colNum==col2,]
   out<-interleave(r1$width[order(r1$rowNum,decreasing=F)],r2$width[order(r2$rowNum,decreasing=F)])
   names(out)<-rep(r1$rowNum[order(r1$rowNum,decreasing=F)],each=2)
   return(out)
}

printColPair(AmpliconLengths,1,2)
printColPair(AmpliconLengths,3,4)
printColPair(AmpliconLengths,5,6)
printColPair(AmpliconLengths,7,8)
printColPair(AmpliconLengths,9,10)
printColPair(AmpliconLengths,11,12)



### make final list of chosen ############
chosenX<-read.csv("./XdcPrimers/DM_DE_promoter_test_primers_WS250_subset.csv",header=TRUE,stringsAsFactors=FALSE)
chosenX<-chosenX[chosenX$finalSelected=="y",]
#names(chosenX)<-c(names(chosenX)[1],"fragID",names(chosenX)[2:18])
#names(chosenX)[5]<-"Amplicon"
names(chosenX)[6]<-"orientation"
#txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% chosenX$Gene_WB_ID,]
i<-match(chosenX$Gene_WB_ID,mcols(genes)$gene_id)
chosenX$strand<-as.vector(strand(genes[i,]))


finalChosenA<-read.csv("./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                       header=TRUE,stringsAsFactors=FALSE)
finalChosenA<-finalChosenA[finalChosenA$finalSelected=="y",]
names(finalChosenA)[6]<-"NameFromBEDfile"
names(chosenX)[2]<-"fragID"
#names(finalChosenA)[5]<-"PrimerID"
names(finalChosenA)[7]<-"Amplicon"
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% finalChosenA$Gene_WB_ID,]
i<-match(finalChosenA$Gene_WB_ID,mcols(genes)$gene_id)
finalChosenA$strand<-as.vector(strand(genes[i,]))


chosenX<-data.frame(chosenX[,c(1:2)],"decile"=ranksX,"strand"=chosenX$strand,chosenX[,c(3:(dim(chosenX)[2]-1))])
finalChosenA<-data.frame(finalChosenA[,c(1:3)],"strand"=finalChosenA$strand,finalChosenA[,c(4,6,7,5,8:(dim(finalChosenA)[2]-1))])

chosen<-rbind(chosenX,finalChosenA)

primerNames=paste0(c(paste0("X",1:48,"_f"),paste0("A",1:48,"_f")), "_&_",c(paste0("X",1:48,"_r"), paste0("A",1:48,"_r")))

chosen<-data.frame(primerNames,chosen[,c(2:(dim(chosen)[2]))])

source("/home/jenny/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction.R")
conversionTable<-convertGeneNames(chosen$Gene_WB_ID,inputType="wormbase_gene")

i<-match(chosen$Gene_WB_ID,conversionTable$wormbase_gene)
chosen<-data.frame(chosen[,c(1:4)],"publicID"=conversionTable[i,1],chosen[,c(5:(dim(chosen)[2]))])
write.csv(chosen,"finalChosenList.csv",quote=FALSE,row.names=FALSE)


################################333
################## making gr object of final chosen amplicons
amplicons<-read.csv("finalChosenList.csv",stringsAsFactors=FALSE)
#amplicons<-read.csv("finalChosenList_tss_dc.csv",stringsAsFactors=FALSE)
startEnd<-sapply(strsplit(as.character(amplicons$Amplicon),":",fixed=T), '[[',2)
chr<-sapply(strsplit(as.character(amplicons$Amplicon),":",fixed=T), '[[',1)
start<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',1))
end<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',2))
strand<-amplicons$strand
orientation<- ifelse(amplicons$orientation=="fw","+","-")
amplicons$dosageComp<-rep(c("dc","ndc"),each=48)
amplicons$dosageComp[grep("ndc",amplicons$fragID)]<-"ndc"
amplicons$dosageComp[grep("_dc",amplicons$fragID)]<-"dc"



#y<-readRDS("GRangesOfAmplicons.RDS")

### create BED file
amplicons.bed<-data.frame(chrom=chr,chromStart=start,chromEnd=end,name=amplicons$Gene_WB_ID,score=".",
                          strand=orientation)
write.table(amplicons.bed,file="amplicons.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

tss<-readRDS("../../TSS/scripts/confirmedTSS_Chen&Kreus_WS235.Rds")
tss<-tss[mcols(tss)$WBGeneID %in% amplicons$Gene_WB_ID,]

i<-match(amplicons$Gene_WB_ID,tss$WBGeneID)
amplicons$TSSdata<-start(tss[i,])

# compare to annotated TSS
library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% amplicons$Gene_WB_ID,]
#tss<-sort(tss)
#genes<-sort(genes)
i<-match(amplicons$Gene_WB_ID,mcols(genes)$gene_id)

amplicons$TSSwb<-ifelse(strand(genes[i])=="+",start(genes[i]),end(genes[i]))
#test start positions
sum(ifelse(amplicons$strand=="+",amplicons$TSSwb-amplicons$TSSdata,amplicons$TSSdata-amplicons$TSSwb)<0)
# in 30 cases the annotated is upstream of the data TSS. replace those is a "most5p" column
amplicons$most5p<-amplicons$TSSdata
i<-which(ifelse(amplicons$strand=="+",amplicons$TSSwb-amplicons$TSSdata,amplicons$TSSdata-amplicons$TSSwb)<0)
amplicons$most5p[i]<-amplicons$TSSwb[i]

write.csv(amplicons,"finalChosenList_tss_dc.csv",quote=FALSE,row.names=FALSE)

amplicons<-read.csv("finalChosenList_tss_dc.csv",stringsAsFactors=FALSE)

tss.bed<-data.frame(chrom=chr,chromStart=amplicons$TSSdata,chromEnd=amplicons$TSSdata,
                    name=amplicons$Gene_WB_ID,score=0,strand=strand)
#tss.bed<-data.frame(chrom=chr,chromStart=start(tss),chromEnd=end(tss),name=mcols(tss)$WBGeneID,score=".",strand=strand(tss))
write.table(tss.bed,file="tss.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

most5p.bed<-data.frame(chrom=chr,chromStart=amplicons$most5p,chromEnd=amplicons$most5p,
                       name=amplicons$Gene_WB_ID,score=0,strand=strand)
#genes.bed<-data.frame(chrom=as.vector(seqnames(genes)),chromStart=start(genes),chromEnd=end(genes),name=mcols(genes)$gene_id,score=".",strand=as.vector(strand(genes)))
write.table(most5p.bed,file="most5p.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)





library(dplyr)
metadata<-data.frame(amplicons$publicID,amplicons$Gene_WB_ID,amplicons$primerNames,amplicons$dosageComp,
                     amplicons$orientation,amplicons$TSSdata,amplicons$TSSwb,amplicons$most5p)#,
#amplicons$fragID, amplicons$NameFromBEDfile)
names(metadata)<-gsub("amplicons.","",names(metadata))
library(GenomicRanges)
GRamplicons<-GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand,
                     mcols=metadata)
saveRDS(GRamplicons,"GRangesOfAmplicons.RDS")




source("/home/jenny/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction.R")
conversionTable<-convertGeneNames(mcols(tss)$gene_id,inputType="wormbase_gene_seq_name")

# txpts<-transcripts(txdb)
# txpts.bed<-data.frame(chrom=as.vector(seqnames(txpts)),chromStart=start(txpts),chromEnd=end(txpts),name=mcols(txpts)$tx_name,score=".",strand=as.vector(strand(txpts)))
# write.table(txpts.bed,file="txpts.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# exons<-exons(txdb)
# exons.bed<-data.frame(chrom=as.vector(seqnames(exons)),chromStart=start(exons),chromEnd=end(exons),name=mcols(exons)$exon_id,score=".",strand=as.vector(strand(exons)))
# write.table(exons.bed,file="exons.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)


#read in tissue data from "Comparative Analysis of the Transcriptome across Distant Species" doi:10.1038/nature13424

exprn<-read.csv("../../expressionMatching/worm_gene_tissues_compressed.csv",header=TRUE, stringsAsFactors=FALSE)
#exprnAmplicons<-exprn[exprn$Gene %in% mcols(tss)$assignedGeneName,]

i<-match(mcols(tss)$assignedGeneName, exprn$Gene)

exprnAmplicons<-exprn[i,]
exprnOther<-exprn[-i,]

boxplot(exprnAmplicons$BroadlyExpr_Score,exprnOther$BroadlyExpr_Score,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, main="BroadlyExpr_Score",col="light blue")
abline(h=median(exprnOther$BroadlyExpr_Score,na.rm=TRUE),lty=2)

#MEsamples<-read.csv("../../expressionMatching/modEncode_samples_compressed.csv",stringsAsFactors=FALSE,header=TRUE)

colMedians<-apply(exprn[,3:15],2,median)

numTissuesAmplicons<-apply(exprnAmplicons[,3:15]>colMedians,1,sum)

numTissuesOthers<-apply(exprnOther[,3:15]>colMedians,1,sum)

boxplot(numTissuesAmplicons,numTissuesOthers,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, main="numTissues",col="light blue")

uniqueTissues<-c(7,8,12)
colMedians<-apply(exprn[,uniqueTissues],2,median)

numTissuesAmplicons<-apply(exprnAmplicons[,uniqueTissues]>colMedians,1,sum)

numTissuesOthers<-apply(exprnOther[,uniqueTissues]>colMedians,1,sum)

boxplot(numTissuesAmplicons,numTissuesOthers,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, main="numTissues",col="light blue")

par(mfrow=c(1,2))
barplot(table(numTissuesAmplicons),main="Number of tissues - Amplicons",xlab="number of tissues",ylab="number of genes")
barplot(table(numTissuesOthers),main="Number of tissues - OtherGenes",xlab="number of tissues",ylab="number of genes")


flp27<-"C25H3.5"
