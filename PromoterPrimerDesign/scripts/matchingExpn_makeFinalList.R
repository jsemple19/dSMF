#2017-05-15
#script to continue processing amplicon table and add dc status, TSS and expression
#breadth data


setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")
## make final list of chosen ############
chosenX<-read.csv("./XdcPrimers/DM_DE_promoter_test_primers_WS250_subset.csv",header=TRUE,stringsAsFactors=FALSE)
chosenX<-chosenX[chosenX$finalSelected=="y",]
#names(chosenX)<-c(names(chosenX)[1],"fragID",names(chosenX)[2:18])
#names(chosenX)[5]<-"Amplicon"
names(chosenX)[6]<-"orientation"

library(GenomicFeatures)
if(!exists("txdb")) {
   txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
}
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

################## recreate ranksX object##################
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
#determine the cutoffs for gene expression to define deciles in this dataset
qt<-quantile(expressed$L3_N2_L3.1,probs=seq(0,1,0.1))
ranks<-sapply(expressed$L3_N2_L3.1,function(x) {max(which(x>qt),na.rm=TRUE)})
ranks[which(ranks=="-Inf")]<-1 #deal with -Inf for low values equal to lowest quantile

### now look at the distribution of my chosen genes on X among the deciles:
#filter table for genes with primers on X
chosenXexpn<-wormG[match(chosenX$Gene_WB_ID, wormG$WBgeneID),]
#check which decile they fall into
ranksX<-sapply(chosenXexpn$L3_N2_L3.1,function(x) {max(which(x>qt))})
########################

chosenX<-data.frame(chosenX[,c(1:2)],"decile"=ranksX,"strand"=chosenX$strand,chosenX[,c(3:(dim(chosenX)[2]-1))])
finalChosenA<-data.frame(finalChosenA[,c(1:3)],"strand"=finalChosenA$strand,finalChosenA[,c(4,6,7,5,8:(dim(finalChosenA)[2]-1))])

chosen<-rbind(chosenX,finalChosenA)

primerNames=paste0(c(paste0("X",1:48,"_f"),paste0("A",1:48,"_f")), "_&_",c(paste0("X",1:48,"_r"), paste0("A",1:48,"_r")))

chosen<-data.frame(primerNames,chosen[,c(2:(dim(chosen)[2]))])

source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
conversionTable<-convertGeneNames(chosen$Gene_WB_ID,inputType="WBgeneID")

i<-match(chosen$Gene_WB_ID,conversionTable$WBgeneID)
chosen<-data.frame(chosen[,c(1:4)],"publicID"=conversionTable[i,"publicID"],chosen[,c(5:(dim(chosen)[2]))])
write.csv(chosen,"finalChosenList.csv",quote=FALSE,row.names=FALSE)

##################### add DC classification and TSS data #######################
amplicons<-read.csv("finalChosenList.csv",stringsAsFactors=FALSE)
amplicons$orientation<- ifelse(amplicons$orientation=="fw","+","-")
amplicons$dosageComp<-rep(c("dc","ndc"),each=48)
amplicons$dosageComp[grep("ndc",amplicons$fragID)]<-"ndc"
amplicons$dosageComp[grep("_dc",amplicons$fragID)]<-"dc"


#load TSS data
TSS<-read.table("~/Documents/MeisterLab/dSMF/TSS/scripts/ChenKreusSaitoTSS_2827_maxTSS.txt",stringsAsFactors=FALSE)
colnames(TSS)<-c("WBGeneID","chr","start","end","strand","chenTSS","kruesiTSS","saitoTSS","maxTSS", "most5p")

i<-match(amplicons$Gene_WB_ID,TSS$WBGeneID)
amplicons<-cbind(amplicons,TSS[i,c("chenTSS","kruesiTSS","saitoTSS","maxTSS", "most5p")])

# compare to annotated TSS
library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
if(!exists("txdb")) {
   txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
}
genes<-genes(txdb)
genes<-genes[mcols(genes)$gene_id %in% amplicons$Gene_WB_ID,]

i<-match(amplicons$Gene_WB_ID,mcols(genes)$gene_id)

amplicons$maxTSS<-as.numeric(amplicons$maxTSS)
amplicons$TSSwb<-ifelse(strand(genes[i])=="+",start(genes[i]),end(genes[i]))
#test start positions
sum(ifelse(amplicons$strand=="+",amplicons$TSSwb-amplicons$maxTSS,amplicons$maxTSS-amplicons$TSSwb)<0)
# in 34 cases the annotated is upstream of the data TSS.

write.csv(amplicons,"finalChosenList_tss_dc.csv",quote=FALSE,row.names=FALSE)



############################ add expression breadth data ################

source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
conversionTable<-convertGeneNames(amplicons$Gene_WB_ID,inputType="WBgeneID")

#read in tissue data from "Comparative Analysis of the Transcriptome across Distant
#Species" doi:10.1038/nature13424. Gerstein et al

exprn<-read.csv("../../expressionMatching/worm_gene_tissues_compressed.csv",
                header=TRUE, stringsAsFactors=FALSE)
#add chromosome info to genes
exprn$WBgeneID<-convertGeneNames(exprn$Gene,inputType="seqID",outputType="WBgeneID")
genes<-genes(txdb)
exprn<-exprn[exprn$WBgeneID %in% names(genes),]
exprn$chr<-as.character(seqnames(genes)[na.omit(match(exprn$WBgeneID,names(genes)))])

i<-match(conversionTable$seqID, exprn$Gene)

exprnAmplicons<-exprn[i,]
exprnOther<-exprn[-i,]

exprnXamp<-exprnAmplicons[exprnAmplicons$chr=="X",]
exprnAamp<-exprnAmplicons[exprnAmplicons$chr!="X",]
exprnXall<-exprnOther[exprnOther$chr=="X",]
exprnAall<-exprnOther[exprnOther$chr!="X",]

pdf("plotsExprnBreadth.pdf",paper="a4",height=11,width=8)

par(mfrow=c(2,1))
boxplot(exprnAmplicons$BroadlyExpr_Score,
        exprnOther$BroadlyExpr_Score,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, main="Broadly Expressed Score",col="light blue")
abline(h=median(exprnOther$BroadlyExpr_Score,na.rm=TRUE),lty=2)

boxplot(exprnXamp$BroadlyExpr_Score, exprnAamp$BroadlyExpr_Score,
        exprnXall$BroadlyExpr_Score, exprnAall$BroadlyExpr_Score,
        names=paste(rep(c("amplicons","other genes"),each=2),
                    rep(c("chrX","Autosome"),2),sep=" "), cex.axis=0.7,
        main="Broadly Expressed Score", col=rep(c("purple","blue"),2))
abline(h=median(exprnAall$BroadlyExpr_Score,na.rm=TRUE),lty=2)


#MEsamples<-read.csv("../../expressionMatching/modEncode_samples_compressed.csv",stringsAsFactors=FALSE,header=TRUE)

colMedians<-apply(exprn[,3:15],2,median)

numTissuesAmplicons<-apply(exprnAmplicons[,3:15]>colMedians,1,sum)
numTissuesOthers<-apply(exprnOther[,3:15]>colMedians,1,sum)

boxplot(numTissuesAmplicons,numTissuesOthers,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, main="All Tissues (Gerstein)",col="light blue")

numTissuesAmpX<-apply(exprnXamp[,3:15]>colMedians,1,sum)
numTissuesAmpA<-apply(exprnAamp[,3:15]>colMedians,1,sum)
numTissuesAllX<-apply(exprnXall[,3:15]>colMedians,1,sum)
numTissuesAllA<-apply(exprnAall[,3:15]>colMedians,1,sum)

boxplot(numTissuesAmpX, numTissuesAmpA, numTissuesAllX, numTissuesAllA,
        names=paste(rep(c("amplicons","other genes"),each=2),
                    rep(c("chrX","Autosome"),2),sep=" "), cex.axis=0.7,
        main="All Tissus (Gerstein)", col=rep(c("purple","blue"),2))
abline(h=median(numTissuesAllA,na.rm=TRUE),lty=2)


uniqueTissues<-c(7,8,12)
colMedians<-apply(exprn[,uniqueTissues],2,median)

numTissuesAmplicons<-apply(exprnAmplicons[,uniqueTissues]>colMedians,1,sum)
numTissuesOthers<-apply(exprnOther[,uniqueTissues]>colMedians,1,sum)

boxplot(numTissuesAmplicons,numTissuesOthers,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, ylab="Tissues with expression",
        main="Unique Tissues (Gerstein)",col="light blue")

numTissuesAmpX<-apply(exprnXamp[,uniqueTissues]>colMedians,1,sum)
numTissuesAmpA<-apply(exprnAamp[,uniqueTissues]>colMedians,1,sum)
numTissuesAllX<-apply(exprnXall[,uniqueTissues]>colMedians,1,sum)
numTissuesAllA<-apply(exprnAall[,uniqueTissues]>colMedians,1,sum)

boxplot(numTissuesAmpX, numTissuesAmpA, numTissuesAllX, numTissuesAllA,
        names=paste(rep(c("amplicons","other genes"),each=2),
                    rep(c("chrX","Autosome"),2),sep=" "), cex.axis=0.7,
        main="Unique Tissus (Gerstein)", col=rep(c("purple","blue"),2),
        ylab="Tissues with expression")
abline(h=median(numTissuesAllA,na.rm=TRUE),lty=2)


par(mfrow=c(2,2))
barplot(table(numTissuesAmplicons),main="Number of tissues - Amplicons",
        xlab="number of tissues",ylab="number of genes")
barplot(table(numTissuesOthers),main="Number of tissues - OtherGenes",
        xlab="number of tissues",ylab="number of genes")

amplicons<-cbind(amplicons,Gerstein3unique=numTissuesAmplicons)

######### read in Spencer2012 expn breadth data
expnBreadth<-readRDS("../../tissueExpn/scripts/tissueExpressionBreadth.RDS")


#add chromosome info to genes
expnBreadth<-expnBreadth[expnBreadth$WBgeneID %in% names(genes),]
expnBreadth$chr<-as.character(seqnames(genes)[na.omit(match(expnBreadth$WBgeneID,names(genes)))])

i<-match(amplicons$Gene_WB_ID,expnBreadth$WBgeneID)

ampExpnBreadth<-expnBreadth[i,]
allExpnBreadth<-expnBreadth[-i,]
#subset by chromosome location
expbXamp<-ampExpnBreadth[ampExpnBreadth$chr=="X",]
expbAamp<-ampExpnBreadth[ampExpnBreadth$chr!="X",]
expbXall<-allExpnBreadth[allExpnBreadth$chr=="X",]
expbAall<-allExpnBreadth[allExpnBreadth$chr!="X",]

par(mfrow=c(2,1))
boxplot(ampExpnBreadth$larvalTissues,allExpnBreadth$larvalTissues,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, ylab="tissues with expression",
        main="Larval tissues (Spencer)",col="light blue")
abline(h=median(expnBreadth$larvalTissues))
boxplot(ampExpnBreadth$allTissues,allExpnBreadth$allTissues,varwidth=TRUE,
        names=c("amplicons","other genes"),notch=TRUE, ylab="tissues with expression",
        main="All tissues (Spencer)",col="light blue")
abline(h=median(expnBreadth$allTissues))

boxplot(expbXamp$larvalTissues, expbAamp$larvalTissues,
        expbXall$larvalTissues, expbAall$larvalTissues,
        names=paste(rep(c("amplicons","other genes"),each=2),
                    rep(c("chrX","autosome"),2),sep=" "), cex.axis=0.7,
        main="Larval tissues (Spencer)", col=rep(c("purple","blue"),2))
abline(h=median(expbAall$larvalTissues,na.rm=TRUE),lty=2)

boxplot(expbXamp$allTissues, expbAamp$allTissues,
        expbXall$allTissues, expbAall$allTissues,
        names=paste(rep(c("amplicons","other genes"),each=2),
                    rep(c("chrX","autosome"),2),sep=" "), cex.axis=0.7,
        main="All tissues (Spencer)", col=rep(c("purple","blue"),2))
abline(h=median(expbAall$allTissues,na.rm=TRUE),lty=2)


par(mfrow=c(2,2))
barplot(table(ampExpnBreadth$larvalTissues),main="Larval tissues (Spencer) - Amplicons",
        xlab="number of tissues",ylab="number of genes",cex.main=1)
barplot(table(allExpnBreadth$larvalTissues),main="Larval tissues (Spencer) - OtherGenes",
        xlab="number of tissues",ylab="number of genes",cex.main=1)
barplot(table(ampExpnBreadth$allTissues),main="All tissues (Spencer) - Amplicons",
        xlab="number of tissues",ylab="number of genes",cex.main=1)
barplot(table(allExpnBreadth$allTissues),main="All tissues (Spencer) - OtherGenes",
        xlab="number of tissues",ylab="number of genes",cex.main=1)
dev.off()

amplicons<-cbind(amplicons,SpencerLarval11=ampExpnBreadth$larvalTissues,
      SpencerAll25=ampExpnBreadth$allTissues,fromUnfilt=ampExpnBreadth$fromUnfilt)

write.csv(amplicons,"finalChosenList_tss_dc_tissues.csv",quote=FALSE,row.names=FALSE)
#######################
