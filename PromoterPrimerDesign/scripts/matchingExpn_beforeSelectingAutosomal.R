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
i<-match(chosenA$Gene_WB_ID, expressed$WBgeneID)
chosenAexpn<-expressed[i,]
#check which decile they fall into
ranksA<-sapply(chosenAexpn$L3_N2_L3.1,function(x) {max(which(x>qt))})
ranksA[ranksA=="-Inf"]<-1
#look at the distribution
table(ranksA)
############old:
#decile:            -Inf  1    2    3    4    5    6    7    8    9   10
#number of genes:     7   13   36   40   67   90  124  144  199  192  236
###########new:
#decile:           -Inf    1    2    3    4    5    6    7    8    9   10
#number of genes:    14   29   51   48   93  111  151  172  233  257  352
#removeing -INF
ranksA
#1   2   3   4   5   6   7   8   9  10
#43  51  48  93 111 151 173 232 258 351

x <-barplot(table(ranksA)[2:11],xlab="quantile of gene expression",ylab="number of genes",
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
oldChosenA<-read.csv("./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks_subset.csv",
                                   header=TRUE,stringsAsFactors=FALSE)
oldChosenA<-oldChosenA[oldChosenA$FinalSelected=="y","Gene_WB_ID"]
chosenA$finalSelected[chosenA$Gene_WB_ID %in% oldChosenA]<-"y"

write.csv(chosenA,"./AdcPrimers/DM_DE_promoter_test_primers_WS250_withRanks.csv")


