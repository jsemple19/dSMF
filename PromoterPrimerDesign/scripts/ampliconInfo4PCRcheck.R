#script to get list of amplicons, their length and provide an output in the order they appear on a gel to help
# with checking amplicon PCR results

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

