#2017-05-16
#script to make GRanges and bed files from final amplicon data table

amplicons<-read.csv("finalChosenList_tss_dc_tissues.csv",stringsAsFactors=FALSE)

tss.bed<-data.frame(chrom=chr,chromStart=amplicons$TSSdata,chromEnd=amplicons$TSSdata,
                    name=amplicons$Gene_WB_ID,score=0,strand=strand)
#tss.bed<-data.frame(chrom=chr,chromStart=start(tss),chromEnd=end(tss),name=mcols(tss)$WBGeneID,score=".",strand=strand(tss))
write.table(tss.bed,file="tss.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

most5p.bed<-data.frame(chrom=chr,chromStart=amplicons$most5p,chromEnd=amplicons$most5p,
                       name=amplicons$Gene_WB_ID,score=0,strand=strand)
#genes.bed<-data.frame(chrom=as.vector(seqnames(genes)),chromStart=start(genes),chromEnd=end(genes),name=mcols(genes)$gene_id,score=".",strand=as.vector(strand(genes)))



write.table(most5p.bed,file="most5p.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

################################333
################## making gr object of final chosen amplicons
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


library(dplyr)
metadata<-data.frame(amplicons$publicID,amplicons$Gene_WB_ID,amplicons$primerNames,amplicons$dosageComp,
                     amplicons$orientation,amplicons$TSSdata,amplicons$TSSwb,amplicons$most5p)#,
#amplicons$fragID, amplicons$NameFromBEDfile)
names(metadata)<-gsub("amplicons.","",names(metadata))
library(GenomicRanges)
GRamplicons<-GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand,
                     mcols=metadata)
saveRDS(GRamplicons,"GRangesOfAmplicons.RDS")





