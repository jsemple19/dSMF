#2017-05-16
#script to make GRanges and bed files from final amplicon data table
################## making gr object of final chosen amplicons
setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")
amplicons<-read.csv("finalChosenList_tss_dc_tissues.csv",stringsAsFactors=FALSE)

startEnd<-sapply(strsplit(as.character(amplicons$Amplicon),":",fixed=T), '[[',2)
chr<-sapply(strsplit(as.character(amplicons$Amplicon),":",fixed=T), '[[',1)
start<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',1))
end<-as.numeric(sapply(strsplit(startEnd,"-",fixed=T), '[[',2))
strand<-amplicons$strand
orientation<- ifelse(amplicons$orientation=="fw","+","-")

library(dplyr)
metadata<-data.frame(amplicons$publicID,amplicons$Gene_WB_ID,amplicons$primerNames,amplicons$dosageComp,
                     amplicons$orientation,amplicons$maxTSS,amplicons$most5p,amplicons$TSSwb,
                     amplicons$SpencerLarval11)#,
#amplicons$fragID, amplicons$NameFromBEDfile)
names(metadata)<-gsub("amplicons.","",names(metadata))
library(GenomicRanges)
GRamplicons<-GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand,
                     mcols=metadata)
names(mcols(GRamplicons))<-gsub("mcols.","",names(mcols(GRamplicons)))
names(GRamplicons)<-mcols(GRamplicons)$Gene_WB_ID
saveRDS(GRamplicons,"GRangesOfAmplicons.RDS")

################################333
amplicons<-readRDS("GRangesOfAmplicons.RDS")

### create BED file
amplicons.bed<-data.frame(chrom=chr,chromStart=start,chromEnd=end,name=amplicons$Gene_WB_ID,score=".",
                          strand=orientation)
write.table(amplicons.bed,file="amplicons_stranded.bed", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)


##### get seqs and create DNAstring object

GRamplicons<-readRDS("GRangesOfAmplicons.RDS")
mcols(GRamplicons)<-data.frame(chr=as.character(seqnames(GRamplicons)),start=start(GRamplicons),
                              end=end(GRamplicons), strand=as.character(strand(GRamplicons)),
                              mcols(GRamplicons))
library("BSgenome.Celegans.UCSC.ce11")
ampliconSeqs<-getSeq(BSgenome.Celegans.UCSC.ce11,GRamplicons)
mcols(ampliconSeqs)<-mcols(GRamplicons)
names(ampliconSeqs)<-mcols(GRamplicons)$Gene_WB_ID

saveRDS(ampliconSeqs,"ampliconsSeqs.RDS")


######### create GR of amplicon TSSs
setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")
amplicons<-read.csv("finalChosenList_tss_dc_tissues.csv",stringsAsFactors=FALSE)

GRamplicons<-readRDS("GRangesOfAmplicons.RDS")

GRampTSS<-GRamplicons
start(GRampTSS)<-mcols(GRamplicons)$maxTSS
end(GRampTSS)<-mcols(GRamplicons)$maxTSS

saveRDS(GRampTSS,"GRampliconTSS.RDS")
export.bed(GRampTSS,"ampliconTSS.bed")
