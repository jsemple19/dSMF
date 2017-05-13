setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")
library("BSgenome.Celegans.UCSC.ce11")

###### get amplicon seqs with gene strand (not orientation strand)
#if you want orientation strand look at getAmpliconSeqs_4methDigest.R
amplicons<-readRDS("GRangesOfAmplicons.RDS")
names(mcols(amplicons))<-gsub("mcols.","",names(mcols(amplicons)))

genome<-BSgenome.Celegans.UCSC.ce11

mySeqs<-getSeq(genome,amplicons)#[i,])
#names(mySeqs)<-whichGenes
names(mySeqs)<-mcols(amplicons)$Gene_WB_ID
mcols(mySeqs)<-DataFrame(chr=seqnames(amplicons), start=start(amplicons), 
                        end=end(amplicons), strand=strand(amplicons),
                        mcols(amplicons))
saveRDS(mySeqs,"ampliconsSeqs.RDS")

writeXStringSet(mySeqs[1],paste0("./ampliconSeqs/allAmplicons.fa"))
for(i in 2:length(mySeqs)) {
   seqName<-names(mySeqs)[i]
   writeXStringSet(mySeqs[i],paste0("./ampliconSeqs/allAmplicons.fa"),append=TRUE)
}
