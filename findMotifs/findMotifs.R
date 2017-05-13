setwd("~/Documents/MeisterLab/dSMF/findMotifs")
amplicons<-readRDS("../PromoterPrimerDesign/scripts/ampliconsSeqs.RDS")

library(MotifDb)
library(seqLogo)
library(MotIV)
values(MotifDb)$organism
head(MotifDb)

ce<-query(MotifDb,'celegans')
tata<-query(MotifDb,'TATA')
sp1<-query(MotifDb,'SP1')
inr<-query(MotifDb,'Inr')
#from Jaspar http://jaspar.genereg.net/cgi-bin/jaspar_db.pl
inr<-rbind(c( 49,   0, 288,  26,  77,  67,  45,  50),
       c( 48, 303,   0,  81,  95, 118,  85,  96),
       c(69,   0,   0, 116,   0,  46,  73,  56),
       c(137,   0,  15,  80, 131,  72, 100, 101))
row.names(inr)=c("A","C","G","T")
colnames(inr)=c(1:8)
inr<-inr/colSums(inr)
seqLogo(inr)

seqLogo(as.list(tata)[[5]])
tataMotif<-tata[[5]]
tataRel2wbTSS<-c()
for (i in 1:length(amplicons)) {
   i=84
   matches<-matchPWM(tataMotif,amplicons[[i]])
   metaData<-mcols(amplicons[i])
   ampStrand=as.character(metaData$strand)
   reldataTSS<-abs(getRelTSS(ampStrand,metaData$TSSdata, metaData$start,metaData$end))+1
   relPosMotif2dataTSS<-getPosRel2TSS(ampStrand,reldataTSS,start(matches))
   #tataRel2dataTSS<-rbind(tataRel2dataTSS,getPosRel2TSS(ampStrand,reldataTSS,start(matches))
   relwbTSS<-abs(getRelTSS(ampStrand,metaData$TSSwb, metaData$start,metaData$end))+1
   relPosMotif2wbTSS<-getPosRel2TSS(ampStrand,relwbTSS,start(matches))
   i<-which.min(abs(-40-c(relPosMotif2dataTSS,relPosMotif2wbTSS)))
   chosenTata<-
}

findPositionalPattern<-function(TSS, expPos, width,pattern) {
   rangePos<-resize(TSS+expPos,width,fix="center",ignore.strand=FALSE)
   #library("BSgenome.Celegans.UCSC.ce11")
   #library(Biostrings)
   genome<-BSgenome.Celegans.UCSC.ce11
   mySeqs<-getSeq(genome,rangePos)
   posPattern<-vmatchPattern(pattern,mySeqs)
   return(posPattern)
}
findPositionalPattern(allTSS,0,width=3,inr)

getRelTSS<-function(strand,tss,ampStart,ampEnd) {
   if (strand=="+"){
      relDistance<-tss-ampStart
   } else {
      relDistance<-ampEnd-tss
   }
   return(relDistance)
}

getPosRel2TSS<-function(strand,tss,pos) {
   if (strand=="+"){
      relDistance<-pos-tss
   } else {
      relDistance<-tss-pos
   }
   return(relDistance)
}

seqLogo(as.list(inr)[1])
inr<-DNAString("CA")
for (i in 1:96) {
   a<-matchPattern(inr,amplicons[[i]])
   print(a)
}

allTSS<-readRDS("../TSS/scripts/confirmedTSS_Chen&Kreus_WS235.Rds")

promoters<-resize(allTSS,101,fix="end",ignore.strand=FALSE)
promoters<-resize(promoters,151,fix="start",ignore.strand=FALSE)

library("BSgenome.Celegans.UCSC.ce11")
library(Biostrings)
genome<-BSgenome.Celegans.UCSC.ce11

promSeqs<-getSeq(genome,promoters)
names(promSeqs)<-mcols(promoters)$WBGeneID
i<-which(duplicated(names(promSeqs)))
promSeqs<-promSeqs[-i]

numProms=300
set.seed(numProms)
chosen<-sample(1:length(promSeqs),numProms)
writeXStringSet(promSeqs[chosen[1]],paste0("./",numProms,"Promoters.fa"))
for(i in 2:numProms) {
   seqName<-names(promSeqs)[chosen[i]]
   writeXStringSet(promSeqs[chosen[i]],paste0("./",numProms,"Promoters.fa"),append=TRUE)
}


writeXStringSet(promSeqs[1],paste0("./allPromoters.fa"))
for(i in 2:length(promSeqs)) {
   seqName<-names(promSeqs)[i]
   writeXStringSet(promSeqs[i],paste0("./allPromoters.fa"),append=TRUE)
}

