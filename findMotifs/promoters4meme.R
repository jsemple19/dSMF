#script to extract promoter sequences and write them as fasta files to use as input for meme and meme-chip

#################### writing promoter input files for MEME ############################
allTSS<-readRDS("../TSS/scripts/confirmedTSS_Chen&Kreus_WS235.Rds")

ss<-resize(allTSS,3,fix="center",ignore.strand=FALSE)
promoters<-resize(allTSS,101,fix="end",ignore.strand=FALSE)
promoters<-resize(promoters,151,fix="start",ignore.strand=FALSE)

library("BSgenome.Celegans.UCSC.ce11")
library(Biostrings)
genome<-BSgenome.Celegans.UCSC.ce11

ssSeqs<-getSeq(genome,ss)

promSeqs<-getSeq(genome,promoters)
names(promSeqs)<-mcols(promoters)$WBGeneID
i<-which(duplicated(names(promSeqs)))
promSeqs<-promSeqs[-i]

numProms=300
set.seed(numProms)
chosen<-sample(1:length(promSeqs),numProms)

# to write a minimal set for MEME
writeXStringSet(promSeqs[chosen[1]],paste0("./",numProms,"Promoters.fa"))
for(i in 2:numProms) {
   seqName<-names(promSeqs)[chosen[i]]
   writeXStringSet(promSeqs[chosen[i]],paste0("./",numProms,"Promoters.fa"),append=TRUE)
}

# to write full set for MEME-ChIP
writeXStringSet(promSeqs[1],paste0("./allPromoters.fa"))
for(i in 2:length(promSeqs)) {
   seqName<-names(promSeqs)[i]
   writeXStringSet(promSeqs[i],paste0("./allPromoters.fa"),append=TRUE)
}

