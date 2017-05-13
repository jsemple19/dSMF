setwd("~/Documents/MeisterLab/dSMF/findMotifs")
library(MotifDb)
library(seqLogo)
library(GenomicRanges)
library(rtracklayer)
#library(MotIV)

source("./findMotif_functions.R")

# read in amplicons
amplicons<-readRDS("../PromoterPrimerDesign/scripts/ampliconsSeqs.RDS")

# get vectors of TSSs with different definitions
relTSSdata<-ifelse(mcols(amplicons)$strand=="+",
                   mcols(amplicons)$TSSdata-mcols(amplicons)$start+1,
                   mcols(amplicons)$end-mcols(amplicons)$TSSdata+1)

relTSSmost5p<-ifelse(mcols(amplicons)$strand=="+",
                     mcols(amplicons)$most5p-mcols(amplicons)$start+1,
                     mcols(amplicons)$end-mcols(amplicons)$most5p+1)

relTSSwb<-ifelse(mcols(amplicons)$strand=="+",
                 mcols(amplicons)$TSSwb-mcols(amplicons)$start+1,
                 mcols(amplicons)$end-mcols(amplicons)$TSSwb+1)

notMost5p<-relTSSwb<relTSSdata


########## read in motifs from motifDb and files
tata<-query(MotifDb,'TATA') # number 5 is similar
seqLogo(as.list(tata)[[5]])
tata<-tata[[5]]
#seqLogo(tata)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSdata,-30,rc=FALSE)
sum(mMatches$relPos<(-25) & mMatches$relPos>(-45),na.rm=TRUE)
#13 -30:-45
#21 -25:-45
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -30",col="grey")
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSwb,-30,rc=FALSE)
sum(mMatches_wb$relPos<(-30) & mMatches_wb$relPos>(-45),na.rm=TRUE)
#4 
hist(mMatches_wb$relPos,breaks=50, main="position of TATA closest to -30",col="grey",xlim=c(-250,300))
which(mMatches_wb$rePos<(-500))
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)


#sp1<-query(MotifDb,'SP1')
#seqLogo(sp1[[2]])
sp1<-readMemeChipMotif("./fromMemeChip/sp1.txt")
seqLogo(sp1)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,sp1,relTSSdata,-45,rc=TRUE)
sum(mMatches$relPos<(-40) & mMatches$relPos>(-50),na.rm=TRUE)
#7
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,sp1,relTSSwb,-45,rc=TRUE)
sum(mMatches_wb$relPos<(-40) & mMatches_wb$relPos>(-50),na.rm=TRUE)
#3 
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)


tblocks<-readMemeChipMotif("./fromMemeChip/tblockMotif.txt")
seqLogo(tblocks)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSdata,-50,rc=FALSE)
sum(mMatches$relPos<(20) & mMatches$relPos>(-200),na.rm=TRUE)
#6
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSwb,-50,rc=FALSE)
sum(mMatches_wb$relPos<(20) & mMatches_wb$relPos>(-200),na.rm=TRUE)
#29
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)


#inr<-query(MotifDb,'Inr') # not goot likeness
#from Jaspar http://jaspar.genereg.net/cgi-bin/jaspar_db.pl
inr<-rbind(c( 49,   0, 288,  26,  77,  67,  45,  50),
           c( 48, 303,   0,  81,  95, 118,  85,  96),
           c(69,   0,   0, 116,   0,  46,  73,  56),
           c(137,   0,  15,  80, 131,  72, 100, 101))
row.names(inr)=c("A","C","G","T")
colnames(inr)=c(1:8)
inr<-inr/colSums(inr)
seqLogo(inr)

#inr='TCAKTY'
#mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSdata,0,rc=FALSE,pattern=TRUE)
#sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#65

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSdata,0,rc=FALSE)
sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#69
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSwb,0,rc=FALSE)
sum(mMatches_wb$relPos<(-10) & mMatches_wb$relPos>(10),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)


motifFile<-"./fromMemeChip/motif_3_counts.txt" #in 26 promoters
bisl1a<-readMemeChipMotif(motifFile)
seqLogo(bisl1a)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1a,relTSSdata,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,bisl1a,relTSSwb,-50,rc=TRUE)
sum(mMatches_wb$relPos<(-150) & mMatches_wb$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)



motifFile<-"./fromMemeChip/motif_4_counts.txt" #in 9 promoters
bisl1d<-readMemeChipMotif(motifFile)
seqLogo(bisl1d)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1d,relTSSdata,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,bisl1d,relTSSwb,-50,rc=TRUE)
sum(mMatches_wb$relPos<(-150) & mMatches_wb$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)



motifFile<-"./fromMemeChip/motif_9_counts.txt" #in 7 promoters
bisl1b<-readMemeChipMotif(motifFile)
seqLogo(bisl1b)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1b,relTSSdata,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches_wb<-searchAmplicons4pwm(amplicons,ampMeta,bisl1b,relTSSwb,-50,rc=TRUE)
sum(mMatches_wb$relPos<(-150) & mMatches_wb$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches_wb)

#motifFile<-"./fromMemeChip/motif_10_counts.txt"
#seqLogo(readMemeChipMotif(motifFile))

saveRDS(amplicons,file="./amplicons_motifs.RDS")
write.csv(mcols(amplicons),file="./amplicons_motifs.csv",row.names=FALSE)


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

