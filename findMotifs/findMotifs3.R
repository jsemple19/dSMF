#2017-05-17
#script identify the presence of some positional motifs in amplicons

setwd("~/Documents/MeisterLab/dSMF/findMotifs")
library(MotifDb)
library(seqLogo)
library(GenomicRanges)
library(rtracklayer)
#library(MotIV)

source("./findMotifs_functions.R")

# read in amplicons
amplicons<-readRDS("../PromoterPrimerDesign/scripts/ampliconsSeqs.RDS")

# get vectors of TSSs with different definitions
relTSSmax<-ifelse(mcols(amplicons)$strand=="+",
                  mcols(amplicons)$maxTSS-mcols(amplicons)$start+1,
                  mcols(amplicons)$end-mcols(amplicons)$maxTSS+1)

relTSSmost5p<-ifelse(mcols(amplicons)$strand=="+",
                     mcols(amplicons)$most5p-mcols(amplicons)$start+1,
                     mcols(amplicons)$end-mcols(amplicons)$most5p+1)

relTSSwb<-ifelse(mcols(amplicons)$strand=="+",
                 mcols(amplicons)$TSSwb-mcols(amplicons)$start+1,
                 mcols(amplicons)$end-mcols(amplicons)$TSSwb+1)

notMost5p<-relTSSwb<relTSSmax
#34
sum(relTSSwb<relTSSmost5p)
#30

########## read in motifs from motifDb and files

######################## inr ###################
#from Jaspar http://jaspar.genereg.net/cgi-bin/jaspar_db.pl
inr<-rbind(c( 49,   0, 288,  26,  77,  67,  45,  50),
           c( 48, 303,   0,  81,  95, 118,  85,  96),
           c(69,   0,   0, 116,   0,  46,  73,  56),
           c(137,   0,  15,  80, 131,  72, 100, 101))
row.names(inr)=c("A","C","G","T")
colnames(inr)=c(1:8)
inr<-inr/colSums(inr)
seqLogo(inr)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmax,0,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of Inr closest to +1",col="grey")
sum(mMatches$relPos>(-5) & mMatches$relPos<(5),na.rm=TRUE)
#61
mMatches[which(mMatches$relPos<(-5) | mMatches$relPos>(5)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmost5p,0,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of Inr closest to +1",col="grey")
sum(mMatches$relPos>(-5) & mMatches$relPos<(5),na.rm=TRUE)
#52
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


############################ TATA ####################################
tata<-readMemeChipMotif("./motifs/tata_tbp.txt") #MA0108.2 from meme-combined
seqLogo(tata)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSmax,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#24 -20:-35
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey")
mMatches[which(mMatches$relPos>(-20) | mMatches$relPos<(-35)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSmost5p,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#18 -20:-35
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey")
mMatches[which(mMatches$relPos>(-20) | mMatches$relPos<(-35)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


########################## cg rich ##################################
cgrich<-readMemeChipMotif("./motifs/cgrich_meme-2.txt")
seqLogo(cgrich)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmax,-45,rc=TRUE)
sum(mMatches$relPos<(-40) & mMatches$relPos>(-60),na.rm=TRUE)
#15 (-40)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -45",col="grey",xlim=c(-100,0))
mMatches[which(mMatches$relPos>(-40) | mMatches$relPos<(-60)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmost5p,-45,rc=TRUE)
sum(mMatches$relPos<(-40) & mMatches$relPos>(-60),na.rm=TRUE)
#16 (-40)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -45",col="grey",xlim=c(-100,0))
mMatches[which(mMatches$relPos>(-40) | mMatches$relPos<(-60)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)




saveRDS(amplicons,file="./amplicons_motifs.RDS")
write.csv(mcols(amplicons),file="./amplicons_motifs.csv",row.names=FALSE)
