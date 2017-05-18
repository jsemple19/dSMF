#script exploring the presence of positonal motifs in the aplicons. 

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
#mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmax,0,rc=FALSE,pattern=TRUE)
#sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#65

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmax,0,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of Inr closest to +1",col="grey")
sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#72
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmost5p,0,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of Inr closest to +1",col="grey")
sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#64
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSwb,0,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of Inr closest to +1",col="grey")
sum(mMatches$relPos>(-10) & mMatches$relPos<(10),na.rm=TRUE)
#57
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


inr<-readMemeChipMotif("./motifs/inr_hbp1.txt")
seqLogo(inr)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmax,0,rc=TRUE)
sum(mMatches$relPos>(-10) & mMatches$relPos>(10),na.rm=TRUE)
#22
hist(mMatches$relPos,breaks=50, main="position of inr closest to +1",col="grey",xlim=c(-100,0))

inr<-readMemeChipMotif("./motifs/inr_tec-1.txt")
seqLogo(inr)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,inr,relTSSmax,0,rc=TRUE)
sum(mMatches$relPos>(-10) & mMatches$relPos>(10),na.rm=TRUE)
#35
hist(mMatches$relPos,breaks=50, main="position of inr closest to +1",col="grey",xlim=c(-100,0))

############################ TATA ####################################
########## from motifDb
tata<-query(MotifDb,'TATA') # number 5 is similar
seqLogo(as.list(tata)[[5]])
tata<-tata[[5]]
#seqLogo(tata)
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


mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSwb,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#16 -20:-35 
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey",xlim=c(-250,300))
mMatches[which(mMatches$relPos>(-20) | mMatches$relPos<(-35)),]<-NA
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

###### meme-3
tata<-readMemeChipMotif("./motifs/tata_meme-3.txt")
seqLogo(tata)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSmax,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#16 -20:-35
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey")

tata<-readMemeChipMotif("./motifs/tata_tbp.txt")
seqLogo(tata)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSmax,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#24 -20:-35
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey")

tata<-readMemeChipMotif("./motifs/tata_dreme-1.txt")
seqLogo(tata)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tata,relTSSmax,-28,rc=FALSE)
sum(mMatches$relPos<(-20) & mMatches$relPos>(-35),na.rm=TRUE)
#18 -20:-35
hist(mMatches$relPos,breaks=50, main="position of TATA closest to -28",col="grey")

######################## SP1 or CG rich -50 region ########################
#sp1<-query(MotifDb,'SP1')
#seqLogo(sp1[[2]])
sp1<-readMemeChipMotif("./motifs/sp1.txt")
seqLogo(sp1)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,sp1,relTSSmax,-45,rc=TRUE)
sum(mMatches$relPos<(-35) & mMatches$relPos>(-60),na.rm=TRUE)
#9
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,sp1,relTSSwb,-45,rc=TRUE)
sum(mMatches$relPos<(-35) & mMatches$relPos>(-60),na.rm=TRUE)
#4
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

cgrich<-readMemeChipMotif("./motifs/cgrich_ascl2.txt")
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmax,-60,rc=TRUE)
sum(mMatches$relPos<(-40) & mMatches$relPos>(-80),na.rm=TRUE)
#15 (-70, periodic?)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -60",col="grey",xlim=c(-100,0))

cgrich<-readMemeChipMotif("./motifs/cgrich_RSC30.txt")
seqLogo(cgrich)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmax,-60,rc=TRUE)
sum(mMatches$relPos<(-40) & mMatches$relPos>(-80),na.rm=TRUE)
#15 (a abit more upstream at -40..-70)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -60",col="grey",xlim=c(-100,0))


cgrich<-readMemeChipMotif("./motifs/cgrich_meme-2.txt")
seqLogo(cgrich)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmax,-45,rc=TRUE)
sum(mMatches$relPos<(-35) & mMatches$relPos>(-60),na.rm=TRUE)
#15 (-40)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -45",col="grey",xlim=c(-100,0))

cgrich<-readMemeChipMotif("./motifs/cgrich_periodic_sut1.txt")
seqLogo(cgrich)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,cgrich,relTSSmax,-45,rc=TRUE)
sum(mMatches$relPos<(-35) & mMatches$relPos>(-60),na.rm=TRUE)
#18 (a bit more downstream closer to 40)
hist(mMatches$relPos,breaks=50, main="position of cg rich closest to -45",col="grey",xlim=c(-100,0))

########################## tblocks  #########################################
###### not a good method to look at periodic motifs!!!
tblocks<-readMemeChipMotif("./motifs/tblockMotif.txt")
seqLogo(tblocks)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSmax,30,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of tblocks closest to 30",col="grey",xlim=c(-100,200))
sum(mMatches$relPos<(200) & mMatches$relPos>(-200),na.rm=TRUE)
#42
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSmost5p,30,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of tblocks closest to 30",col="grey",xlim=c(-100,200))
sum(mMatches$relPos<(200) & mMatches$relPos>(-200),na.rm=TRUE)
#42
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSwb,30,rc=FALSE)
hist(mMatches$relPos,breaks=50, main="position of tblocks closest to 30",col="grey",xlim=c(-100,200))
sum(mMatches$relPos<(200) & mMatches$relPos>(-200),na.rm=TRUE)
#28
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


tblocks<-readMemeChipMotif("./motifs/tblocks_dreme-2.txt")
seqLogo(tblocks)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSmax,50,rc=TRUE)
sum(mMatches$relPos>(-50) & mMatches$relPos<(100),na.rm=TRUE)
#94
hist(mMatches$relPos,breaks=50, main="position of tblocks closest to 50",col="grey")

tblocks<-readMemeChipMotif("./motifs/tblocks_meme-1.txt")
seqLogo(tblocks)
mMatches<-searchAmplicons4pwm(amplicons,ampMeta,tblocks,relTSSmax,50,rc=TRUE)
sum(mMatches$relPos>(-50) & mMatches$relPos<(100),na.rm=TRUE)
#45
hist(mMatches$relPos,breaks=50, main="position of tblocks closest to 50",col="grey")
#######################


motifFile<-"./motifs/motif_3_counts.txt" #in 26 promoters
bisl1a<-readMemeChipMotif(motifFile)
seqLogo(bisl1a)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1a,relTSSmax,-50,rc=TRUE)
sum(mMatches$relPos>(-150) & mMatches$relPos<(50),na.rm=TRUE)
hist(mMatches$relPos,breaks=50, main="position of motif3 closest to -50",col="grey")
#3
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)


motifFile<-"./motifs/motif_4_counts.txt" #in 9 promoters
bisl1d<-readMemeChipMotif(motifFile)
seqLogo(bisl1d)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1d,relTSSmax,-50,rc=TRUE)
hist(mMatches$relPos,breaks=50, main="position of motif3 closest to -50",col="grey")

sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1d,relTSSwb,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)



motifFile<-"./motifs/motif_9_counts.txt" #in 7 promoters
bisl1b<-readMemeChipMotif(motifFile)
seqLogo(bisl1b)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1b,relTSSmax,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

mMatches<-searchAmplicons4pwm(amplicons,ampMeta,bisl1b,relTSSwb,-50,rc=TRUE)
sum(mMatches$relPos<(-150) & mMatches$relPos>(50),na.rm=TRUE)
#0
mcols(amplicons)<-cbind(mcols(amplicons),mMatches)

#motifFile<-"./motifs/motif_10_counts.txt"
#seqLogo(readMemeChipMotif(motifFile))

saveRDS(amplicons,file="./amplicons_motifs.RDS")
write.csv(mcols(amplicons),file="./amplicons_motifs.csv",row.names=FALSE)


