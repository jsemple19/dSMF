# 2016-11-17
# getTSS.R
# using output of convertTSSdata2WS235.R which lifted Chen and Kreus data over to WS235

setwd("/SharedDocuments/MeisterLab/dSMF/TSS/scripts/")
library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

source("./getTSS_functions.R")


#######################

# finding TSS that are least variable
KreusTSS<-readRDS("../dataSets/Kreus2013/KreusTSS_WS235.Rds")
TSScols<-c("DCCmutantEmbryo_TSS","WildTypeEmbryo_TSS","WildTypeStarvedL1_TSS","WildTypeL3_TSS","SameTSS")
# convert to numeric the TSSs
TSSpos<-as.data.frame(mcols(KreusTSS)[,TSScols][,1:4])
#calculate the SD of the TSS position
TSS_sd<-apply(TSSpos,1,sd,na.rm=TRUE)
# if only in one dataset, output is NA. If identical in two, output is 0.0000. Looking at histograme - maybe 
# use cutoff of SD=20
hist(TSS_sd,breaks=100,main="Std dev of TSS position (Kreus 2013)",xlab="SD of TSS")
hist(TSS_sd,breaks=10000,xlim=c(0,30),main="Std dev of TSS position (Kreus 2013)",xlab="SD of TSS")

# genes with TSS in at least 2 conditions which is in the same place
TSS_sameIn2<-TSS_sd==0
#1907 genes

# genes with TSS in at least 2 conditions which is close
TSS_closeIn2<-TSS_sd<=20
# 3047 genes

# genes with TSS identified in only 1 condition
TSS_In1<-is.na(TSS_sd)
# 2353 genes

# get index for genes on each strand
fwStrand<-which(strand(KreusTSS)=="+")
rcStrand<-which(strand(KreusTSS)=="-")

# get most upstream TSS (depending on strand)
TSSpos<-cbind(TSSpos,"chosenTSS"=0)
TSSpos[fwStrand,"chosenTSS"]<-apply(TSSpos[,1:4],1,min,na.rm=TRUE)[fwStrand]
TSSpos[rcStrand,"chosenTSS"]<-apply(TSSpos[,1:4],1,max,na.rm=TRUE)[rcStrand]

#WBgeneStarts<-ifelse(strand(KreusTSS)=="+",start(KreusTSS)+1,end(KreusTSS)-1)

#TSSpos<-cbind(TSSpos,"wbStart"=0)
#wb5p<-promoters(KreusTSS,upstream=0,downstream=1)

# create a new GRanges object for TSS
constantTSS<-GRanges(seqnames(KreusTSS),ranges=IRanges(start=TSSpos[,"chosenTSS"],
                    end=TSSpos[,"chosenTSS"]), strand=strand(KreusTSS))
# add gene names and TSS in different conditions as metadata columns
mcols(constantTSS)<-cbind(mcols(KreusTSS)[,c("WBGeneID","GeneName")],data.frame(TSSpos[,1:4]))

# subset to only TSSs with little positional variation (or only identified in a single condition)
KreusTSS1a<-constantTSS[(TSS_sameIn2 | TSS_closeIn2 | TSS_In1),]
#5400 (down from 6353)

KreusTSS1<-resize(KreusTSS1a,width=50,fix="center")

#GR5p<-promoters(KreusTSS,upstream=20000,downstream=100)
#wb5p<-promoters(KreusTSS,upstream=0,downstream=1)

#mcols(GR5p)<-cbind(mcols(KreusTSS)[,c("WBGeneID","GeneName")],data.frame(TSSpos[,1:4]))

#GR5p<-GR5p[(TSS_sameIn2 | TSS_closeIn2 | TSS_In1),]

#reduce(GR5p,constantTSS1)


############# import Chen 2013 data ###############
ChenTSS<-readRDS("../dataSets/Chen2013/ChenTSS_WS235.Rds")

#take only wormbase_tss and raft_to_wormbase_tss assignments
ChenTSS<-ChenTSS[mcols(ChenTSS)$assignmentType %in% c("wormbase_tss","raft_to_wormbase_tss"),]

# remove rows with NA in wormbase gene id (no physical location in wormbase... "WbgeneName")
ChenTSS<-ChenTSS[-which(is.na(mcols(ChenTSS)$WbgeneName)),]
# 75 TSSs removed


hist(width(ChenTSS),breaks=100)
quantile(width(ChenTSS),c(0.90,0.95,0.99))
# 90%    95%    99% 
# 107.00 135.00 196.63 
sum(width(ChenTSS)<=100)/length(width(ChenTSS))
# 0.8830272
sum(width(ChenTSS)<=150)/length(width(ChenTSS))
# 0.9666353
sum(width(ChenTSS)<=200)/length(width(ChenTSS))
# 0.9907364

#118 tSS removed if >200
ChenTSS<-ChenTSS[width(ChenTSS)<=200,]

ol<-findOverlaps(KreusTSS1,ChenTSS)

KreusTSS_ol<-KreusTSS1[queryHits(ol)]
ChenTSS_ol<-ChenTSS[subjectHits(ol)]
#2781 TSSs
# if use Kreus width 50: 3053 TSSs

hist(start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition,breaks=20,
     xlab="KruesTSSposition-ChenTSSposition",
     main="Distance between Kreus and Chen TSS")

#abline(v=c(-50,0),col="red")
abline(v=c(-75,25),col="red")

hist(apply(cbind(start(KreusTSS_ol),mcols(ChenTSS_ol)$modePosition),1,sd),breaks=20,
     xlab="sd(KreusTSSposition,ChenTSSposition)",
     main="std dev between Kreus and Chen TSS")

KreusTSS_ol1<-KreusTSS_ol[(start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)>=-75 & 
                            (start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)<=25,]

#2576 ranges if -50 to 0
#2938 ranges if -75 to 25

ChenTSS_ol1<-ChenTSS_ol[(start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)>=-75 & 
                            (start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)<=25,]

# do they refer to the same genes?
sum((mcols(KreusTSS_ol1)$WBGeneID)!=mcols(ChenTSS_ol1)$WbgeneName)
# 100 genes differ (109 differ if -75 to 25)
i<-which((mcols(KreusTSS_ol1)$WBGeneID)!=mcols(ChenTSS_ol1)$WbgeneName)
tail(mcols(KreusTSS_ol1)$WBGeneID[i])
tail(mcols(ChenTSS_ol1)$WbgeneName[i])
# remove rows where genes differ.
KreusTSS_ol1<-KreusTSS_ol1[-i,]
ChenTSS_ol1<-ChenTSS_ol1[-i,]
#2476 ranges if -50 to 0
# 2829 ranges if -75 to 25

table(seqnames(ChenTSS_ol1))
# chrI  chrII chrIII  chrIV   chrV   chrX 
# 419    411    377    384    521    364 
# for -50 to 0

#chrI  chrII chrIII  chrIV   chrV   chrX 
#471    478    433    442    589    416
# for -75 to 0

#any duplicate genes with more than one TSS?
length(unique(mcols(ChenTSS_ol1)$WbgeneName,mcols(KreusTSS_ol1)$WBGeneID))
# 2476 ranges. no!
# 2829 ranges for -75 to 25. no!

### get upstream most TSS in between Chen and Kreus data:

# get index for genes on each strand
fwStrand<-which(strand(ChenTSS_ol1)=="+")
rcStrand<-which(strand(ChenTSS_ol1)=="-")

# get most upstream TSS (depending on strand)
mostUpstreamTSS<-rep(0,length(ChenTSS_ol1))
#create DF with both TSSs:
bothTSS<-data.frame("ChenTSS"=mcols(ChenTSS_ol1)$modePosition,
                    "KreusTSS"=start(KreusTSS_ol1))
mostUpstreamTSS[fwStrand]<-apply(bothTSS,1,min,na.rm=TRUE)[fwStrand]
mostUpstreamTSS[rcStrand]<-apply(bothTSS,1,max,na.rm=TRUE)[rcStrand]

ConsensusTSS<-GRanges(seqnames=seqnames(ChenTSS_ol1), 
                      ranges=IRanges(start=mostUpstreamTSS,
                                     end=mostUpstreamTSS),
                      strand=strand(ChenTSS_ol1))

mcols(ConsensusTSS)<-cbind(mcols(KreusTSS_ol1),
                           mcols(ChenTSS_ol1)[,c("modePosition",
                                                 "assignmentType",
                                                 "assignedGeneName")])

saveRDS(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.Rds")
export(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.gff",format="gff",version="3")
export(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.bed",format="bed")
