# 2017-05-02
# modified to get modes of TSSs after seeing alignments with motifs are wrong
# 2016-11-17
# getTSS.R
# using output of convertTSSdata2WS235.R which lifted Chen and Kreus data over to WS235

setwd("~/Documents/MeisterLab/dSMF/TSS/scripts/")
#library(QuasR)
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

#resize KreusTSS1 to get a larger range for overlaps
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

# overlap with reduced Kreus TSS1
ol<-findOverlaps(KreusTSS1,ChenTSS)

KreusTSS_ol<-KreusTSS1[queryHits(ol)]  #3053 from 5400
ChenTSS_ol<-ChenTSS[subjectHits(ol)]   #3053 from 12545
#2781 TSSs
# if use Kreus width 50: 3053 TSSs

KreusTSS_ol1<-resize(KreusTSS_ol,width=1,fix="center")
KreusTSS_ol1<-shift(KreusTSS_ol1,shift=1)

hist(start(KreusTSS_ol1)-mcols(ChenTSS_ol)$modePosition,breaks=20,
     xlab="KruesTSSposition-ChenTSSposition",
     main="Distance between Kreus and Chen TSS")

#abline(v=c(-50,0),col="red")
abline(v=c(-50,50),col="red")

hist(apply(cbind(start(KreusTSS_ol1),mcols(ChenTSS_ol)$modePosition),1,sd),breaks=20,
     xlab="sd(KreusTSSposition,ChenTSSposition)",
     main="std dev between Kreus and Chen TSS")

#copy single bp TSS GRanges for Kreus data to a new object to preserve it for next two comparison
KreusTSS_ol<-KreusTSS_ol1

#take only ranges where the TSS is no further than 50 bp
KreusTSS_ol1<-KreusTSS_ol[(start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)>=-50 &
                            (start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)<=50,]

#2938 if use -50 to 50

#do this on chen data too to preserve metadata
ChenTSS_ol1<-ChenTSS_ol[(start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)>=-50 &
                            (start(KreusTSS_ol)-mcols(ChenTSS_ol)$modePosition)<=50,]


# do they refer to the same genes?
sum((mcols(KreusTSS_ol1)$WBGeneID)!=mcols(ChenTSS_ol1)$WbgeneName)
# 100 genes differ (109 differ if -50 to 50)
i<-which((mcols(KreusTSS_ol1)$WBGeneID)!=mcols(ChenTSS_ol1)$WbgeneName)
head(sort(mcols(KreusTSS_ol1)$WBGeneID[i]))
head(sort(mcols(ChenTSS_ol1)$WbgeneName[i]))
# remove rows where genes differ.
KreusTSS_ol1<-KreusTSS_ol1[-i,]
ChenTSS_ol1<-ChenTSS_ol1[-i,]
#2476 ranges if -50 to 0
# 2829 ranges if -50 to 50

table(seqnames(ChenTSS_ol1))
# chrI  chrII chrIII  chrIV   chrV   chrX
# 419    411    377    384    521    364
# for -50 to 0

#chrI  chrII chrIII  chrIV   chrV   chrX
#471    478    433    442    589    416
# for -50 to 50

#any duplicate genes with more than one TSS?
length(unique(mcols(ChenTSS_ol1)$WbgeneName,mcols(KreusTSS_ol1)$WBGeneID))
# 2476 ranges. no!
# 2829 ranges for -50 to 50. no!

################################3
#get Mode of TSS
#downloaded TSS wig files for WS250 from:
#ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/annotation/TSS/
source("../../usefulFunctions.R")
source("../../findMotifs/findMotifs_functions.R")
library(Biostrings)
library(Gviz)
#library(grid)
#library(gridExtra)

#create seqinfo object for conversion of wig to bigwig
SI<-seqinfo(BSgenome.Celegans.UCSC.ce11)
SI<-ucsc2wb(SI)
genome(SI)<-"WS250"

#convert wig to bigwig binary format (do only once)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Chen.Forward.wig",SI)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Chen.Reverse.wig",SI)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Kruesi.Forward.wig",SI)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Kruesi.Reverse.wig",SI,clip=TRUE)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Saito.Forward.wig",SI)
#wigToBigWig("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Saito.Reverse.wig",SI)

#copy gene names to keep track of id
names(ChenTSS_ol1)<-mcols(ChenTSS_ol1)$WbgeneName
names(KreusTSS_ol1)<-mcols(KreusTSS_ol1)$WBGeneID
#find overlap between chen range and kreus artifical 100 bp range and then reduce to create
#a single probably TSS range per gene in which to look at wiggle data
maxoltemp<-c(granges(ChenTSS_ol1),granges(resize(KreusTSS_ol1,width=101,fix="center")))
maxol<-reduce(maxoltemp)
#findoverlaps to copy accross gene names
tempol<-findOverlaps(maxol,maxoltemp)
names(maxol)<-unique(names(maxoltemp[subjectHits(tempol)]))
#chnage seqnames to WB genome format
seqlevels(maxol)<-gsub("chr","",seqlevels(maxol))

#import WS250 wiggle (bigwig) files of TSS data from three studies
chenFbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Chen.Forward.bw",format="BigWig",
                which=maxol)
chenRbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Chen.Reverse.bw",format="BigWig",
                which=maxol)
kreusFbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Kruesi.Forward.bw",format="BigWig",
                which=maxol)
kreusRbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Kruesi.Reverse.bw",format="BigWig",
                which=maxol)
saitoFbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Saito.Forward.bw",format="BigWig",
                which=maxol)
saitoRbw<-import("../dataSets/c_elegans.PRJNA13758.WS250.TSS.wig/TSS_Saito.Reverse.bw",format="BigWig",
                which=maxol)

# create GR with a sum of all three data sets (forward strand)
uF<-unique(c(granges(chenFbw),granges(kreusFbw),granges(saitoFbw)))
mcols(uF)[match(chenFbw,uF),"chenFbw"]<-mcols(chenFbw)$score
mcols(uF)[match(kreusFbw,uF),"kreusFbw"]<-mcols(kreusFbw)$score
mcols(uF)[match(saitoFbw,uF),"saitoFbw"]<-mcols(saitoFbw)$score
mcols(uF)$allF<-apply(mcols(uF),1,sum,na.rm=TRUE)

# create GR with a sum of all three data sets (reverse strand)
uR<-unique(c(granges(chenRbw),granges(kreusRbw),granges(saitoRbw)))
mcols(uR)[match(chenRbw,uR),"chenRbw"]<-mcols(chenRbw)$score
mcols(uR)[match(kreusRbw,uR),"kreusRbw"]<-mcols(kreusRbw)$score
mcols(uR)[match(saitoRbw,uR),"saitoRbw"]<-mcols(saitoRbw)$score
mcols(uR)$allR<-apply(mcols(uR),1,sum,na.rm=TRUE)

#modes<-apply(mcols(uR[subjectHits(ol)]),2,which.max)

# construct PWM for Inr motif
inr<-rbind(c( 49,   0, 288,  26,  77,  67,  45,  50),
           c( 48, 303,   0,  81,  95, 118,  85,  96),
           c(69,   0,   0, 116,   0,  46,  73,  56),
           c(137,   0,  15,  80, 131,  72, 100, 101))
row.names(inr)=c("A","C","G","T")
colnames(inr)=c(1:8)
inr<-inr/colSums(inr)
# get DNA sequences for putative TSS regions
maxolseqs<-getSeq(BSgenome.Celegans.UCSC.ce11,wb2ucsc(maxol))

# ############ plotting TSSs with inr motifs to help choose between most5p and maxTSS
# #prepare universal tracks
# gtrack<-GenomeAxisTrack()
# options(ucscChromosomeNames=FALSE)
# dtrack1<-DataTrack(range=chenFbw,type="histogram",name="Chen F")
# dtrack2<-DataTrack(range=chenRbw,type="histogram",name="Chen R")
# dtrack3<-DataTrack(range=kreusFbw,type="histogram",name="Kruesi F")
# dtrack4<-DataTrack(range=kreusRbw,type="histogram",name="Kruesi R")
# dtrack5<-DataTrack(range=saitoFbw,type="histogram",name="Saito F")
# dtrack6<-DataTrack(range=saitoRbw,type="histogram",name="Saito R")
#
#
# txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
# GeneRegionTrack(txdb)
# txTr<-GeneRegionTrack(txdb,showFeatureId=TRUE,name="txdb")
# #genes<-genes(txdb)
# #plotTracks(txTr)
# #feature(txTr)
#
# getMotifAnnotationTrack<-function(motif,DNAseq,gr) {
#    matches<-matchPWM(motif,DNAseq)
#    if(length(matches)==0) {
#       atrack<-AnnotationTrack()
#       return(atrack)
#    }
#    motifName<-getOriginalVarName(motif)
#    names(matches)<-paste0(motifName,"_",seq_along(matches))
#    regionMeta<-DataFrame(start=start(gr), end=end(gr), strand=strand(gr), chr=seqnames(gr))
#    genomeCoord_gr<-matches2GR(regionMeta,matches)
#    #motifName<-getOriginalVarName(motif)
#    atrack<-AnnotationTrack(genomeCoord_gr,shape="box",fill="light blue",name="Inr")
#    #print(genomeCoord_gr)
#    return(atrack)
# }
#
# pdf("./allPromoterTSS&inr_271_.pdf",width=6,height=6)
# #par(mfrow=c(2,2))
# #layout(matrix(c(1,2,3,4),nrow=2,byrow=TRUE))
# # run through TSS regions and create plots for each region
# for (i in 272:length(maxol)){
#    #find inr motifs in this region
#    atrack_inr<-getMotifAnnotationTrack(inr,maxolseqs[[i]],maxol[i])
#    #atrack_tata<-getMotifAnnotationTrack(tata,maxolseqs[[i]],maxol[i])
#
#    #check which strand the gene is on for different TSS counts
#    if(as.character(strand(maxol[i]))=="+") {
#       ol<-subsetByOverlaps(uF,maxol[i])
#       modes<-unlist(apply(mcols(ol),2,which.max))
#       atrack_TSS<-AnnotationTrack(ol[modes,],id=as.character(names(modes)),name="TSS",
#                                   showId=TRUE,fill="black")
#       plotTracks(list(gtrack,dtrack1,dtrack3,dtrack5,atrack_TSS,atrack_inr,txTr),
#                  chromosome=seqnames(maxol[i]),from=start(maxol[i])-5,to=end(maxol[i])+5)
#    } else {
#       ol<-subsetByOverlaps(uR,maxol[i])
#       modes<-unlist(apply(mcols(ol),2,which.max))
#       atrack_TSS<-AnnotationTrack(ol[modes,],id=as.character(names(modes)),name="TSS",
#                                   showId=TRUE,fill="black")
#       plotTracks(list(gtrack,dtrack2,dtrack4,dtrack6,atrack_TSS,atrack_inr,txTr),
#                  chromosome=seqnames(maxol[i]),from=start(maxol[i])-5,to=end(maxol[i])+5)
#    }
# }
# dev.off()
#



###### Store all mode for each TSS, as well as mode of sum of datasets and the most 5p ###########
chenTSS<-rep(NA,length(maxol))
kreusTSS<-rep(NA,length(maxol))
saitoTSS<-rep(NA,length(maxol))
maxTSS<-rep(NA,length(maxol))
most5p<-rep(NA,length(maxol))

for (i in 1:length(maxol)){
   #check which strand the gene is on for different TSS counts
   if(as.character(strand(maxol[i]))=="+") {
      ol<-subsetByOverlaps(uF,maxol[i])
      modes<-unlist(apply(mcols(ol),2,which.max))
      #need to check for NAs because some datasets have no data in range
      chenTSS[i]<-getStartTSS(modes,"chenFbw")
      kreusTSS[i]<-getStartTSS(modes,"kreusFbw")
      saitoTSS[i]<-getStartTSS(modes,"saitoFbw")
      maxTSS[i]<-getStartTSS(modes,"allF")
      most5p[i]<-max(start(ol[modes],na.rm=TRUE))
   } else {
      ol<-subsetByOverlaps(uR,maxol[i])
      modes<-unlist(apply(mcols(ol),2,which.max))
      #need to check for NAs because some datasets have no data in range
      chenTSS[i]<-getStartTSS(modes,"chenRbw")
      kreusTSS[i]<-getStartTSS(modes,"kreusRbw")
      saitoTSS[i]<-getStartTSS(modes,"saitoRbw")
      maxTSS[i]<-getStartTSS(modes,"allR")
      most5p[i]<-max(start(ol[modes],na.rm=TRUE))
   }
}


###################################
mcols(maxol)<-DataFrame(chenTSS=chenTSS,kruesiTSS=kreusTSS,saitoTSS=saitoTSS,maxTSS=maxTSS, most5p=most5p)

saveRDS(maxol,"ChenKreusSaitoTSS_2827.RDS")
export.bed(maxol,"ChenKreusSaitoTSS_2827.bed")
WBGeneID<-names(maxol)
chr<-seqnames(maxol)
start<-start(maxol)
end<-end(maxol)
strand<-strand(maxol)
df<-data.frame(WBGeneID=WBGeneID,chr=chr,start=start,end=end,strand=strand,mcols(maxol))
write.table(df,"ChenKreusSaitoTSS_2827.txt",row.names=FALSE,col.names=TRUE)
sum(rowSums(abs(df[,6:10]-df$maxTSS))<5,na.rm=TRUE)
#1053
sum(rowSums(abs(df[,6:10]-df$maxTSS))<1,na.rm=TRUE)
#872
highConfidence<-na.omit(df[rowSums(abs(df[,6:10]-df$maxTSS))<1,])
write.table(highConfidence,"ChenKreusSaitoTSS_highConf_872.txt",row.names=FALSE,col.names=FALSE)
highConfidence_GR<-with(highConfidence, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand))
names(highConfidence_GR)<-highConfidence$WBGeneID
mcols(highConfidence_GR)<-highConfidence[6:10]
saveRDS(highConfidence_GR,"ChenKreusSaitoTSS_highConf_872.RDS")


#### write text files with sequences surrounding TSSs to use for motif searches
library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)

# get sequences for all not-too-variable TSSs
maxTSS<-fastaFromGR(granges(maxol),mcols(maxol)$maxTSS,150,150,"./allTSSmax_-150_150.fa")
most5pTSS<-fastaFromGR(granges(maxol),mcols(maxol)$most5p,150,150,"./allTSSmost5p_-150_150.fa")

sum(start(maxTSS)==start(most5pTSS)) #1805 are the same
sum(abs(start(maxTSS)-start(most5pTSS))<=10) #2227 are within 10bp of eachother
#2428 are within 20 bp of eachother
#2788 are within 50 bp of eachother
#2826 are within 100 bp of eachother

# get sequences for all TSSs that are identical in all datasets
perfectTSS<-fastaFromGR(granges(highConfidence_GR),mcols(highConfidence_GR)$maxTSS,150,150,"./perfectTSS_-150_150.fa")

# ################### no longer used ###########################
# ### get upstream most TSS in between Chen and Kreus data:
# #note, because in previous script i had not re-shrunk the KreusTSSs I selected many
# #Kreus TSS that were placed artifically upstream.
# # note, only 2037 of my wig file modes match the Chen et al. "modePosition" call in the paper
#
# # get index for genes on each strand
# fwStrand<-which(strand(ChenTSS_ol1)=="+")
# rcStrand<-which(strand(ChenTSS_ol1)=="-")
#
# # get most upstream TSS (depending on strand)
# mostUpstreamTSS<-rep(0,length(ChenTSS_ol1))
# #create DF with both TSSs:
# bothTSS<-data.frame("ChenTSS"=mcols(ChenTSS_ol1)$modePosition,
#                     "KreusTSS"=start(KreusTSS_ol1))
# mostUpstreamTSS[fwStrand]<-apply(bothTSS,1,min,na.rm=TRUE)[fwStrand]
# mostUpstreamTSS[rcStrand]<-apply(bothTSS,1,max,na.rm=TRUE)[rcStrand]
#
# ConsensusTSS<-GRanges(seqnames=seqnames(ChenTSS_ol1),
#                       ranges=IRanges(start=mostUpstreamTSS,
#                                      end=mostUpstreamTSS),
#                       strand=strand(ChenTSS_ol1))
#
# mcols(ConsensusTSS)<-cbind(mcols(KreusTSS_ol1),
#                            mcols(ChenTSS_ol1)[,c("modePosition",
#                                                  "assignmentType",
#                                                  "assignedGeneName")])
#
# #2829 confirmed similar TSSs
# saveRDS(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.Rds")
# export(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.gff",format="gff",version="3")
# export(ConsensusTSS,"confirmedTSS_Chen&Kreus_WS235.bed",format="bed")
