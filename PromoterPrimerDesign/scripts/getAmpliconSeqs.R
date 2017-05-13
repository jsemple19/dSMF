setwd("~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts")
library("BSgenome.Celegans.UCSC.ce11")

amplicons<-readRDS("GRangesOfAmplicons.RDS")
names(mcols(amplicons))<-gsub("mcols.","",names(mcols(amplicons)))

genome<-BSgenome.Celegans.UCSC.ce11
#whichGenes<-c("WBGene00015955","WBGene00007811")
#i<-match(whichGenes,mcols(amplicons)$Gene_WB_ID)
amplicons.strand<-amplicons
strand(amplicons.strand)<-ifelse(mcols(amplicons)$orientation=="fw","+","-")

mySeqs<-getSeq(genome,amplicons.strand)#[i,])
#names(mySeqs)<-whichGenes

####################### to find digest fragment sizes ###########3

names(mySeqs)<-mcols(amplicons.strand)$Gene_WB_ID
BstUI<-unlist(vmatchPattern('CGCG',mySeqs))
BstUI.genes<-unique(names(BstUI))
length(BstUI.genes)
#40

HhaI<-unlist(vmatchPattern('GCGC',mySeqs))
HhaI.genes<-unique(names(HhaI))
length(HhaI.genes)
#51

length(union(BstUI.genes,HhaI.genes))
#59

TaqI<-unlist(vmatchPattern('TCGA',mySeqs))
TaqI.genes<-unique(names(TaqI))
length(TaqI.genes)
#76 (4000 or 20000 units)

MspI<-unlist(vmatchPattern('CCGG',mySeqs))
MspI.genes<-unique(names(MspI))
length(MspI.genes)
#43

HypCH4V<-unlist(vmatchPattern('TGCA',mySeqs))
HypCH4V.genes<-unique(names(HypCH4V))
length(HypCH4V.genes)
#80
#few units supplied (100 or 500)


AluI<-unlist(vmatchPattern('AGCT',mySeqs))
AluI.genes<-unique(names(AluI))
length(AluI.genes)
#59 (1000 or 5000 units)

length(union(TaqI.genes,AluI.genes))
#87

frags<-vector("list",96)

fragsFromDigest<-function(amplicons,restSites) {
   frags<-list()
   for (gene in 1:length(amplicons)){
      geneName<-as.character(mcols(amplicons)$Gene_WB_ID[gene])
      if(geneName %in% names(restSites)) {
         i<-which(names(restSites)==geneName)
         sites<-c(0,width(amplicons[gene,]),start(restSites[i]),end(restSites[i]))
         sites<-sort(sites)
         frags[[geneName]]<-paste(as.character(sort(diff(sites)[diff(sites)>5])), sep=",")
      } else {
         frags[[geneName]]<-width(amplicons[gene])
      }
   }
   return(frags)
}

frags<-fragsFromDigest(amplicons,TaqI)

frags<-unlist(lapply(frags,paste,collapse=","))

mcols(amplicons)$rowNum=rep(c("A","B","C","D","E","F","G","H"),each=12)
mcols(amplicons)$colNum=rep(1:12,8)
mcols(amplicons)$wellNum=paste0(mcols(amplicons)$rowNum,mcols(amplicons)$colNum)

plate<-data.frame(WBID=names(frags),frags,mcols(amplicons))
names(plate)

plate[plate$colNum==1,c("wellNum","frags")]
plate[plate$colNum==2,c("wellNum","frags")]

plate[plate$colNum==5,c("wellNum","frags")]
plate[plate$colNum==6,c("wellNum","frags")]

plate[plate$colNum==7,c("wellNum","frags")]
plate[plate$colNum==8,c("wellNum","frags")]


############## funcitons modified from Arnaud's functions in useful_functionsV1.r
#converts all C to T in any context
bisConv=function(in.seq){
   C.pos=vmatchPattern('C', in.seq)
   do.call(c,
           lapply(seq(length(in.seq)),function(i){
              Cind=start(C.pos[[i]])
              convind=Cind
              DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
           })
   )
}

#converts all C to T in  except for CG context
bisConvCG=function(in.seq){
   C.pos=vmatchPattern('C', in.seq)
   CG.pos=vmatchPattern('CG', in.seq)
   do.call(c,
           lapply(seq(length(in.seq)),function(i){
              Cind=start(C.pos[[i]])
              CGind=start(CG.pos[[i]])
              convind=Cind[!Cind %in% CGind]
              DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
           })
   )
}

#concerts all Cto T in a non CG or GCcontext
bisConvCGGC=function(in.seq){
   C.pos=vmatchPattern('C', in.seq)
   CG.pos=vmatchPattern('CG', in.seq)
   GC.pos=vmatchPattern('GC', in.seq)
   
   do.call(c,
           lapply(seq(length(in.seq)),function(i){
              Cind=start(C.pos[[i]])
              CGind=start(CG.pos[[i]])
              GCind=start(GC.pos[[i]])
              convind1=!Cind %in% CGind  
              convind2=! Cind %in% (GCind+1)
              convind= Cind[convind1 & convind2]
              DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
           })
   )
}

writeFasta<-function(mySeqs,mySeqsBS,mySeqsBSnotCG) {
   for(i in seq(length(mySeqs))) {
      seqName<-names(mySeqs)[i]
      writeXStringSet(mySeqs[i],paste0("./ampliconSeqs/",seqName,".fa"))
      writeXStringSet(mySeqsBS[i],paste0("./ampliconSeqs/",seqName,".fa"),append=TRUE)
      writeXStringSet(mySeqsBSnotCG[i],paste0("./ampliconSeqs/",seqName,".fa"),append=TRUE)
   }
}


mySeqsBS<-bisConv(mySeqs)
names(mySeqsBS)<-paste0("BS_",names(mySeqs))
mySeqsBSnotCG<-bisConvCG(mySeqs)
names(mySeqsBSnotCG)<-paste0("BSnotCG_",names(mySeqs))

writeFasta(mySeqs,mySeqsBS,mySeqsBSnotCG)


#########################################################################3
###### get amplicon seqs with gene strand
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
