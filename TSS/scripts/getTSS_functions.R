library(GenomicRanges)
library(rtracklayer)
#library(genomation)

###################
# functions to liftover genome coordinates for metadata columns 
constructGRforMcol<-function(myGRanges,myMcol) {
  ind<-which(!is.na(mcols(myGRanges)[,myMcol]))
  myIRanges<-IRanges(start=mcols(myGRanges)[ind,myMcol],width=1)
  tempGR<-GRanges(seqnames=seqnames(myGRanges)[ind],
                  ranges=myIRanges, strand=strand(myGRanges)[ind])
  return(tempGR)
}

doLiftOver<-function(myGRanges,chainFile="/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain") {
  myChain<-import.chain(chainFile)
  tempGR<-unlist(liftOver(myGRanges,myChain))
  return(tempGR)
}

liftOverMcols<-function(myGRanges,mcolList,chainFile="/SharedDocuments/MeisterLab/GenomeVer/chainFiles/ce10ToCe11.over.chain") {
  for (mcolName in mcolList) {
    ind<-which(!is.na(mcols(myGRanges)[,mcolName]))
    tempGR<-constructGRforMcol(myGRanges[ind,],mcolName)
    tempGR<-doLiftOver(tempGR,chainFile)
    mcols(myGRanges)[ind,mcolName]<-start(tempGR)
  }
  return(myGRanges)
}
#######################################